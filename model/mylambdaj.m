%--------------------------------
%This function calculates LAMBDAj
%--------------------------------
function [LAMBDAj MFNOPTIMALTARIFFj]=mylambdaj(j,TARGETTARIFFj)
%Preliminary calculations
DATA=load('DATA');
TARIFF=DATA.TARIFF;
N=size(TARIFF,2);
S=size(TARIFF,1)./N;
TEMP=reshape(TARIFF',[N N S]);
TARIFFs=permute(TEMP,[2 1 3]);
%Imposing upper bounds for close to prohibitive tariffs
PROHIBj=(TARGETTARIFFj>2.25);
UB=max(TARGETTARIFFj+0.03,2.25); %Slightly larger bound than prohibitive TARGETTARIFFj to ensure that boundary LAMBDAj is computed
%Initializing guess
GUESSBEF=[TARGETTARIFFj,[1:S]',ones(S,1)];
MFNOPTIMALTARIFFGUESSj=0.1*ones(S,1);
LAMBDAPOS=ones(S,1);
TEST=-1;
ITER=0;
MAXITER=250;
PROGRESS=zeros(S,7,MAXITER);
%Calculating LAMBDAj
while (TEST==-1 && ITER<MAXITER)
    %Checking current guess
    LAMBDA=ones(N,S);
    LAMBDA(j,:)=GUESSBEF(:,3)';
    [MFNOPTIMALTARIFFj GOVERNMENTWELFAREHAT]=mymfnoptimaltariffj(j,TARIFFs,LAMBDA,0,UB,MFNOPTIMALTARIFFGUESSj);
    CHECK=sortrows([TARGETTARIFFj,[1:S]',MFNOPTIMALTARIFFj,PROHIBj,UB],[1 2]);
    PPRELIM=find(CHECK(:,4)==0,1,'last');
    SHOULDBEPRELIM=CHECK(:,1)-mean(CHECK(1:PPRELIM,1))+mean(CHECK(1:PPRELIM,3));
    P=min(find((SHOULDBEPRELIM>UB)==0,1,'last'),PPRELIM);
    SHOULDBE=max(min([CHECK(1:P,1)-mean(CHECK(1:P,1))+mean(CHECK(1:P,3));CHECK(P+1:S,5)-0.03],CHECK(:,5)-0.03),0);
    CHECK=[CHECK(:,1:3),SHOULDBE,CHECK(:,3)-SHOULDBE];
    CHECK(:,5)=LAMBDAPOS.*CHECK(:,5);
    %Calculating correction
    CLOSE=0;
    TEMP=zeros(S,6);
    for s=1:S
        for v=1:6
            CUTOFF=0.1/(2^(6-v));
            if max(abs(CHECK(:,5)))>CUTOFF;
                CLOSE=v;
                if CHECK(s,5)<-CUTOFF
                    TEMP(s,v)=1;
                elseif CHECK(s,5)>CUTOFF
                    TEMP(s,v)=-1;
                end
            end
        end
    end
    MULT=[0.00005;0.00005;0.0005;0.005;0.05;0.1];
    if ITER<75
        PICKMULT=max(CLOSE,3);
        PICKTEMP=CLOSE;
    elseif (ITER>=75 && ITER<125)
        PICKMULT=min(CLOSE,3);
        PICKTEMP=min(CLOSE,3);
    elseif (ITER>=125 && ITER<175)
        PICKMULT=min(CLOSE,2);
        PICKTEMP=min(CLOSE,2);
     elseif (ITER>=175 && ITER<=MAXITER)
        PICKMULT=min(CLOSE,4);
        PICKTEMP=CLOSE;
    end
    if CLOSE>0
        CORRECTION=MULT(PICKMULT).*TEMP(:,PICKTEMP);
    end
    CORRECTION=LAMBDAPOS.*CORRECTION;
    %Determining convergence
    if  (CLOSE==0 || (CLOSE==1 && ITER>100) || (CLOSE==2 && ITER>150) || (CLOSE==3 && ITER>200))
        TEST=1;
        CORRECTION=zeros(S,1);
    end
    %Applying correction
    GUESSAFT=sortrows(GUESSBEF,[1 2]);
    LAMBDAPOS=(GUESSAFT(:,3)>0);
    KEEP=GUESSAFT(:,3);
    GUESSAFT(:,3)=LAMBDAPOS.*((GUESSAFT(:,3)+CORRECTION)-mean(GUESSAFT(:,3)+CORRECTION)+1);
    GUESSAFT(:,3)=max(GUESSAFT(:,3),0)./mean(max(GUESSAFT(:,3),0));
    GUESSAFT=sortrows(GUESSAFT,2);
    %Showing progress and finishing up
    PROGRESS(:,:,ITER+1)=[CHECK,CORRECTION,KEEP];
    [CHECK,CORRECTION,KEEP]
    [ITER,CLOSE,GOVERNMENTWELFAREHAT(j)]
    GUESSBEF=GUESSAFT;
    if ITER<10
        MFNOPTIMALTARIFFGUESSj=MFNOPTIMALTARIFFj;
    elseif ITER>=10
        MFNOPTIMALTARIFFGUESSj=SHOULDBE;
        if j==5
            MFNOPTIMALTARIFFGUESSj=TARGETTARIFFj; %Works better for Japan
        end
    end
    ITER=ITER+1;
end
%Picking iteration that minimizes RSS (relevant only if earlier algorithm does not fully converge)
PROGRESS(:,:,ITER+1:MAXITER)=[];
MINRSS=sortrows([reshape(sum(PROGRESS(:,5,:).^2),ITER,1),[1:ITER]'],1);
MINPROGRESS=sortrows(PROGRESS(:,:,MINRSS(1,2)),2);
MFNOPTIMALTARIFFj=MINPROGRESS(:,3);
LAMBDAj=MINPROGRESS(:,7);
end

%This is checked and correct
