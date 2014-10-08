%------------------------------------------
%This function computes cooperative tariffs
%------------------------------------------
function [COOPERATIVETARIFFs GOVERNMENTWELFAREHAT EXPENDITUREHAT WAGEHAT ceq]=mycooperativetariff(LAMBDA,LBj,UBj,LIBONLY,COOPERATIVETARIFFGUESSs)
%Defining global variables
global ARG SCALING TAUs %Defined here
%Preliminary calculations
mycalculations
ARG=LAMBDA-1;
%Calculating guess
if nargin>4
    COOPERATIVEIMPTARIFFGUESSs=zeros(N-1,N,S);
    for j=1:N
        COOPERATIVEIMPTARIFFGUESSs(:,j,:)=COOPERATIVETARIFFGUESSs([1:j-1 j+1:N],j,:);
    end
elseif nargin<=4
    if isequal(LAMBDA,ones(N,S))==1
        INVMARKUP=(SIGMA-1)./SIGMA; %This is motivated by the desire to eliminate markup distortions
        TEMP=0.25*(INVMARKUP./mean(INVMARKUP)-1);
        COOPERATIVEIMPTARIFFGUESSs=repmat(reshape(TEMP,[1 1 S]),[N-1,N,1]);
        COOPERATIVETARIFFGUESSs=zeros(N,N,S);
        for j=1:N
            COOPERATIVETARIFFGUESSs(:,j,:)=[COOPERATIVEIMPTARIFFGUESSs(1:j-1,j,:);zeros(1,1,S);COOPERATIVEIMPTARIFFGUESSs(j:end,j,:)];
        end
    elseif isequal(LAMBDA,ones(N,S))==0
        COOPERATIVEIMPTARIFFGUESSs=zeros(N-1,N,S);
        COOPERATIVETARIFFGUESSs=zeros(N,N,S);
    end 
end
[GGUESS SKIP1 WGUESS SKIP2 SKIP3 XGUESS]=mycounterfactuals(COOPERATIVETARIFFGUESSs,zeros(N,1),LAMBDA);
%Transforming tariff to improve performance of fmincon
SCALING=0.1;
IMPTARIFFs=zeros(N-1,N,S);
for j=1:N
    IMPTARIFFs(:,j,:)=TARIFFs([1:j-1 j+1:N],j,:);
end
TAUs=1+SCALING*IMPTARIFFs;%Brings order of magnitude of tariff step size in line with other step sizes
CGUESSTRANSF=reshape(reshape((1+SCALING*COOPERATIVEIMPTARIFFGUESSs)./TAUs,[(N-1)*N,1,S]),[(N-1)*N*S,1]); %Makes sure optimization runs over tariff changes
LBCTRANSF=reshape(reshape((1+SCALING*LBj)./TAUs,[(N-1)*N,1,S]),[(N-1)*N*S,1]);
if LIBONLY==0
    UBCTRANSF=reshape(reshape((1+SCALING*UBj)./TAUs,[(N-1)*N,1,S]),[(N-1)*N*S,1]);
elseif LIBONLY==1
    UBCTRANSF=ones((N-1)*N*S,1); %Only tariff reductions are allowed
end
GUESS=[GGUESS;XGUESS;WGUESS;CGUESSTRANSF]; %[GOVERNMENTWELFAREHAT;EXPENDITUREHAT;WAGEHAT;COOPERATIVETARIFF(transformed)]
LB=[0.1*ones(N,1);0.1*ones(2*N,1);LBCTRANSF];
UB=[10*ones(N,1);10*ones(2*N,1);UBCTRANSF];
tol=1e-10;
%Computing COOPERATIVETARIFF
options=optimset('algorithm','active-set','display','off','TolFun',tol,'TolX',tol,'TolCon',tol,'MaxFunEvals',inf,'MaxIter',5000);
%options=optimset('algorithm','active-set','TolFun',tol,'TolX',tol,'TolCon',tol,'MaxFunEvals',5000,'MaxIter',inf,'plotfcns',{@optimplotx,@optimplotfunccount,@optimplotfval,@optimplotconstrviolation,@optimplotstepsize, @optimplotfirstorderopt});
x=fmincon(@myfuncoop,GUESS,[],[],[],[],LB,UB,@myconcoop,options);
[c,ceq]=myconcoop(x);
GOVERNMENTWELFAREHAT = x(1:N);
EXPENDITUREHAT=x(N+1:2*N);
WAGEHAT=x(2*N+1:3*N);
COOPERATIVETARIFFTRANSF=x(3*N+1:end);
TEMP=(reshape(COOPERATIVETARIFFTRANSF,[N-1,N,S]).*TAUs-1)./SCALING;
COOPERATIVETARIFFs=zeros(N,N,S);
for n=1:N
    COOPERATIVETARIFFs(:,n,:)=[TEMP(1:n-1,n,:);zeros(1,1,S);TEMP(n:end,n,:)];
end
if LAMBDA==ones(N,S)
    if LBj==0
        RESTRICTEDCOOPERATIVETARIFFBASs=COOPERATIVETARIFFs;
        save('RESTRICTEDCOOPERATIVETARIFFBASs','RESTRICTEDCOOPERATIVETARIFFBASs')
    elseif LBj<0
        UNRESTRICTEDCOOPERATIVETARIFFBASs=COOPERATIVETARIFFs;
        save('UNRESTRICTEDCOOPERATIVETARIFFBASs','UNRESTRICTEDCOOPERATIVETARIFFBASs')
    end
elseif LAMBDA~=ones(N,S)
    if LBj==0
        RESTRICTEDCOOPERATIVETARIFFPOLs=COOPERATIVETARIFFs;
        save('RESTRICTEDCOOPERATIVETARIFFPOLs','RESTRICTEDCOOPERATIVETARIFFPOLs')
    elseif LBj<0
        UNRESTRICTEDCOOPERATIVETARIFFPOLs=COOPERATIVETARIFFs;
        save('UNRESTRICTEDCOOPERATIVETARIFFPOLs','UNRESTRICTEDCOOPERATIVETARIFFPOLs')
    end
end
end
        
%This is checked and correct