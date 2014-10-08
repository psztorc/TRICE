%--------------------------------------------------------
%This function computes the optimal tariff subject to MFN
%--------------------------------------------------------
function [MFNOPTIMALTARIFFj GOVERNMENTWELFAREHAT EXPENDITUREHAT WAGEHAT]=mymfnoptimaltariffj(j,TARIFFOTHERs,LAMBDA,LBj,UBj,MFNOPTIMALTARIFFGUESSj)
%Defining global variables
global ARG3 ARG4 ARG5 SCALING MFNTAUj %Defined here
%Preliminary calculations
mycalculations
ARG3=j-1;
ARG4=TARIFFOTHERs-1;
ARG5=LAMBDA-1;
%Calculating guess
if nargin>5
    GUESSOj=MFNOPTIMALTARIFFGUESSj;
elseif nargin<=5
    if LAMBDA(j,:)==ones(1,S)
        GUESSOj=1./(SIGMA-1); %This is motivated by the 2x1 optimal tariff formula
    elseif LAMBDA(j,:)~=ones(1,S)
        GUESSOj=0.1*ones(S,1); %Starting too low is best because it does not depress trade volumes
    end
end
TEMP=repmat(reshape(GUESSOj,[1 1 S]),[N-1,1,1]);
TARIFFCs=TARIFFOTHERs;
TARIFFCs(:,j,:)=[TEMP(1:j-1,:,:);zeros(1,1,S);TEMP(j:end,:,:)];
[GGUESS SKIP1 WGUESS SKIP2 SKIP3 XGUESS]=mycounterfactuals(TARIFFCs,zeros(N,1),LAMBDA);
%Transforming tariff to improve performance of fmincon
SCALING=0.1; %Brings order of magnitude of tariff step size in line with other step sizes
TARIFFj=TARIFFOTHERs(:,j,:);
TARIFFj(j,:,:)=[];
MFNTARIFFj=reshape(mean(TARIFFj,1),S,1);
MFNTAUj=1+SCALING*MFNTARIFFj;
GUESSOj=(1+SCALING*GUESSOj)./MFNTAUj; %Makes sure optimization runs over tariff change
LBO=(1+SCALING*repmat(LBj,S,1))./MFNTAUj;
if size(UBj,1)==1
    UBO=(1+SCALING*repmat(UBj,S,1))./MFNTAUj;
elseif size(UBj,1)>1
    UBO=(1+SCALING*UBj)./MFNTAUj;
end
GUESS=[GGUESS;XGUESS;WGUESS;GUESSOj]; %[GOVERNMENTWELFAREHAT;EXPENDITUREHAT;WAGEHAT;MFNOPTIMALTARIFFj(transformed)]
LB=[0.1*ones(N,1);0.1*ones(2*N,1);LBO];
UB=[10*ones(N,1);10*ones(2*N,1);UBO];
%Computing OPTIMALTARIFFj
tol=1e-13;
options=optimset('algorithm','active-set','display','off','TolFun',tol,'TolX',tol,'TolCon',tol,'MaxFunEvals',inf,'MaxIter',5000);
%options=optimset('algorithm','active-set','display','off','TolFun',tol,'TolX',tol,'TolCon',tol,'MaxFunEvals',inf,'MaxIter',5000,'plotfcns',{@optimplotx,@optimplotfunccount,@optimplotfval,@optimplotconstrviolation,@optimplotstepsize, @optimplotfirstorderopt});
x=fmincon(@myfunopt,GUESS,[],[],[],[],LB,UB,@mymfnconopt,options);
GOVERNMENTWELFAREHAT = x(1:N);
EXPENDITUREHAT=x(N+1:2*N);
WAGEHAT=x(2*N+1:3*N);
MFNOPTIMALTARIFFj=(x(3*N+1:end).*MFNTAUj-1)./SCALING;
end
        
%This is checked and correct