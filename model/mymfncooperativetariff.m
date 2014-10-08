%------------------------------------------
%This function computes cooperative tariffs
%------------------------------------------
function [MFNCOOPERATIVETARIFFs GOVERNMENTWELFAREHAT EXPENDITUREHAT WAGEHAT ceq]=mymfncooperativetariff(LAMBDA,LBj,UBj,LIBONLY,MFNCOOPERATIVETARIFFGUESSs)
%Defining global variables
global ARG SCALING TAU %Defined here
%Preliminary calculations
mycalculations
ARG=LAMBDA-1;
%Calculating guess
if nargin>4
    MFNCOOPERATIVEIMPTARIFFGUESS=zeros(S,N);
    for j=1:N
        MFNCOOPERATIVEIMPTARIFFGUESS(:,j,:)=reshape(mean(MFNCOOPERATIVETARIFFGUESSs([1:j-1 j+1:N],j,:),1),S,1);
    end
elseif nargin<=4
    if isequal(LAMBDA,ones(N,S))==1
        INVMARKUP=(SIGMA-1)./SIGMA; %This is motivated by the desire to eliminate markup distortions
        TEMP=0.25*(INVMARKUP./mean(INVMARKUP)-1);
        MFNCOOPERATIVEIMPTARIFFGUESS=repmat(TEMP,1,N);
        TEMP=repmat(reshape(TEMP,[1 1 S]),[N-1,N,1]);
        MFNCOOPERATIVETARIFFGUESSs=zeros(N,N,S);
        for j=1:N
            MFNCOOPERATIVETARIFFGUESSs(:,j,:)=[TEMP(1:j-1,j,:);zeros(1,1,S);TEMP(j:end,j,:)];
        end
    elseif isequal(LAMBDA,ones(N,S))==0
        MFNCOOPERATIVEIMPTARIFFGUESS=zeros(S,N);
        MFNCOOPERATIVETARIFFGUESSs=zeros(N,N,S);
    end
end
[GGUESS SKIP1 WGUESS SKIP2 SKIP3 XGUESS]=mycounterfactuals(MFNCOOPERATIVETARIFFGUESSs,zeros(N,1),LAMBDA);
%Transforming tariff to improve performance of fmincon
SCALING=0.1;
MFNIMPTARIFF=zeros(S,N);
for j=1:N
    MFNIMPTARIFF(:,j)=reshape(mean(TARIFFs([1:j-1 j+1:N],j,:),1),S,1);
end
TAU=1+SCALING*MFNIMPTARIFF;%Brings order of magnitude of tariff step size in line with other step sizes
CGUESSTRANSF=reshape((1+SCALING*MFNCOOPERATIVEIMPTARIFFGUESS)./TAU,N*S,1); %Makes sure optimization runs over tariff changes
LBCTRANSF=reshape((1+SCALING*LBj)./TAU,N*S,1);
if LIBONLY==0
    UBCTRANSF=reshape((1+SCALING*UBj)./TAU,N*S,1);
elseif LIBONLY==1
    UBCTRANSF=ones(N*S,1); %Only tariff reductions are allowed
end
GUESS=[GGUESS;XGUESS;WGUESS;CGUESSTRANSF]; %[GOVERNMENTWELFAREHAT;EXPENDITUREHAT;WAGEHAT;MFNCOOPERATIVETARIFF(transformed)]
LB=[0.1*ones(N,1);0.1*ones(2*N,1);LBCTRANSF];
UB=[10*ones(N,1);10*ones(2*N,1);UBCTRANSF];
tol=1e-10;
%Computing MFNCOOPERATIVETARIFFs
options=optimset('algorithm','active-set','display','off','TolFun',tol,'TolX',tol,'TolCon',tol,'MaxFunEvals',inf,'MaxIter',5000);
%options=optimset('algorithm','active-set','TolFun',tol,'TolX',tol,'TolCon',tol,'MaxFunEvals',inf,'MaxIter',5000,'plotfcns',{@optimplotx,@optimplotfunccount,@optimplotfval,@optimplotconstrviolation,@optimplotstepsize, @optimplotfirstorderopt});
x=fmincon(@myfuncoop,GUESS,[],[],[],[],LB,UB,@mymfnconcoop,options);
[c,ceq]=mymfnconcoop(x);
GOVERNMENTWELFAREHAT = x(1:N);
EXPENDITUREHAT=x(N+1:2*N);
WAGEHAT=x(2*N+1:3*N);
MFNCOOPERATIVETARIFFTRANSF=x(3*N+1:end);
TEMP=(reshape(MFNCOOPERATIVETARIFFTRANSF,S,N).*TAU-1)./SCALING;
MFNCOOPERATIVETARIFFs=zeros(N,N,S);
for j=1:N
    TEMP1=repmat(reshape(TEMP(:,j),[1 1 S]),[N-1 1 1]);
    MFNCOOPERATIVETARIFFs(:,j,:)=[TEMP1(1:j-1,1,:);zeros(1,1,S);TEMP1(j:end,1,:)];
end
if LAMBDA==ones(N,S)
    if LBj==0
        RESTRICTEDMFNCOOPERATIVETARIFFBASs=MFNCOOPERATIVETARIFFs;
        save('RESTRICTEDMFNCOOPERATIVETARIFFBASs','RESTRICTEDMFNCOOPERATIVETARIFFBASs')
    elseif LBj<0
        UNRESTRICTEDMFNCOOPERATIVETARIFFBASs=MFNCOOPERATIVETARIFFs;
        save('UNRESTRICTEDMFNCOOPERATIVETARIFFBASs','UNRESTRICTEDMFNCOOPERATIVETARIFFBASs')
    end
elseif LAMBDA~=ones(N,S)
    if LBj==0
        RESTRICTEDMFNCOOPERATIVETARIFFPOLs=MFNCOOPERATIVETARIFFs;
        save('RESTRICTEDMFNCOOPERATIVETARIFFPOLs','RESTRICTEDMFNCOOPERATIVETARIFFPOLs')
    elseif LBj<0
        UNRESTRICTEDMFNCOOPERATIVETARIFFPOLs=MFNCOOPERATIVETARIFFs;
        save('UNRESTRICTEDMFNCOOPERATIVETARIFFPOLs','UNRESTRICTEDMFNCOOPERATIVETARIFFPOLs')
    end
end
end
        
%This is checked and correct