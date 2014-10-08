%-----------------------------------------
%This function computes the optimal tariff
%-----------------------------------------
function [OPTIMALTARIFFj GOVERNMENTWELFAREHAT EXPENDITUREHAT WAGEHAT]=myoptimaltariffj(j,TARIFFOTHERs,LAMBDA,LBj,UBj,OPTIMALTARIFFGUESSj)
%Defining global variables
global ARG3 ARG4 ARG5 SCALING TAUj %Defined here
%Preliminary calculations
mycalculations
ARG3=j-1;
ARG4=TARIFFOTHERs-1;
ARG5=LAMBDA-1;
%Calculating guess
if nargin>5
    OGUESS=OPTIMALTARIFFGUESSj;
elseif nargin<=5
    if LAMBDA(j,:)==ones(1,S)
        OGUESS=1./(reshape(repmat(SIGMAs,[N-1,1,1]),[(N-1)*S,1,1])-1); %This is motivated by the 2x1 optimal tariff formula
    elseif LAMBDA(j,:)~=ones(1,S)
        OGUESS=0.1*ones((N-1)*S,1); %Starting too low is best because it does not depress trade volumes
    end
end
TEMP=reshape(OGUESS,[N-1,1,S]);
TARIFFCs=TARIFFOTHERs;
TARIFFCs(:,j,:)=[TEMP(1:j-1,:,:);zeros(1,1,S);TEMP(j:end,:,:)];
[GGUESS SKIP1 WGUESS SKIP2 SKIP3 XGUESS]=mycounterfactuals(TARIFFCs,zeros(N,1),LAMBDA);
%Transforming tariff to improve performance of fmincon
SCALING=0.1;
IMPTARIFFj=reshape(TARIFFOTHERs([1:j-1 j+1:N],j,:),(N-1)*S,1);
TAUj=1+SCALING*IMPTARIFFj; %Brings order of magnitude of tariff step size in line with other step sizes
OGUESSTRANSF=(1+SCALING*OGUESS)./TAUj; %Makes sure optimization runs over tariff change
LBOTRANSF=(1+SCALING*repmat(LBj,(N-1)*S,1))./TAUj;
UBOTRANSF=(1+SCALING*repmat(UBj,(N-1)*S,1))./TAUj;
GUESS=[GGUESS;XGUESS;WGUESS;OGUESSTRANSF]; %[GOVERNMENTWELFAREHAT;EXPENDITUREHAT;WAGEHAT;OPTIMALTARIFFj(transformed)]
LB=[0.1*ones(N,1);0.1*ones(2*N,1);LBOTRANSF];
UB=[10*ones(N,1);10*ones(2*N,1);UBOTRANSF];
%Computing OPTIMALTARIFFj
tol=1e-13;
options=optimset('algorithm','active-set','display','off','TolFun',tol,'TolX',tol,'TolCon',tol,'MaxFunEvals',inf,'MaxIter',5000);
%options=optimset('algorithm','active-set','display','off','TolFun',tol,'TolX',tol,'TolCon',tol,'MaxFunEvals',inf,'MaxIter',5000,'plotfcns',{@optimplotx,@optimplotfunccount,@optimplotfval,@optimplotconstrviolation,@optimplotstepsize, @optimplotfirstorderopt});
x=fmincon(@myfunopt,GUESS,[],[],[],[],LB,UB,@myconopt,options);
GOVERNMENTWELFAREHAT = x(1:N);
EXPENDITUREHAT=x(N+1:2*N);
WAGEHAT=x(2*N+1:3*N);
OPTIMALTARIFFj=(x(3*N+1:end).*TAUj-1)./SCALING;
end
        
%This is checked and correct