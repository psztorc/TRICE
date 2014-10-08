%--------------------------------------------------------------
%This function solves for the effects of counterfactual tariffs
%--------------------------------------------------------------
function [GOVERNMENTWELFAREHAT WELFAREHAT WAGEHAT TRADECs LOBBYWELFAREHAT EXPENDITUREHAT]=mycounterfactuals(TARIFFCs,NXC,LAMBDA)
%Defining global variables
global N %Defined in mycalculations
global ARG1 ARG2 %Defined here 
%Calculating parameters
ARG1=TARIFFCs-1;
ARG2=NXC-1;
%Solving for equilibrium
y0=ones(2*N,1);
options=optimset('display','off');
y=fsolve('mysolvedconditions',y0,options);
%Retrieving results and computing WELFAREHAT
EXPENDITUREHAT=y(1:N);
WAGEHAT=y(N+1:2*N);
[EQUILIBRIUMCONDITION NUMERAIRECONDITION PROFITSHATs AGGREGATEPRICEINDEXHAT TRADECs]=myconditions(EXPENDITUREHAT,WAGEHAT,TARIFFCs,NXC);
[GOVERNMENTWELFAREHAT LOBBYWELFAREHAT]=mycongov(LAMBDA,TRADECs,TARIFFCs,AGGREGATEPRICEINDEXHAT);
WELFAREHAT=EXPENDITUREHAT./AGGREGATEPRICEINDEXHAT;
end

%This is checked and correct