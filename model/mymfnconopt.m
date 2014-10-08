%-------------------------------------------------
%These are the constraints from mymfnoptimaltariffj.m
%-------------------------------------------------
function [c,ceq]=mymfnconopt(x)
%Defining global variables
global N S %Defined in mycalculations 
global ARG3 ARG4 ARG5 SCALING MFNTAUj %Defined in mymfnoptimaltariffj
%Preparing calculation of [c,ceq]
j=ARG3+1;
TARIFFOTHERs=ARG4+1;
LAMBDA=ARG5+1;
GOVERNMENTWELFAREHAT = x(1:N);
EXPENDITUREHAT=x(N+1:2*N);
WAGEHAT=x(2*N+1:3*N);
MFNOPTIMALTARIFFj=(x(3*N+1:end).*MFNTAUj-1)./SCALING;
TEMP=repmat(reshape(MFNOPTIMALTARIFFj,[1 1 S]),[N-1,1,1]);
TARIFFCs=TARIFFOTHERs;
TARIFFCs(:,j,:)=[TEMP(1:j-1,:,:);zeros(1,1,S);TEMP(j:end,:,:)];
[EQUILIBRIUMCONDITION NUMERAIRECONDITION PROFITSHATs AGGREGATEPRICEINDEXHAT TRADECs]=myconditions(EXPENDITUREHAT,WAGEHAT,TARIFFCs,zeros(N,1));
EQUILIBRIUMCONDITION(1)=[];
GOVERNMENTWELFAREHATCON=mycongov(LAMBDA,TRADECs,TARIFFCs,AGGREGATEPRICEINDEXHAT);
ceq1=[EQUILIBRIUMCONDITION;NUMERAIRECONDITION];
ceq2=GOVERNMENTWELFAREHAT-GOVERNMENTWELFAREHATCON;
%Calculating [c,ceq]
c=[];
ceq=[ceq1;ceq2];
end

%This is checked and correct