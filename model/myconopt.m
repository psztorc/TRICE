%-------------------------------------------------
%These are the constraints from myoptimaltariffj.m
%-------------------------------------------------
function [c,ceq]=myconopt(x)
%Defining global variables
global N S %Defined in mycalculations 
global ARG3 ARG4 ARG5 SCALING TAUj %Defined in myoptimaltariffj
%Preparing calculation of [c,ceq]
j=ARG3+1;
TARIFFOTHERs=ARG4+1;
LAMBDA=ARG5+1;
GOVERNMENTWELFAREHAT = x(1:N);
EXPENDITUREHAT=x(N+1:2*N);
WAGEHAT=x(2*N+1:3*N);
OPTIMALTARIFFj=(x(3*N+1:end).*TAUj-1)./SCALING;
TEMP=reshape(OPTIMALTARIFFj,[N-1,1,S]);
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