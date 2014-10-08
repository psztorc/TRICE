%----------------------------------------------------
%These are the constraints from mycooperativetariff.m
%----------------------------------------------------
function [c,ceq]=myconcoop(x)
%Defining global variables
global N S %Defined in mycalculations
global ARG SCALING TAUs %Defined in mycooperativetariff
%Preparing calculations of [c,ceq]
LAMBDA=ARG+1;
GOVERNMENTWELFAREHAT = x(1:N);
EXPENDITUREHAT=x(N+1:2*N);
WAGEHAT=x(2*N+1:3*N);
COOPERATIVETARIFFTRANSF=x(3*N+1:end);
TEMP=(reshape(COOPERATIVETARIFFTRANSF,[N-1,N,S]).*TAUs-1)./SCALING;
TARIFFCs=zeros(N,N,S);
for n=1:N
    TARIFFCs(:,n,:)=[TEMP(1:n-1,n,:);zeros(1,1,S);TEMP(n:end,n,:)];
end
[EQUILIBRIUMCONDITION NUMERAIRECONDITION PROFITSHATs AGGREGATEPRICEINDEXHAT TRADECs]=myconditions(EXPENDITUREHAT,WAGEHAT,TARIFFCs,zeros(N,1));
EQUILIBRIUMCONDITION(1)=[];
GOVERNMENTWELFAREHATCON=mycongov(LAMBDA,TRADECs,TARIFFCs,AGGREGATEPRICEINDEXHAT);
ceq1=[EQUILIBRIUMCONDITION;NUMERAIRECONDITION];
ceq2=GOVERNMENTWELFAREHAT-GOVERNMENTWELFAREHATCON;
ceq3=[ones(N-1,1) -eye(N-1)]*GOVERNMENTWELFAREHAT;
%Calculating [c,ceq]
c=[];
ceq=[ceq1;ceq2;ceq3];
end

%This is checked and correct