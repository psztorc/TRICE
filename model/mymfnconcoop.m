%-------------------------------------------------------
%These are the constraints from mymfncooperativetariff.m
%-------------------------------------------------------
function [c,ceq]=mymfnconcoop(x)
%Defining global variables
global N S %Defined in mycalculations
global ARG SCALING TAU %Defined in mycooperativetariff
%Preparing calculations of [c,ceq]
LAMBDA=ARG+1;
GOVERNMENTWELFAREHAT = x(1:N);
EXPENDITUREHAT=x(N+1:2*N);
WAGEHAT=x(2*N+1:3*N);
MFNCOOPERATIVETARIFFTRANSF=x(3*N+1:end);
TEMP=(reshape(MFNCOOPERATIVETARIFFTRANSF,S,N).*TAU-1)./SCALING;
TARIFFCs=zeros(N,N,S);
for j=1:N
    TEMP1=repmat(reshape(TEMP(:,j),[1 1 S]),[N-1 1 1]);
    TARIFFCs(:,j,:)=[TEMP1(1:j-1,1,:);zeros(1,1,S);TEMP1(j:end,1,:)];
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