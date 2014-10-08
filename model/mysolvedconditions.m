%-------------------------------------------------------------
%This function picks the conditions for which MATLAB "fsolves"
%-------------------------------------------------------------
function SOLVEDCONDITIONS=mysolvedconditions(x)
%Defining global variables
global N %Defined in mycalculations
global ARG1 ARG2 %Defined in mycounterfactual
%Miscellaneous calculations
TARIFFCs=ARG1+1;
NXC=ARG2+1;
%Defining fsolved variables
EXPENDITUREHAT=x(1:N);
WAGEHAT=x(N+1:2*N);
%Compiling fsolved conditions
[EQUILIBRIUMCONDITION NUMERAIRECONDITION]=myconditions(EXPENDITUREHAT,WAGEHAT,TARIFFCs,NXC);
EQUILIBRIUMCONDITION(1)=[]; %It does not matter which condition is dropped
SOLVEDCONDITIONS=[EQUILIBRIUMCONDITION;NUMERAIRECONDITION];
end

%This is checked and correct