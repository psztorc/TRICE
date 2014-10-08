function MINUSGOVERNMENTWELFAREHATj=myfunopt(x)
%Defining global variables
global ARG3 %this is defined in myoptimaltariffj
%Preparing calculation of MINUSGOVERNMENTWELFAREHATj
j=ARG3+1;
%Calculating MINUSGOVERNMENTWELFAREHATj
MINUSGOVERNMENTWELFAREHATj=-x(j);
end

%This is checked and correct