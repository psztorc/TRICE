%----------------------------------------------------------
%----------------------------------------------------------
%This program illustrates a few things you might want to do
%----------------------------------------------------------
%----------------------------------------------------------

%----------------------
%Understanding the data
%----------------------

clear all
close all
clc

mycalculations %This program performs a number of frequently needed calculations
TRADEs %TRADEs(i,j,s) is the value of trade flowing from country i to country j in industry s; see regions.csv and sectors.csv for the list of regions and industries
TARIFFs %TARIFFs(i,j,s) is the tariff applied by country j against industry s imports from country i
SIGMA %These are the elasticities of substitution from Table 1 of the paper
NX %These are aggregate net exports which are basically zero indicating that this is the purged version of the raw data (as explained in the paper)

%-------------------------------------------------
%Computing the effects of exogenous tariff changes
%-------------------------------------------------

%China increases its tariffs against the US by 20 percentage points
clear all
close all
clc

i=7; %This is the US
j=2; %This is China

mycalculations
TARIFFCs=TARIFFs;
TARIFFCs(i,j,:)=TARIFFCs(i,j,:)+0.2; %These are now the counterfactual tariffs
NXC=zeros(N,1); %NXC are counterfactual aggregate net exports. I use this to purge the raw data from aggregate trade imbalances as described in the main text
LAMBDA=LAMBDABAS; %Select LAMBDABAS if you don't want the lobbying weights, and LAMBDAPOL otherwise
[GOVERNMENTWELFAREHAT,WELFAREHAT,WAGEHAT,TRADECs,LOBBYWELFAREHAT,EXPENDITUREHAT]=mycounterfactuals(TARIFFCs,NXC,LAMBDA);
100*(WELFAREHAT-1); %These are the percentage changes in welfare

%-------------------------
%Computing optimal tariffs
%-------------------------

%The US imposes optimal tariffs
clear all
close all
clc

j=7; %This is the US

mycalculations
TARIFFOTHERs=TARIFFs; %Here you select the tariffs the US is optimally responding to (the ones pertaining to the US will be replaced). Here, I have just chosen factual tariffs
LAMBDA=LAMBDABAS; %Select LAMBDABAS if you don't want the lobbying weights, and LAMBDAPOL otherwise
LBj=0; %Here you set the lower bound of tariffs to limit iterations
UBj=10; %Here you set the upper bound of tariffs to limit iterations

%Now you can either compute optimal tariffs under the restriction of MFN or
%unconstrained optimal tariffs. As you can see in the paper, they are very
%similar so just working with MFN tariffs might be OK. In any case, I
%suggest computing MFN tariffs first and then using them as guess for the
%unconstrained ones to improve accuracy, as I illustrate below.

[MFNOPTIMALTARIFFj GOVERNMENTWELFAREHAT EXPENDITUREHAT WAGEHAT]=mymfnoptimaltariffj(j,TARIFFOTHERs,LAMBDA,LBj,UBj); %this takes 18 seconds on my computer
100*MFNOPTIMALTARIFFj %This is the US optimal MFN tariff by industry in percent
100*(GOVERNMENTWELFAREHAT-1) %This is the percentage change in government welfare (and welfare if LAMBDA=LAMBDABAS)

OPTIMALTARIFFGUESSj=reshape(repmat(MFNOPTIMALTARIFFj',N-1,1),(N-1)*S,1); %This expands MFNOPTIMALTARIFFj so that it can be used as a guess for the unconstrained case (you will understand the format from the output below)
[OPTIMALTARIFFj,GOVERNMENTWELFAREHAT,EXPENDITUREHAT,WAGEHAT]=myoptimaltariffj(j,TARIFFOTHERs,LAMBDA,LBj,UBj,OPTIMALTARIFFGUESSj); %This takes 322 second on my computer
100*OPTIMALTARIFFj %This is the optimal tariff in percent ordered by country and industry with the US always left out
100*(GOVERNMENTWELFAREHAT-1) %This is the percentage change in government welfare (and welfare if LAMBDA=LAMBDABAS)
TEMP=reshape(OPTIMALTARIFFj,[N-1,1,S]); %Now it is split up by industry
TARIFFCs=TARIFFOTHERs;
TARIFFCs(:,j,:)=[TEMP(1:j-1,:,:);zeros(1,1,S);TEMP(j:end,:,:)]; %Now you have the tariff matrix with US optimal tariffs for use, say, in mycounterfcatuals()
[GOVERNMENTWELFAREHAT,WELFAREHAT,WAGEHAT,TRADECs,LOBBYWELFAREHAT,EXPENDITUREHAT]=mycounterfactuals(TARIFFCs,zeros(N,1),LAMBDA);
100*(WELFAREHAT-1) %These are the percentage changes in welfare
