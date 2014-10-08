%----------------------------------------------------
%This program performs frequently needed calculations
%----------------------------------------------------
%Defining global variables
global TRADEs TARIFFs N S ALPHAs GAMMAs DELTAs SIGMAs_1_1_S SIGMAs_N_1_S SIGMAs_N_N_S MUs EXPENDITURE PROFITSs EXPENDITUREs LABORINCOME_div_EXPENDITURE NX_div_EXPENDITURE PROFITSs_div_EXPENDITUREs LABORINCOME %Defined here
%Loading data
load DATA
%Performing calculations
N=size(TRADE,2);
S=size(TRADE,1)./N;
TEMP=reshape(TRADE',[N N S]);
TRADEs=permute(TEMP,[2 1 3]);
TEMP=reshape(TARIFF',[N N S]);
TARIFFs=permute(TEMP,[2 1 3]);
WEIGHTs=TRADEs.*(1-repmat(eye(N,N),[1 1 S]))./repmat(sum(TRADEs.*(1-repmat(eye(N,N),[1 1 S])),1),[N 1 1]);
MEANTARIFF=permute(sum(WEIGHTs.*TARIFFs),[3 2 1]);
SIGMAs_1_1_S=reshape(SIGMA,[1 1 S]);
SIGMAs_N_1_S=repmat(SIGMAs_1_1_S,[N 1 1]);
SIGMAs_N_N_S=repmat(SIGMAs_N_1_S,[1 N 1]);
ALPHA=TRADE./repmat(sum(TRADE,2),1,N);
TEMP=reshape(ALPHA',[N N S]);
ALPHAs=permute(TEMP,[2 1 3]);
TAUTRADEs=(1+TARIFFs).*TRADEs;
GAMMAs=TAUTRADEs./repmat(sum(TAUTRADEs,1),[N 1 1]);
SIGMAs=reshape(SIGMA,[1 1 S]);
TRADEADJs=repmat((SIGMAs-1)./SIGMAs,[N N 1]).*TRADEs;
DELTAs=sum(TRADEADJs,2)./repmat(sum(sum(TRADEADJs,3),2),[1 1 S]);
TAUTRADEs=(1+TARIFFs).*TRADEs;
MUs=permute(sum(TAUTRADEs,1),[2 1 3])./repmat(sum(sum(TAUTRADEs,1),3)',[1 1 S]);
EXPENDITURE=sum(TRADE.*(1+TARIFF),1)';
EXPENDITUREs=repmat(EXPENDITURE,[1 1 S]);
PROFITSs=(1./SIGMAs_N_1_S).*sum(TRADEs,2);
NX=sum(sum(TRADEs,2)-permute(sum(TRADEs,1),[2 1 3]),3);
TARIFFREVENUE=sum(TARIFF.*TRADE,1)';
LABORINCOME=EXPENDITURE+NX-TARIFFREVENUE-sum(PROFITSs,3);
LABORINCOME_div_EXPENDITURE=LABORINCOME./EXPENDITURE;
NX_div_EXPENDITURE=NX./EXPENDITURE;
PROFITSs_div_EXPENDITUREs=PROFITSs./EXPENDITUREs;

%This is checked and correct

