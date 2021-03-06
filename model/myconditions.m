%--------------------------------------------------
%This function compiles important model conditions
%-------------------------------------------------
function [EQUILIBRIUMCONDITION NUMERAIRECONDITION PROFITSHATs AGGREGATEPRICEINDEXHAT TRADECs]=myconditions(EXPENDITUREHAT,WAGEHAT,TARIFFCs,NXC)
%Defining global variables
global TRADEs TARIFFs N S ALPHAs GAMMAs DELTAs SIGMAs_N_1_S SIGMAs_N_N_S MUs EXPENDITURE EXPENDITUREs LABORINCOME_div_EXPENDITURE PROFITSs_div_EXPENDITUREs LABORINCOME %Defined in mycalculations
%Miscellaneous calculations
TAUHATs=(1+TARIFFCs)./(1+TARIFFs);
WAGEHATs_1_1_S=repmat(WAGEHAT,[1 1 S]);
WAGEHATs_1_N_S=repmat(WAGEHAT,[1 N S]);
EXPENDITUREHATs=repmat(EXPENDITUREHAT,[1 1 S]);
%Computing PRICEINDEXHATADJs=PRICEINDEXHATs.^(SIGMAs-1)
PRICEINDEXHATADJs=(permute(sum(GAMMAs.*((WAGEHATs_1_N_S.*TAUHATs).^(1-SIGMAs_N_N_S)),1),[2 1 3])).^(-1);
%Computing PROFITSHATs
TEMP1=(repmat(WAGEHATs_1_1_S,[1 N 1]).^(1-SIGMAs_N_N_S)).*(TAUHATs.^(-SIGMAs_N_N_S)).*(repmat(permute(PRICEINDEXHATADJs,[2 1 3]),[N 1 1]));
TEMP2=TEMP1.*repmat(permute(EXPENDITUREHATs,[2 1 3]),[N 1 1]);
PROFITSHATs=sum(ALPHAs.*TEMP2,2);
%Computing EQUILIBRIUMCONDITION1: WAGEHAT
EQUILIBRIUMCONDITION1=WAGEHAT-sum(DELTAs.*PROFITSHATs,3);
%Computing EQUILIBRIUMCONDITION2: EXPENDITUREHAT
EQUILIBRIUMCONDITION2=EXPENDITUREHAT-(LABORINCOME_div_EXPENDITURE.*WAGEHAT-NXC./EXPENDITURE+sum(PROFITSs_div_EXPENDITUREs.*PROFITSHATs,3))./(1-sum(sum(((TRADEs./repmat(permute(EXPENDITUREs,[2 1 3]),[N 1 1])).*TARIFFCs).*TEMP1,1),3)');
%Computing EQUILIBRIUMCONDITION
EQUILIBRIUMCONDITION=[EQUILIBRIUMCONDITION1;EQUILIBRIUMCONDITION2];
%Calculating AGGREGATEPRICEINDEXHAT
PRICEINDEXHATs=PRICEINDEXHATADJs.^(1./(SIGMAs_N_1_S-1));
AGGREGATEPRICEINDEXHAT=prod(PRICEINDEXHATs.^MUs,3);
%Calculating NUMERAIRECONDITION
NUMERAIRECONDITION=sum((LABORINCOME./repmat(sum(LABORINCOME),N,1)).*WAGEHAT)-1;
%Calculating TRADECs
TRADECs=TRADEs.*TEMP2;
end

%This is checked and correct