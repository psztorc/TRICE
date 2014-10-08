%--------------------------------------------------------------
%This function decomposes welfare effects into three components
%--------------------------------------------------------------
function [TOTEFFECT PROFEFFECT TRADEEFFECT]=mydecomposition(TARIFFCs,NXC,RHO)
mycalculations
COEFF1s=TRADEs./repmat(EXPENDITURE',[N 1 S]);
COEFF2s=PROFITSs./EXPENDITUREs;
COEFF3as=TARIFFs./repmat(EXPENDITURE',[N 1 S]);
COEFF3bs=TARIFFs.*TRADEs./repmat(EXPENDITURE',[N 1 S]);
[GOVERNMENTWELFAREHAT WELFAREHAT WAGEHAT TRADECs]=mycounterfactuals(TARIFFCs,NXC,RHO);
DELTAWAGE=WAGEHAT-1;
DELTAPROFITSs=sum(TRADECs,2)./sum(TRADEs,2)-1;
TOTEFFECT=100*sum(sum(COEFF1s.*(repmat(DELTAWAGE',[N 1 S])-repmat(DELTAWAGE,[1 N S])),1),3)';
PROFEFFECT=100*sum(COEFF2s.*(DELTAPROFITSs-repmat(DELTAWAGE,[1 1 S])),3);
TRADEEFFECT=100*sum(sum((COEFF3as.*(TRADECs-TRADEs)-COEFF3bs.*repmat(DELTAWAGE,[1,N,S])),1),3)';
end

%This is checked and correct