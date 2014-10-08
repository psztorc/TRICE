%-------------------------------------------------------
%This function computes the change in government welfare
%-------------------------------------------------------
function [GOVERNMENTWELFAREHAT LOBBYWELFAREHAT]=mycongov(LAMBDA,TRADECs,TARIFFCs,AGGREGATEPRICEINDEXHAT)
%Defining global variables
global TRADEs TARIFFs SIGMAs_1_1_S N S DELTAs %Defined in mycalculations 
%Preparing calculation of GOVERNMENTWELFAREHAT
LAMBDAs=reshape(LAMBDA,[N 1 S]);
EXPENDITURELOBBYs=sum(TRADEs,2)+DELTAs.*repmat(sum(sum(TARIFFs.*TRADEs,1),3)',[1 1 S]);
EXPENDITUREADJ=sum(LAMBDAs.*EXPENDITURELOBBYs,3);
TRADEADJCs=repmat((SIGMAs_1_1_S-1)./SIGMAs_1_1_S,[N N 1]).*TRADECs;
DELTACs=sum(TRADEADJCs,2)./repmat(sum(sum(TRADEADJCs,3),2),[1 1 S]);
EXPENDITURELOBBYCs=sum(TRADECs,2)+DELTACs.*repmat(sum(sum(TARIFFCs.*TRADECs,1),3)',[1 1 S]);
EXPENDITURELOBBYHATs=EXPENDITURELOBBYCs./EXPENDITURELOBBYs;
EXPENDITUREADJC=sum(LAMBDAs.*EXPENDITURELOBBYCs,3);
EXPENDITUREADJHAT=EXPENDITUREADJC./EXPENDITUREADJ;
%Calculating GOVERNMENTWELFAREHAT and LOBBYWELFAREHAT
GOVERNMENTWELFAREHAT=EXPENDITUREADJHAT./AGGREGATEPRICEINDEXHAT;
LOBBYWELFAREHAT=reshape(EXPENDITURELOBBYHATs./repmat(AGGREGATEPRICEINDEXHAT,[1 1 S]),[N S 1]);
end

%This is checked and correct