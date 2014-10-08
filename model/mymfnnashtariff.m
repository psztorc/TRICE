%------------------------------------------------------
%This function computes the Nash tariffs subject to MFN
%------------------------------------------------------
function MFNNASHTARIFFs=mymfnnashtariff(LAMBDA,LBj,UBj,MFNNASHTARIFFGUESSs)
%Preliminary calculations
N=size(MFNNASHTARIFFGUESSs,1);
S=size(MFNNASHTARIFFGUESSs,3);
TARIFFOTHERBEFs=MFNNASHTARIFFGUESSs;
MFNIMPTARIFFOTHERBEF=zeros(S,N);
for j=1:N
    TEMP=MFNNASHTARIFFGUESSs(:,j,:);
    TEMP(j,:,:)=[];
    MFNIMPTARIFFOTHERBEF(:,j)=reshape(mean(TEMP,1),S,1);
end
%Choosing iteration parameters and initiating
maxit=24;
t=1;
tol=1e-2;
error=2*tol;
%Iterating over optimal tariffs
while ((t<=maxit)&&(error>tol))
    TARIFFOTHERAFTs=TARIFFOTHERBEFs;
    MFNIMPTARIFFOTHERAFT=MFNIMPTARIFFOTHERBEF;
    matlabpool open
    parfor j=1:N
        MFNOPTIMALTARIFFj=mymfnoptimaltariffj(j,TARIFFOTHERBEFs,LAMBDA,LBj,UBj,MFNIMPTARIFFOTHERBEF(:,j));
        TEMP=repmat(reshape(MFNOPTIMALTARIFFj,[1,1,S]),[N-1,1,1]);
        TARIFFOTHERAFTs(:,j,:)=[TEMP(1:j-1,:,:);zeros(1,1,S);TEMP(j:end,:,:)];
        MFNIMPTARIFFOTHERAFT(:,j)=MFNOPTIMALTARIFFj;
    end
    matlabpool close
    error=abs(mean(mean(mean(TARIFFOTHERAFTs-TARIFFOTHERBEFs))));
    TARIFFOTHERBEFs=TARIFFOTHERAFTs;
    MFNIMPTARIFFOTHERBEF=MFNIMPTARIFFOTHERAFT;
    t=t+1;
end
MFNNASHTARIFFs=TARIFFOTHERAFTs;
end

%This is checked and correct