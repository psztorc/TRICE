%---------------------------------------
%This function computes the Nash tariffs
%---------------------------------------
function NASHTARIFFs=mynashtariff(LAMBDA,LBj,UBj,NASHTARIFFGUESSs)
%Preliminary calculations
N=size(NASHTARIFFGUESSs,1);
S=size(NASHTARIFFGUESSs,3);
TARIFFOTHERBEFs=NASHTARIFFGUESSs;
OPTIMALTARIFFGUESS=zeros((N-1)*S,N);
for j=1:N
    TEMP=NASHTARIFFGUESSs(:,j,:);
    TEMP(j,:,:)=[];
    OPTIMALTARIFFGUESS(:,j)=reshape(TEMP,(N-1)*S,1);
end
%Choosing iteration parameters and initiating
maxit=24;
t=1;
tol=1e-2;
error=2*tol;
%Iterating over optimal tariffs
while ((t<=maxit)&&(error>tol))
    TARIFFOTHERAFTs=TARIFFOTHERBEFs;
    matlabpool open
    parfor j=1:N
        OPTIMALTARIFFj=myoptimaltariffj(j,TARIFFOTHERBEFs,LAMBDA,LBj,UBj,OPTIMALTARIFFGUESS(:,j));
        TEMP=reshape(OPTIMALTARIFFj,[N-1,1,S]);
        TARIFFOTHERAFTs(:,j,:)=[TEMP(1:j-1,:,:);zeros(1,1,S);TEMP(j:end,:,:)];
    end
    matlabpool close
    error=abs(mean(mean(mean(TARIFFOTHERAFTs-TARIFFOTHERBEFs))));
    TARIFFOTHERBEFs=TARIFFOTHERAFTs;
    t=t+1;
end
NASHTARIFFs=TARIFFOTHERAFTs;
end

%This is checked and correct