
%----------------------
% OSSA Model
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

clear all
close all
clc


  Shocks = (0:10)/100;
% Shocks = (25:65)/100;
% Shocks = [0.5014162, 0.5330983, 0.5408727, 0.4320262, 0.5654997, 0.5451430, 0.5540263];

L = length(Shocks);
C = 7; %seven countries/regions


% 1 Brazil
% 2 China
% 3 EU
% 4 India
% 5 Japan
% 6 ROW
% 7 US



mycalculations

% create a place to dump the results
BigResults = zeros(C,L,35,(C-1));  %35 happens to be max { rows( nchoosek(1:7,1:7) ) }
BigWhoIsIn = ones(35,C,(C-1)) * -1; % later, we will remove all -1's;

% separately, keep all the tarrifs
dimTARIFFs = size(TARIFFs);
AllTariffs = zeros(dimTARIFFs(1),dimTARIFFs(2),dimTARIFFs(3),L);

%for ClubSize = 1 % special edit - for simplicity (when only considering
%club sizes of one
    
for ClubSize = 1:(C-1) %clubs can be all or one
    
    ClubsOfThisClubSize = nchoosek(1:C,ClubSize);
    dimClubs = size(ClubsOfThisClubSize); %rows and columns of these clubs 
    qClubs = dimClubs(1);
    
    WhoIsIn = zeros(qClubs,C); % declare for use - all out 
    
    for ClubIndex = 1:qClubs
        Club = ClubsOfThisClubSize(ClubIndex,:);
 
        % Keep track of who is in
        for n = 1:C % ...for each "i" trading partner...
            if ismember(n, Club)
                WhoIsIn(ClubIndex,n) = 1;
            end
        end
        
        BigWhoIsIn(ClubIndex,:,ClubSize) = WhoIsIn(ClubIndex,:);
        
        for Shock = 1:L   %for each shock...

            %have the importer shock the other countries by increasing their
            %tarrifs by Shocks(n)
            TARIFFCs=TARIFFs;  %reset tarrifs
            
            for n = 1:C % ...for each "i" trading partner...        
                if not(ismember(n, Club)) %If they aren't in the climate club, they're getting taxed (rows)
                    for o = 1:C %...for each importer...     
                        if ismember(o, Club) %only club-members are doing the taxing (columns)
                            % TARIFFCs(n,o,:)=TARIFFCs(n,o,:)+Shocks(Shock);
                            % % additive method
                            TARIFFCs(n,o,:)= (1+TARIFFCs(n,o,:))*(1+Shocks(Shock))-1;

                        end
                    end
                end
            end
    

        %welfare calculations
        NXC=zeros(N,1); %NXC are counterfactual aggregate net exports. I use this to purge the raw data from aggregate trade imbalances as described in the main text
        LAMBDA=LAMBDABAS; %Select LAMBDABAS if you don't want the lobbying weights, and LAMBDAPOL otherwise
        [GOVERNMENTWELFAREHAT,WELFAREHAT,WAGEHAT,TRADECs,LOBBYWELFAREHAT,EXPENDITUREHAT]=mycounterfactuals(TARIFFCs,NXC,LAMBDA);
        BigResults(:,Shock,ClubIndex,ClubSize) = WELFAREHAT;  %country-welfare-rows, tarrif-effects (columns), by The Club Itself, ClubSize 

        % save also the tarrifs
        AllTariffs(:,:,:,Shock) = TARIFFCs;
        
        end    
    end     
    

    
end

save('results\bigresults.mat', 'BigResults') ;

save('results\shocks.mat', 'Shocks') ;

save('results\basetariffs.mat', 'TARIFFs') ;

save('results\alltariffs.mat', 'AllTariffs') ;
    
save('results\whoisin.mat', 'BigWhoIsIn') ;
    
