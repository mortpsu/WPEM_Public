function [sector_losses_ALL] =CalcSullivanLosses(IterationCount, UnservedData)


% STEP 0: Load Sullivan Parameters, coefficients, and definitions
% load('timeslice32.mat','timeslice32');
% load('Day_Season_Multiplier.mat','Day_Season_Multiplier');
% load ('CaliIndustryMultiplier.mat','CaliIndustryMultiplier');
% load('ROIndustryMultiplier.mat','ROIndustryMultiplier');
load('SullivanParameters.mat','timeslice32','ind_coeff_season','ind_coeff_region'); % residential seasonal (res_coeff_season) is in the .mat


% Restore this to run on ACI
%cd results
%filename = ['Iteration', int2str(IterationCount)];
%cd(filename)
%filename = ['UnservedEnergyData', int2str(IterationCount), '.csv'];

% Use a test file to debug - turn off for ACI
%%%filename = 'UnservedEnergyDataTEST.csv';


%%% STEP 1: Get hourly nonserved energy and demand data from last PSM run
%NSEFILEDATA = readmatrix(filename);

%%% USE INPUT ARGUMENT TO FUNCTION - DO NOT READ FILE
% define number of economic regions
num_econr = size(UnservedData,2)/2;  % == number of economic regions in PSM
num_drem_econr = size(ind_coeff_region,3); % == number of economic regions in DREM

medemand_dr  = UnservedData(:,1:num_econr);
nsedemand_dr = UnservedData(:,num_econr+1:end);

% pre-allocate output matrix arrary (S x econr - S == sectors - 1)
sector_losses_ALL = zeros(size(ind_coeff_region,1),num_drem_econr);

%%% loop
for i = 1:1:num_drem_econr
    
    %%% STEP 2: Go through hours chronologically, and look for sequences of
    %%% non-zeros

    % save avg percentage nse of each event in a list of events
    % store in a cell array: rows=sullivan timeslice, cols=# hrs in event
	[event_data]    = Count_Nonzero_Events(medemand_dr(:,i), nsedemand_dr(:,i), timeslice32);
    
    %%% STEP 3: Sum up the total losses by sector
    [sector_losses] = CalcSectoralLosses(event_data, ind_coeff_season, ind_coeff_region(:,:,i));
    
    %%% STEP 4: prepare and return final table of losses by sector and region
    %%% for IMPLAN
    sector_losses_ALL(:,i) = sector_losses;
    
end

%%% STEP 2: Go through hours chronologically, and look for sequences of
%%% non-zeros

% save avg percentage nse of each event in a list of events
% store in a cell array: rows=sullivan timeslice, cols=# hrs in event
% [event_data_CA] = Count_Nonzero_Events(demand_CA, nse_CA, timeslice32);
% [event_data_RO] = Count_Nonzero_Events(demand_RO, nse_RO, timeslice32);


% consolidate into only 1, 4, and 8 hour events
% [event_data_148_CA] = Consolidate_Sullivan_Events(event_data_CA);
% [event_data_148_RO] = Consolidate_Sullivan_Events(event_data_RO);


%%% STEP 3: Sum up the total losses by sector
%[sector_losses_CA] = CalcSectoralLosses(event_data_CA, Day_Season_Multiplier, CA_Ind_Coeffs);
%[sector_losses_RO] = CalcSectoralLosses(event_data_RO, Day_Season_Multiplier, RO_Ind_Coeffs);


%%% STEP 4: prepare and return final table of losses by sector and region
%%% for IMPLAN
%sector_losses_ALL = [sector_losses_CA sector_losses_RO zeros(length(sector_losses_CA),1)];


end
