%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write a script to read in the output file and then run the Economic 
% Model and then the Convergence script.
% 
% I.   The output from the job arrays will come here.
% II.  Then we will be passing a GDX file to Economic Model, and 
% III. then we will pass through the the Economic Model and PSM model 
%      to the convergence scripts.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z=Trial_Coupling(IterationCount)

	% --------------------------------------------------------------------
    % define working directory for iteration
    % --------------------------------------------------------------------
	
	% 		define working directory name: ./IterationX
	filename = ['Iteration',int2str(IterationCount)];
	
	% 		move to psm detailed results directory: ./codedir/Results_Each_Iter
	cd Results_Each_Iter
	
	% working directory: ./codedir/Results_Each_Iter/IterationX
	cd(filename)
	
	% --------------------------------------------------------------------
    % export csv of iteration count
    % --------------------------------------------------------------------
	
	filename = ['Count', int2str(IterationCount), '.csv']; 		
	csvwrite(filename,IterationCount);                       
	
	% --------------------------------------------------------------------
    % import main results csv
    % --------------------------------------------------------------------
	
	% 		loop through each week and import results	
	for i=1:1:52
		
		%		 import csv file
		disp(i);							     % display week i	
		filename = ['week', int2str(i), '.csv']; % define filename   
		z1       = csvread(filename);            % read csv file
		
		% 		define vectors
		NonServedEnergy(i)    = z1(1);           % non-served electricity - total - week
		TotalCostCali(i)      = z1(2);           % total cost of electricity generation in CA - total - week
		TotalCostRO(i)        = z1(3);           % total cost of electricity generation in RO - total - week
		GenCali(i)            = z1(4);           % electricity generation in CA - total - week
		GenRO(i)              = z1(5);           % electricity generation in RO - total - week
		CaliCoal(i)           = z1(6);           % coal electricity generation in CA - total - week
		CaliGas(i)            = z1(7);           % natural gas CT, CCGT, steam, turbine electricity generation in CA  - total - week
		ROCoal(i)             = z1(8);           % coal electricity generation in ROWECC - total - week
		ROGas(i)              = z1(9);           % natural gas CT, CCGT, steam, turbine electricity generation in ROWECC  - total - week
		FuelUsageCoalCali(i)  = z1(32);          % coal fuel usage CA - total - week
		FuelUsageGasCali(i)   = z1(33);          % natural gas (CT and CCGT) fuel usage CA - total - week
		FuelUsageCoalRO(i)    = z1(34);          % coal fuel usage ROWEE - total - week
		FuelUsageGasRO(i)     = z1(35);          % natural gas (CT and CCGT) fuel usage ROWECC - total - week
		UnservedEnergyCali(i) = z1(36);          % non-served electricity in CA - total - week  - added 5
		UnservedEnergyRO(i)   = z1(37);          % non-served electricity in RO - total - week  - added 5
		
	end
	
	% --------------------------------------------------------------------
    % define annual values for costs, generation, and fuel usage
    % --------------------------------------------------------------------

	% 		total annual cost
	YearlyCostCali = sum(TotalCostCali);
	YearlyCostRO   = sum(TotalCostRO);
	
	%		total annual electricity generation
	YearlyGenCali  = sum(GenCali);
	YearlyGenRO    = sum(GenRO);
	
	% 		total annual electricity generation for fuel types coal and 
	% 		natural gas (CT and CCGT)
	YearlyCaliCoal = sum(CaliCoal);
	YearlyCaliGas  = sum(CaliGas);
	YearlyROCoal   = sum(ROCoal);
	YearlyROGas    = sum(ROGas);

	% 		total annual fuel usage for fuel types coal and 
	%		natural gas (CT and CCGT)
	YearlyFuelCoalCali = sum(FuelUsageCoalCali);
	YearlyFuelGasCali  = sum(FuelUsageGasCali);
	YearlyFuelCoalRO   = sum(FuelUsageCoalRO);
	YearlyFuelGasRO    = sum(FuelUsageGasRO);
	
	%		total annual electricity price
	CaliAvgElectricPrice = YearlyCostCali/YearlyGenCali;
	ROAvgElectricPrice   = YearlyCostRO/YearlyGenRO;

	% --------------------------------------------------------------------
    % import second set of main results csv files
    % --------------------------------------------------------------------

	% 		loop through each week and import results	
	for week =1:1:52
	
		% 		csv file for marginal cost in CA of bus i for hour d
		filename = ['MarginalCostCali', int2str(week), '.csv'];
		MarginalCostCali{week} = csvread(filename);
		
		% 		csv file for marginal cost in ROWECC of bus i for hour d
		filename = ['MarginalCostRO', int2str(week), '.csv'];
		MarginalCostRO{week} = csvread(filename);
		
		% 		csv file for electricity demand in CA of bus i for hour d
		filename = ['MultipliedDemandCali', int2str(week), '.csv'];
		MultipliedDemandCali{week} = csvread(filename);
		
		% 		csv file for electricity demand in ROWECC of bus i for hour d
		filename = ['MultipliedDemandRO', int2str(week), '.csv'];
		MultipliedDemandRO{week} = csvread(filename);
		
		
		% 		csv file for non-served electricity at bus i for hour d
		filename = ['UnservedEnergywhole', int2str(week), '.csv'];
		UnservedEnergy{week} = csvread(filename);
		
		% 		csv file for electricity generation in CA for hour d
		filename = ['CaliGenperhour', int2str(week), '.csv'];
		CaliforniaGenperhour{week} = csvread(filename);
		
		% 		csv file for electricity generation in CA for hour d
		filename = ['ROGenperhour', int2str(week), '.csv'];
		RoweccGenperhour{week} = csvread(filename);
		
		% 		csv file for shed energy at bus i for hour d
		filename = ['ShedEnergywhole', int2str(week), '.csv'];
		ShedEnergy{week} = csvread(filename);
			
		% 		csv file for electricity demand in CA - total for hour d	
		filename = ['NewDemandwholeCali', int2str(week), '.csv'];
		z11{week} = csvread(filename);
		
		% 		csv file for electricity demand in ROWECC - total for hour d
		filename = ['NewDemandwholeRO', int2str(week), '.csv'];
		z12{week} = csvread(filename);

		% 		csv file for non-served electricity in CA for hour d
		filename = ['NewEnergywholeunservedCali', int2str(week), '.csv'];
		z13{week} = csvread(filename);
		
		% 		csv file for non-served electricity in ROWECC for hour d
		filename = ['NewEnergywholeunservedRO', int2str(week), '.csv'];
		z14{week}=csvread(filename);     

	end

	% --------------------------------------------------------------------
    % return to main coupled model directory
    % --------------------------------------------------------------------
	
	% 		move to psm detailed results directory: ./codedir/Results_Each_Iter
	cd .. 
	
	% 		move to main coupled model directory: ./codedir
	cd ..

	
	% --------------------------------------------------------------------
    % calculate electricity productivity parameter for the regional 
	% economic model (REM)
    % --------------------------------------------------------------------

	% 		convert import data from cell to matrix format
	MarginalCostCaliYear    = cell2mat(MarginalCostCali);     % MC(i,d) \in CA
	MarginalCostROYear      = cell2mat(MarginalCostRO);       % MC(i,d) \in RO
	UnservedEnergywholeyear = cell2mat(UnservedEnergy);       % NSE(i,d)
	ShedEnergywholeyear     = cell2mat(ShedEnergy);           % SE(i,d)
	CaliforniaGenwholeyear  = cell2mat(CaliforniaGenperhour); % GEN(d)  \in CA
	RoweccGenwholeyear      = cell2mat(RoweccGenperhour);     % GEN(d)  \in CA

	
	%		define index for columns with NSE
	linearIndexes2 = find(ShedEnergywholeyear>0);       % define the linear indices
	[~, columns2]  = ind2sub(size(ShedEnergywholeyear), linearIndexes2); %convert the linear indices to matrix notation 
																	     % but keep only the column notations

	%		remove columns with NSE from marginal cost and generation
	MarginalCostCaliYear(:,columns2)    =[];
	MarginalCostROYear(:,columns2)      =[];
	CaliforniaGenwholeyear(:,columns2)  =[];
	RoweccGenwholeyear(:,columns2)      =[];
	UnservedEnergywholeyear(:,columns2) =[];
	
	% --------------------------------------------------------------------
    % import crosswalk for buses to economic regions
    % --------------------------------------------------------------------

	State_Bus_Cali = csvread('Bus_State_Cali.csv'); % col1 == bus id' col2 == economic region id
	State_Bus_RO   = csvread('Bus_State_RO.csv');   % col1 == bus id' col2 == economic region id

	% --------------------------------------------------------------------
	%		define the rows correspond to California or Rowecc and then 
	%		change the prices of onlt those rows to the maximum -- WHY?
	% --------------------------------------------------------------------
	
	%		define index for rows and columns with NSE
	linearIndexes   = find(UnservedEnergywholeyear>0);
	[rows, columns] = ind2sub(size(UnservedEnergywholeyear), linearIndexes);
	
	%		define index for buses in CA
	[M,~] = ismember((rows'),State_Bus_Cali(:,1));
	
	%		define rows and correspond columns
	RowsCali_Inter = rows(M);
	ColsCali       = columns(M);
	
	%		define list of buses in CA -- ?
	[~,RowsCali]=ismember(RowsCali_Inter,State_Bus_Cali(:,1));

	%		define index for buses in ROWECC
	[M,~]=ismember((rows'),State_Bus_RO(:,1));
	
	%		define rows and correspond columns
	RowsRO_Inter = rows(M);
	ColsRO       = columns(M);
	[~,RowsRO]   = ismember(RowsRO_Inter,State_Bus_RO(:,1));

	%		define the maximum marinal price (MP)?
	MaxCalMP = 268;
	MaxROMP  = 268;

	%		replace columns with maximum MP for marginal cost and generation
	MarginalCostCaliYear(RowsCali,ColsCali) = MaxCalMP;
	MarginalCostROYear(RowsRO,ColsRO)       = MaxROMP;


	% --------------------------------------------------------------------
    % calculate electricity productivity parameter for the regional 
	% economic model (REM)
    % --------------------------------------------------------------------
	
	% 		define current reference prices for calculating percentage 
	% 		change from the baseline
	RefMarginalCostCali = 36.218;
	RefMarginalCostRO   = 34.535;

% 	% 		define current reference prices for calculating percentage 
%	% 		change from the baseline demand response (DR) specification
%   RefMarginalCostCali = 37.78;
%   RefMarginalCostRO   = 35.32;


% 	% 		define MORT/VIJAY's reference prices for calculating percentage 
%	% 		change from the baseline demand response (DR) specification
%   RefMarginalCostCali = 36.1970;
%   RefMarginalCostRO   = 34.5088;

% 	% 		define OLDER reference prices for calculating percentage 
%	% 		change from the baseline demand response (DR) specification
%   RefMarginalCostCali = 39.7470;
%   RefMarginalCostRO   = 39.1234;

	%		calculate the mean of marginal cost for each hour: MC_mean(d)
	MarginalPriceCalimean = sum(MarginalCostCaliYear,1)/44; 
	MarginalPriceROmean   = sum(MarginalCostROYear,1)/268;  
	% numerator is sum of mc over buses denominator is the number of bus in the economic region
	
	
	%		calculate total value of generation by hour sum(d, MC_mean(d)*GEN_(d))
	ValueCali = sum(MarginalPriceCalimean.*CaliforniaGenwholeyear); 
	ValueRO   = sum(MarginalPriceROmean.*RoweccGenwholeyear);

	%		calculate electricity generation weighted average marginal cost
	AvgMarginalCostCali = ValueCali/sum(CaliforniaGenwholeyear);
	AvgMarginalCostRO   = ValueRO/sum(RoweccGenwholeyear);


	%		define baseline reference price                           
	PrevIterMarginalCostCali = RefMarginalCostCali;
	PrevIterMarginalCostRO   = RefMarginalCostRO;
	
	%		calcuate percentage change in electricity generation weighted 
	%       average marginal cost from the baseline    
	%		(electricity productivity impact for REM)
	ChangeCali = (AvgMarginalCostCali-PrevIterMarginalCostCali)/PrevIterMarginalCostCali;
	ChangeRO   = (AvgMarginalCostRO-PrevIterMarginalCostRO)/PrevIterMarginalCostRO;
	
	% --------------------------------------------------------------------
    % calcuate non-electricity (Sullivan) productivity impacts from NSE
	% for REM
    % --------------------------------------------------------------------

	% 		convert import demand and NSE from cell to matrix format
	DemandCali   = cell2mat(z11); % bus x hour
	DemandRO     = cell2mat(z12); % bus x hour
	UnservedCali = cell2mat(z13); % bus x hour
	UnservedRO   = cell2mat(z14); % bus x hour

	%		call script: reshape function bus x hour matrix to single
	%		             vector --?
	DemandCali   = reshape(DemandCali,168*52,1);
	DemandRO     = reshape(DemandRO,168*52,1);
	UnservedCali = reshape(UnservedCali,168*52,1);
	UnservedRO   = reshape(UnservedRO,168*52,1);

	%		combine demand and NSE into a single matrix
	UnservedData = [DemandCali DemandRO UnservedCali UnservedRO];
	
	% 		calculate total annual NSE
	UnservedTotalCA = sum(UnservedCali);
	UnservedTotalRO = sum(UnservedRO);	
	disp([UnservedTotalCA UnservedTotalRO]);

	%		call script: to lable hourly date by hour, date, time of day, 
	%                    and seasonally. script creates a csv file--?
	output = Datetimestamp(UnservedData,IterationCount);

    % 		call script: non-electricity (Sullivan) productivity impacts 
	%                    from NSE
 	[sector_losses_ALL] =CalcSullivanLosses(IterationCount, UnservedData);
	
	% --------------------------------------------------------------------
    % create gams gdx files for REM: MtoCGE.gdx and MtoCGE2.gdx
	%
	% (MtoCGE == electricity productivity impacts == ; 
	%  MtoCGE2 == non-electricity productivity impacts)                          
    % --------------------------------------------------------------------	

	%		change working directory to the REM directory
	cd modelcge

	%		define structure of parameter 
	r.name    = 'r';
	r.uels    = {'CA','ROWECC','ROUS'};
	nele.name = 'nele';
	nele.uels = {'osa', 'grn', 'vna', 'oca', 'pfb', 'cba', 'apa', 'frs', 'con', 'fin', 'fbm', 'bom', 'wpm', 'per', 'cpm', 'cem', 'pfm', 'tec', 'trm', 'cru', 'coa', 'min', 'pub', 'hlt', 'bos', 'ngd', 'tel', 'rtl', 'trn'};
	% nele.uels = {'AGR', 'MIN', 'CONST', 'MANUF', 'TEL_UTL', 'TRD_RTL', 'FIN', 'SRV', 'PUB'}; % for static REM (SREM)

	% define parameter for non-electricity productivity impacts
	Sull_impacts.name = 'Sull_impacts';
	Sull_impacts.type = 'parameter';
	Sull_impacts.val  = sector_losses_ALL;
	Sull_impacts.form = 'full';
	Sull_impacts.uels = {nele.uels,r.uels};
	
	% export non-electricity productivity impacts gdx
	wgdx ('MtoCGE', Sull_impacts);
	
	% export electricity productivity impacts gdx
	iwgdx('MtoCGE2','ChangeCali','ChangeRO');

	% --------------------------------------------------------------------
    % call REM model which in the directory
    % --------------------------------------------------------------------

	%		user defined REM data targets
	% --------------------------------------------------------------------
	
	% define IMPLAN targets
	drem_sector = 30;
	drem_region = 12;         % 10,12,19,25
	drem_region_name = "";  % NULL or "nerc" if the user wants the nerc state regions and drem_region == 19 

	% define maximum time step target
	drem_simtp = 2050;        % 2010, 2100

	% define SSP target
	drem_ssp = 2;             % Switch to select SSP: == 1 SSP1, == 2 SSP2, == 3 SSP3 ,== 4 SSP4,== 5 SSP5,

%	% define switches
%	drem_tfp = 0;            % Switch to turn on tfp calibration: == 1 on; == 0 off

	% define couple model switche
	switch_coupled = 3;       % 0 == NONE , 
							  % 1 == WDRG coupling, 
							  % 2 == PSM electricity impacts (w/o Sullivan) coupling, 
							  % 3 == PSM w/Sullivan impacts coupling, 
							  % 4 == Demand response emulator (DRM) coupling, 
							  % 5 == PSM w/Sullivan impacts and DRM coupling

	% define demand response GCM
	drem_gcm = "GFDL-CM3";
	  % There are 21 GCMs  
		% ACCESS1-0 , bcc-csm1-1    , BNU-ESM      , CanESM2     , CCSM4     ,
		% CESM1-BGC , CNRM-CM5      , CSIRO-Mk3-6-0, GFDL-CM3    , GFDL-ESM2G,
		% GFDL-ESM2M, inmcm4        , IPSL-CM5A-LR , IPSL-CM5A-MR, MIROC5    ,
		% MIROC-ESM , MIROC-ESM-CHEM, MPI-ESM-LR   , MPI-ESM-MR  , MRI-CGCM3,
		% NorESM1-M
		
	% determine naming scheme for log and listing files
	drem_mod        = strjoin(["inter","_",int2str(IterationCount)],""); % user defines change first element		
	
	% 		define drem IMPLAN target name	(no user intervention)
	% --------------------------------------------------------------------

	% based on user defined variables
	implan_data     = strjoin(["implan440-",int2str(drem_sector),"sector-",int2str(drem_region),'region',drem_region_name],"");

	% create drem gams call string
	% --------------------------------------------------------------------
	
	% determine drem switch target for data to call
	target_implan   = strjoin(["--target=",implan_data],"");       %IMPLAN data
	target_region   = strjoin(["--t_r=",int2str(drem_region)],""); %region number  
	target_sector   = strjoin(["--t_s=",int2str(drem_sector)],""); %sector number
	target_lastyear = strjoin(["--t_y=",int2str(drem_simtp)],"");  %termination year 
	target_ssp      = strjoin(["--t_p=",int2str(drem_ssp)],"");    %spp number
	
	% determine	type of coupling drem has the with psm
	if (switch_coupled == 0)                       % No coupled model
	  target_couple = "";
	elseif (switch_coupled == 1)                  % GCAM-Land Coupled Model 
	  target_couple = strjoin(["--s_wd=",int2str(1)],"");  
	elseif (switch_coupled == 2)                  % SM electricity impacts (w/o Sullivan)
	  target_couple = strjoin(["--s_pw=",int2str(1)],""); 
	elseif (switch_coupled == 3)                  % PSM w/Sullivan impacts coupling
	  target_couple = strjoin([["--s_pw=",int2str(1)]," ",["--s_pf=",int2str(1)]],""); 
	elseif (switch_coupled == 4)                  % Demand response emulator (DRM) coupling
	  target_couple = strjoin(["--s_dr=",int2str(1)]," "); 
	elseif (switch_coupled == 5)                  % PSM w/Sullivan impacts and DRM coupling
	  target_couple = strjoin([["--s_pw=",int2str(1)]," ",["--s_pf=",int2str(1)]," ",["--s_dr=",int2str(1)]],"");  
	end

	% determine gcm and mod target
	target_gcm = strjoin(["--t_g=",drem_gcm],"");
	target_mod = strjoin(["--mod=",drem_mod],"");

	% determin listing and log output names based on mod
	target_lst = strjoin(["o=listings/",drem_mod,".lst"],"");
	target_log = strjoin(["lf=logs/",drem_mod,".log lo=2"],"");

	% create gams call string based on user defined variables from above
	gams_call = strjoin(["gams",'"pches_drem"',target_implan,target_region,target_sector,target_lastyear,target_ssp,target_couple,target_gcm,target_mod,target_lst,target_log]," ");

	%		call REM model from gams call string
	% ------------------------------------------------------------------
	system(gams_call);
%	system 'gams "pches_drem" --target=implan440-30sector-12region --t_r=12 --t_s=30 --t_y=2050 --t_p=2 --s_pw=1 --s_pf=1 --s_dr=1 --t_g=GFDL-CM3 --mod=psm lo = 2';

	% call SREM -- old commands
%	system 'gams "soe_mge_11sector_Sullivan_UC" lo = 3';
%	system 'gams "soe_mge_11sector_Sullivan_UC_noSull" lo = 3';


	% --------------------------------------------------------------------
    % import REM results from gdx (CCEtoM.gdx)
    % --------------------------------------------------------------------

	%now we read in the data to be read for solutions.
	%The data is demand, electricity price and armington electricity price

	% percentage change from the baseline:
	%	sectoral production for the domestic market 
	%	region CA for sector electricity at time period 2050
	dummygen.name      = 'd_elyCA';
	dummygen.form      = 'sparse';
	CaliDemandChange   = rgdx('CGEtoM',dummygen);
	clear dummygen;
	z(21) = CaliDemandChange.val;

	% percentage change from the baseline:
	%	sectoral production for the domestic market 
	%	region RO for sector electricity at time period 2050
	dummygen.name      = 'd_elyROWECC';
	dummygen.form      = 'sparse';
	RODemandChange     = rgdx('CGEtoM',dummygen);
	clear dummygen;
	z(22) = RODemandChange.val;

	% percentage change from the baseline:
	%	sectoral production price for the domestic market 
	%	region CA for sector electricity at time period 2050
	dummygen.name      = 'p_elyCA';
	dummygen.form      = 'sparse';
	PriceChangeCali    = rgdx('CGEtoM',dummygen);
	clear dummygen;
	z(23) = PriceChangeCali.val;

	% percentage change from the baseline:
	%	sectoral production price for the domestic market 
	%	region RO for sector electricity at time period 2050
	dummygen.name      = 'p_elyROWECC';
	dummygen.form      = 'sparse';
	PriceChangeRO      = rgdx('CGEtoM',dummygen);
	clear dummygen;
	z(24) = PriceChangeRO.val;

	% percentage change from the baseline:
	%	armington aggregate price 
	%	region CA for sector electricity at time period 2050
	dummygen.name      = 'pa_elyCA';
	dummygen.form      = 'sparse';
	ArmPriceChangeCali = rgdx('CGEtoM',dummygen);
	clear dummygen;
	z(25) = ArmPriceChangeCali.val;

	% percentage change from the baseline:
	%	armington aggregate price 
	%	region RO for sector electricity at time period 2050
	dummygen.name      = 'pa_elyROWECC';
	dummygen.form      = 'sparse';
	ArmPriceChangeRO   = rgdx('CGEtoM',dummygen);
	clear dummygen;
	z(26) = ArmPriceChangeRO.val;

	% --------------------------------------------------------------------
    % change working directory to the psm directory
    % --------------------------------------------------------------------
	
	% 		move to main coupled model directory: ./codedir
	cd ..
	
	% 		move to main coupled model directory: ./codedir/Gams_data
	cd Gams_data

	% --------------------------------------------------------------------
    % define stepsize logic for convergence
    % --------------------------------------------------------------------

	%		define current change in demand from the REM
	DemandChangeCali = (CaliDemandChange.val);            %percentage change in electricity demand (from rem) for CA
	DemandChangeRO   = (RODemandChange.val);              %percentage change in electricity demand (from rem) for RO
	DemandChange     = [DemandChangeCali DemandChangeRO]; %combine percentage changes
	
	%		display current percentage change in electricity demand (from rem)
	disp('demand change from implan');
	disp(DemandChange);
	
	%		define current iteration number
	kkkk = IterationCount;
	
	
	if (kkkk == 1)											% iteration = 1
	   
	   % First time with feedback, usually too large -- explain?
	   DemandChangeCGE = DemandChange;
	   DemandStep      = DemandChange;
	   DemandChange(DemandChange < -0.01) = -0.01;
	   DemandChangeADJ = DemandChange;
	
	else											        % iteration > 1
	   
	   % 		step 1: load previous iteration percentage change in electricity demand (from rem)
	   load('DemandChangeData.mat', 'DemandChangeCGE', 'DemandChangeADJ');
	   DemandChangeCGE = [DemandChangeCGE; DemandChange];

	   % 		step 2: check if converged is reached
		if ((DemandChangeCGE(kkkk-1,1) == DemandChange(1)) && (DemandChangeCGE(kkkk-1,2) == DemandChange(2))) 
	   
		  % DO NOT ADJUST -> CONVERGED           :END
	   
		else % ITERATION >= 2 AND NOT YET CONVERGED
		 
			%		step2a -- CA: define the difference between iterations 
			DemandStep = DemandChangeADJ(kkkk-1,:) - DemandChange;
		 
			if (abs(DemandStep(1)) < 0.001)         % do region 1
			
				% DO NOT ADJUST -> VERY SMALL CHANGE FROM LAST ITERATION
			
			elseif (abs(DemandStep(1)) <= 0.01) % small difference - split the difference			
				DemandStep(1) = (DemandStep(1))/2;			
			elseif ( DemandStep(1) >  0.3)      % CGE wants a bigger demand reduction, take a step	
				DemandStep(1) = 0.05;		
			elseif ( DemandStep(1) >  0.2)      % CGE wants a bigger demand reduction, take a step
				DemandStep(1) = 0.03;	
			elseif ( DemandStep(1) >  0.1)      % CGE wants a bigger demand reduction, take a step
				DemandStep(1) = 0.02;
			elseif ( DemandStep(1) >  0.01)     % CGE wants a bigger demand reduction, take a step
				DemandStep(1) = 0.01;		
			elseif ( DemandStep(1) <  -0.01)    % CGE wants a bigger demand increase, take a step	
				DemandStep(1) = -0.01;
			else                                % moderate change take a small step
				DemandStep(1) = DemandStep(1) / 4;
			end
			
			%		step2b -- CA: define the new pchg in demand based on step
			DemandChange(1) = DemandChangeADJ(kkkk-1,1) - DemandStep(1);


			% 		step2c -- CA: check if stuck cycling -- explain?
			if     ((kkkk > 4) && (DemandChange(1) == DemandChangeADJ(kkkk-2,1)) && (DemandChangeADJ(kkkk-1,1) == DemandChangeADJ(kkkk-3,1)))
				DemandChange(1) = (DemandChangeADJ(kkkk-2,1) + DemandChangeADJ(kkkk-1,1))/2;
			elseif ((kkkk > 10) && (DemandChange(1) < DemandChangeADJ(kkkk-1,1)) && (DemandChangeADJ(kkkk-1,1) > DemandChangeADJ(kkkk-2,1)) && (DemandChangeADJ(kkkk-2,1) < DemandChangeADJ(kkkk-3,1)))
				DemandChange(1) = (DemandChangeADJ(kkkk-2,1) + DemandChangeADJ(kkkk-1,1))/2;
			elseif ((kkkk > 10) && (DemandChange(1) > DemandChangeADJ(kkkk-1,1)) && (DemandChangeADJ(kkkk-1,1) < DemandChangeADJ(kkkk-2,1)) && (DemandChangeADJ(kkkk-2,1) > DemandChangeADJ(kkkk-3,1)))
				DemandChange(1) = (DemandChangeADJ(kkkk-2,1) + DemandChangeADJ(kkkk-1,1))/2;
			end


			%		step2a -- RO: define the difference between iterations 
			if (abs(DemandStep(2)) < 0.001)         % do region 2
			
				% DO NOT ADJUST -> VERY SMALL CHANGE FROM LAST ITERATION
			
			elseif (abs(DemandStep(2)) <= 0.01) % small difference - split the difference
				DemandStep(2) = (DemandStep(2))/2;
			elseif ( DemandStep(2) >  0.3)      % CGE wants a bigger demand reduction, take a step
				DemandStep(2) = 0.05;
			elseif ( DemandStep(2) >  0.2)		% CGE wants a bigger demand reduction, take a step
				DemandStep(2) = 0.03;
			elseif ( DemandStep(2) >  0.1)		% CGE wants a bigger demand reduction, take a step
				DemandStep(2) = 0.02;
			elseif ( DemandStep(2) >  0.01)     % CGE wants a bigger demand reduction, take a step
				DemandStep(2) = 0.01;
			elseif ( DemandStep(2) <  -0.01)    % CGE wants a bigger demand increase, take a step
				DemandStep(2) = -0.01;
			else                                % moderate change take a small step
				DemandStep(2) = DemandStep(2) / 4;
			end
			
			%		step2b -- RO: define the new pchg in demand based on ste
			DemandChange(2) = DemandChangeADJ(kkkk-1,2) - DemandStep(2);

			% 		step2c -- RO: check if stuck cycling -- explain?
			if     ((kkkk > 4) && (DemandChange(2) == DemandChangeADJ(kkkk-2,2)) && (DemandChangeADJ(kkkk-1,2) == DemandChangeADJ(kkkk-3,2)))
				DemandChange(2) = (DemandChangeADJ(kkkk-2,2) + DemandChangeADJ(kkkk-1,2))/2;
			elseif ((kkkk > 10) && (DemandChange(2) < DemandChangeADJ(kkkk-1,2)) && (DemandChangeADJ(kkkk-1,2) > DemandChangeADJ(kkkk-2,2)) && (DemandChangeADJ(kkkk-2,2) < DemandChangeADJ(kkkk-3,2)))
				DemandChange(2) = (DemandChangeADJ(kkkk-2,2) + DemandChangeADJ(kkkk-1,2))/2;
			elseif ((kkkk > 10) && (DemandChange(2) > DemandChangeADJ(kkkk-1,2)) && (DemandChangeADJ(kkkk-1,2) < DemandChangeADJ(kkkk-2,2)) && (DemandChangeADJ(kkkk-2,2) > DemandChangeADJ(kkkk-3,2)))
				DemandChange(2) = (DemandChangeADJ(kkkk-2,2) + DemandChangeADJ(kkkk-1,2))/2;
			end


	   end %iteration > 1

		% 		step 3: save requested and adjusted demand changes
		DemandChangeADJ = [DemandChangeADJ; DemandChange];

	end %iteration = 1

	%		display adjusted percentage change in electricity demand (from rem)
	disp('demand change - Adjusted');
	disp(DemandChange);

	%		save percentage change in electricity demand as Matlab data
	save('DemandChangeData.mat', 'DemandChangeCGE', 'DemandChangeADJ');

	% --------------------------------------------------------------------
    % export adjusted percentage change in electricity demand for next
	% iteration
    % --------------------------------------------------------------------

	%		define header names
	EconHeader     = {'DemandChangeCali','DemandChangeRO'};           %define header name
	commaHeader    = [EconHeader;repmat({','},1,numel(EconHeader))];  %define header name with commas
	commaHeader    = commaHeader(:)';                                 %something
	EcontextHeader = cell2mat(commaHeader);                           %convert header names to matrix
	
	%		export: write header to file in main results 
	filename = ['Input_data', '.csv'];                                %deifine file name
	fid = fopen(filename,'w');                                        %being writing
	fprintf(fid,'%s\n',EcontextHeader');                              %write header name to file
	fclose(fid);                                                      %close file writing

	%		export: write data to end of file (under header)
	dlmwrite(filename,DemandChange,'-append');                        %append data to file
	
    % 		change working directory to main coupled model directory
	cd ..	%  move to main coupled model directory: ./codedir

    % save adjusted percentage again in psm model directory)
	savefilename = sprintf('results/Iteration%d/Input_data.csv', IterationCount); %display iteration main results directory
	copyfile('Gams_data/Input_data.csv', savefilename);                           %save psm model directory

	% --------------------------------------------------------------------
    % export main coupled model results to a csv file
    % --------------------------------------------------------------------

	%		define main results data matrix
	z(1)  = YearlyCaliCoal;       %CA - total annual electricity generation for fuel type coal
	z(2)  = YearlyCaliGas;        %CA - total annual electricity generation for fuel type natural gas (CT and CCGT)
	z(3)  = YearlyROCoal;         %RO - total annual electricity generation for fuel type coal
	z(4)  = YearlyROGas;          %RO - total annual electricity generation for fuel type natural gas (CT and CCGT)
	z(5)  = YearlyGenCali;        %CA - total annual electricity generation
	z(6)  = YearlyGenRO;          %RO - total annual electricity generation
	z(7)  = CaliAvgElectricPrice; %CA - total annual electricity price
	z(8)  = ROAvgElectricPrice;   %RO - total annual electricity price
	z(9)  = AvgMarginalCostCali;  %CA - electricity generation weighted average marginal cost
	z(10) = AvgMarginalCostRO;    %RO - electricity generation weighted average marginal cost
	z(11) = ValueCali;            %CA - total value of generation by hour sum(d, MC_mean(d)*GEN_(d))
	z(12) = ValueRO;              %RO - total value of generation by hour sum(d, MC_mean(d)*GEN_(d))
	z(13) = YearlyFuelCoalCali;   %CA - total annual fuel usage for fuel type coal
	z(14) = YearlyFuelGasCali;    %CA - total annual fuel usage for fuel type natural gas (CT and CCGT)
	z(15) = YearlyFuelCoalRO;     %RO - total annual fuel usage for fuel type coal
	z(16) = YearlyFuelGasRO;      %RO - total annual fuel usage for fuel type natural gas (CT and CCGT)
	z(17) = AvgMarginalCostCali;  %CA - electricity generation weighted average marginal cost
	z(18) = AvgMarginalCostRO;    %RO - electricity generation weighted average marginal cost
	z(19) = ChangeCali;           %CA - percentage change in electricity generation weighted average marginal cost (electricity productivity impact for REM)    
	z(20) = ChangeRO;             %RO - percentage change in electricity generation weighted average marginal cost (electricity productivity impact for REM)    

	%		define second main results matrix
	y = z;

	%		define header names
	EconHeader     = {'CaliCoal','CaliGas','RoweccCoal','RoweccGas','GenCali','GenRO','AvgElectricPriceCali','AvgElectricPriceRO','marginalElectricityPriceCali','marginalElectricityPriceRO','ValueCali','ValueRO','YearlyFuelCoalCali','YearlyFuelGasCali','YearlyFuelCoalRO','YearlyFuelGasRO','AvgMarginalCostCali','AvgMarginalCostRO','ChangeCali','ChangeRO','CaliDemandChange','RODemandChange','PriceChangeCali','PriceChangeRO','ArmPriceChangeCali','ArmPriceChangeRO'};
	commaHeader    = [EconHeader;repmat({','},1,numel(EconHeader))];
	commaHeader    = commaHeader(:)';
	EcontextHeader = cell2mat(commaHeader);
	
	% 		change working directory to current iteration in main coupled model results directory
	cd results		%  move to main coupled model directory: ./codedir/results
	
	filename = ['Iteration', int2str(IterationCount)];
	cd(filename)    %  move to main coupled model directory: ./codedir/results/iterationX
	
	%		export: write header to file -- main results
	filename = ['Shares_Price', int2str(IterationCount), '.csv'];
	fid      = fopen(filename,'w');
	fprintf(fid,'%s\n',EcontextHeader);
	fclose(fid);

	%		export: write data to end of file (under header)  -- main results
	dlmwrite(filename,y,'-append');

	%		export: write header to file and data -- NSE
	filename = ['TotalUnserved', int2str(IterationCount), '.csv'];
	fid = fopen(filename,'w');
	fprintf(fid,'%d,%d\n',UnservedTotalCA,UnservedTotalRO);
	fclose(fid);

	% --------------------------------------------------------------------
    % change working directory back to main coupled model directory
    % --------------------------------------------------------------------
	
	cd .. %./codedir/results
	cd .. %./codedir
	
	% --------------------------------------------------------------------
    % change working directory to current iteration directory in detailed
	% psm results directory
    % --------------------------------------------------------------------
	
	cd Results_Each_Iter %./codedir/Results_Each_Iter

	filename = ['Iteration', int2str(IterationCount)];
	cd(filename)         %  move to main coupled model directory: ./codedir/Results_Each_Iter/iterationX
	
	%		export: write header to file and data -- 
	%		percentage change in electricity generation weighted average 
	%		marginal cost (electricity productivity impact for REM)    		
	Temp     = [ChangeCali ChangeRO];
	filename = ['Change', int2str(IterationCount), '.csv'];
	csvwrite(filename,Temp);
	
	% --------------------------------------------------------------------
    % change working directory back to main coupled model directory
    % --------------------------------------------------------------------
		
	cd .. %./codedir/Results_Each_Iter
	cd .. %./codedir

	% --------------------------------------------------------------------
    % get ACI job id and export to csv file
    % --------------------------------------------------------------------

	myjobid=getenv('PBS_JOBID');
	filename = 'CheckConvgJob.csv';
	csvwrite(filename,myjobid);

	% --------------------------------------------------------------------
    % check for convergence
    % --------------------------------------------------------------------
	%		write a csv file which will then be used to Check for convergence. 
	% 
	%and .
	if (abs(ChangeCali)<0.0001)&&(abs(ChangeRO)<0.0001) % if converges, exit status of true -- 
		fid = fopen('outfile','wt');                    % 		then the bash scripts will terminate
		fprintf(fid, 'True');
		fclose(fid);
		rc=200;
	else                                                % else not converge, exit status of false --
		fid = fopen('outfile','wt');                    % 		then the bash scripts will continue
		fprintf(fid, 'False');
		fclose(fid);
		rc=201;
	end

	quit(rc); % true/false statue for convergance 

end
