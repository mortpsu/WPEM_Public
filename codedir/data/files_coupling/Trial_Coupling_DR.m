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
    % Add directory path for functions -- this cleans up main folder
    % --------------------------------------------------------------------

    addpath(genpath('./scripts/functions_psm/'));

    % --------------------------------------------------------------------
    % remove when on ACI: this is the test a single week
    % --------------------------------------------------------------------
    
%    IterationCount = 1;

	% --------------------------------------------------------------------
    % define working directory for iteration
    % --------------------------------------------------------------------
	
	% 		define working directory name: ./IterationX
	filename = ['Iteration',int2str(IterationCount)];
	
	% 		move to psm detailed results directory: ./codedir/Results_Each_Iter
	cd Results_Each_Iter
	
	%       working directory: ./codedir/Results_Each_Iter/IterationX
	cd(filename)
	
	% --------------------------------------------------------------------
    % export csv of iteration count
    % --------------------------------------------------------------------
	
	filename = ['Count', int2str(IterationCount), '.csv']; 		
	csvwrite(filename,IterationCount);                       
	
	% --------------------------------------------------------------------
    % import main results csv: See documentation for results list
    % --------------------------------------------------------------------

	% 		loop through week w and import main results	
	for week=1:1:52
        
        if (week == 1) % define main results matrix by importing the first week
            
            disp(week);                                    % display week i
            filename = ['week', int2str(week), '.csv'];    % define filename
            z1       = csvread(filename);                  % read csv file
            
        else          % import week i > 1 main results
            
            disp(week);							        	
            filename = ['week', int2str(week), '.csv'];   
            z2       = csvread(filename);               

            % concatenate imported week i to main array 
            z1 = cat(1,z1,z2);
        end
          
	end
	
	% --------------------------------------------------------------------
    % define annual values for costs, generation, and fuel usage
    % --------------------------------------------------------------------
    
    %       define array without economic region id column
    z2 = z1(:,2:end);
    
    %       replace NaN values with zeros
    z2(isnan(z2)) = 0;
    
    %       define rows associated with economic region id (XX) do this all
    %       columns (YY)
    [xx, yy] = ndgrid(z1(:,1),1:size(z2,2));
    
    %       defined annual total for main economic region results 
    %       (add economic region id as first column)
    z3 = [unique(z1(:,1)), accumarray([xx(:) yy(:)],z2(:))];
  
    %       calculate annualized cost of electicity generation ($/MWH) 
    %       at end of main results
    z3(:,end+1) = z3(:,3)./z3(:,4);
    
    %       define final matrix of main results
    z = z3;
    
	% --------------------------------------------------------------------
    % import second set of main results csv files
    % --------------------------------------------------------------------

    %       pre-allocate
    mc_ird   = cell(1,52);
    mdem_ird = cell(1,52);
    nse_ird  = cell(1,52);
    shed_ird = cell(1,52);
    gen_rd   = cell(1,52);
    dem_rd   = cell(1,52);
    nse_rd   = cell(1,52);
    
	% 		loop through each week and import results	
	for week =1:1:52
		
		if (week == 1) 
			% 		csv file for marginal cost  - bus i in economic region r - hour d
			filename = ['MarginalCost', int2str(week), '.csv'];
			mc_ird{week}   = csvread(filename);
			
			% 		csv file for electricity demand - bus i in economic region r - hour d
			filename = ['MultipliedDemand', int2str(week), '.csv'];
			mdem_ird{week} = csvread(filename);
			
			% 		csv file for non-served electricity at bus i for hour d
			filename = ['UnservedEnergywhole', int2str(week), '.csv'];
			nse_ird{week}  = csvread(filename);

			% 		csv file for shed energy - bus i - hour d
			filename = ['ShedEnergywhole', int2str(week), '.csv'];
			shed_ird{week} = csvread(filename);
			
			% 		csv file for electricity generation - economic region r  - hour d
			filename = ['Genperhour', int2str(week), '.csv'];
			gen_rd{week}   = csvread(filename);
				
			% 		csv file for electricity demand - economic region r - hour d
			filename = ['NewDemandwhole', int2str(week), '.csv'];
			dem_rd{week}   = csvread(filename);

			% 		csv file for non-served electricity - economic region r - hour d
			filename = ['NewEnergywholeunserved', int2str(week), '.csv'];
			nse_rd{week}   = csvread(filename); 
		else
			% 		csv file for marginal cost  - bus i in economic region r - hour d
			filename = ['MarginalCost', int2str(week), '.csv'];
			tmp      = csvread(filename);
			mc_ird{week} = tmp(:,3:end);
			
			% 		csv file for electricity demand - bus i in economic region r - hour d
			filename = ['MultipliedDemand', int2str(week), '.csv'];
			tmp      = csvread(filename);
			mdem_ird{week} =tmp(:,3:end);
			
			% 		csv file for non-served electricity at bus i for hour d
			filename = ['UnservedEnergywhole', int2str(week), '.csv'];
			tmp      = csvread(filename);
			nse_ird{week} =tmp(:,3:end);
			
			% 		csv file for shed energy - bus i - hour d
			filename = ['ShedEnergywhole', int2str(week), '.csv'];
			tmp      = csvread(filename);
			shed_ird{week} =tmp(:,3:end);
			
			% 		csv file for electricity generation - economic region r  - hour d
			filename = ['Genperhour', int2str(week), '.csv'];
			tmp      = csvread(filename);
			gen_rd{week} =tmp(:,2:end);
				
			% 		csv file for electricity demand - economic region r - hour d
			filename = ['NewDemandwhole', int2str(week), '.csv'];
			tmp      = csvread(filename);
			dem_rd{week} =tmp(:,2:end);

			% 		csv file for non-served electricity - economic region r - hour d
			filename = ['NewEnergywholeunserved', int2str(week), '.csv'];
			tmp      = csvread(filename);
			nse_rd{week} = tmp(:,2:end);
		end
	end

	% --------------------------------------------------------------------
    % return to main coupled model directory
    % --------------------------------------------------------------------
	
	cd ..   % move to psm detailed results directory: ./codedir/Results_Each_Iter
	cd ..   % move to main coupled model directory:   ./codedir

	% --------------------------------------------------------------------
    % calculate electricity productivity parameter for the regional 
	% economic model (REM)
    % --------------------------------------------------------------------
	
	% 		convert import data from cell to matrix format
	annual_mc_ird    = cell2mat(mc_ird);    % marginal cost(i,r,d) x weeks
	annual_nse_ird   = cell2mat(nse_ird);   % non-served electicity(i,r,d) x weeks
	annual_shed_ird  = cell2mat(shed_ird);  % shed electicity(i,r,d) x weeks
	annual_gen_rd    = cell2mat(gen_rd);    % generation(r,d) x weeks
	
	%		define index for columns 
	lin_idx = find(annual_shed_ird(:,3:end)>0);
    
	[~, mat_idx]  = ind2sub(size(annual_shed_ird(:,3:end)), lin_idx);

	%		remove columns with NSE from marginal cost and generation
    tmp_mat = annual_mc_ird(:,3:end);                  % create temp array of i x d
    tmp_mat(:,mat_idx) = [];                           % remove shed column for d
    annual_mc_ird = [annual_mc_ird(:,1:2),  tmp_mat];  % combine set i,r cols with i x d cols 
    
    tmp_mat = annual_nse_ird(:,3:end);
    tmp_mat(:,mat_idx) = [];
    annual_nse_ird = [annual_nse_ird(:,1:2),  tmp_mat];
    
	tmp_mat = annual_gen_rd(:,2:end);                  % create temp array of r x d
    tmp_mat(:,mat_idx) = [];                           % remove shed column for d
    annual_gen_rd = [annual_gen_rd(:,1),  tmp_mat];    % combine set r col with r x d cols 
    
	% --------------------------------------------------------------------
	%		identify non-served electricity hours and replace the marginal
	%       cost with penilty marginal cost
	% --------------------------------------------------------------------
	
	%		define crosswalks
%	ieconr    = unique(annual_nse_ird(:,1:2),'rows');
%	num_bus   = unique(annual_nse_ird(:,1));
	id_econr = unique(annual_nse_ird(:,2));           % economic regions id
	
	%		define index for rows and columns with NSE
	lin_idx = find(annual_nse_ird(:,3:end)>0);
	[mat_idx_row, mat_idx_col]  = ind2sub(size(annual_nse_ird(:,3:end)), lin_idx);
	
	%		define index for buses in economic region
%	[M,~] = ismember((ieconr'),annual_nse_ird(:,2));
	
	%		define rows and correspond columns for economic regions
%	econr_idx_row = mat_idx_row(M);
%	econr_idx_col = mat_idx_col(M);


	%		define the maximum marginal cost for nse
%  	max_mc_nse_r = [1:1:14;ones(1,14)*268]';	
% 	[M,~] = ismember((ieconr'),max_mc_nse_r(:,2));

	%		replace columns with maximum MP for marginal cost and generation
%   annual_mc_ird(econr_idx_row,econr_idx_col) = max_mc_nse_r(M,2);
%	annual_mc_ird(econr_idx_row,econr_idx_col) = 268;

    tmp_mat = annual_mc_ird(:,3:end);
    tmp_mat(mat_idx_row,mat_idx_col) = 268;            
    annual_mc_ird = [annual_mc_ird(:,1:2), tmp_mat];

	% --------------------------------------------------------------------
    %  define reference for calculate electricity productivity parameter for the regional 
	% economic model (REM)
    % --------------------------------------------------------------------
	
	% 		define current reference prices for calculating percentage 
	% 		change from the baseline
     ref_mc_r = [36.0925, 35.2093, 33.1209, 32.0904, 33.9715, ...
                 35.4242, 34.6744, 34.7662, 32.2148, 34.4561, ...
                 32.9551, 35.4318, 33.3362, 35.7772]';

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

	% --------------------------------------------------------------------
    % calculate electricity productivity parameter for the regional 
	% economic model (REM)
    % --------------------------------------------------------------------
    
    %       define generation matrix array without economic region column
    val_gen_rd = annual_gen_rd(:,2:end);
    
	%		calculate the mean of marginal cost for each hour: MC_mean(d)
	% define array without economic region id column
    val_mc_ird = annual_mc_ird(:,3:end);
	
	% define rows associated with economic region id (XX) do this all
    % columns (YY)
    [xx, yy] = ndgrid(annual_mc_ird(:,2),1:size(val_mc_ird,2));
	
	% mean over rows (do not add economic region column
	mean_mc_rd = accumarray([xx(:) yy(:)],val_mc_ird(:),[],@(x) mean(x));
	
	%		calculate total value of generation by hour sum(d, MC_mean(d)*GEN_(d))
	mean_toc_r = sum(mean_mc_rd.*val_gen_rd,2);

	%		calculate electricity generation weighted average marginal cost
	annual_mean_mc_r = mean_toc_r./sum(val_gen_rd,2);
	
	%		define import previous iteration marginal cost reference price     
 	if (IterationCount == 1)
 		prev_mc_r = ref_mc_r;
	else
	  % load previous iteration main results file
		filename = ['Results_Each_Iter/Iteration', int2str(IterationCount-1),'/Change',int2str(IterationCount-1),'.csv'];
		prev_iter_main_results = csvread(filename,1,0);
		prev_mc_r = prev_iter_main_results(:,end);
 	end
	
	%		calcuate percentage change in electricity generation weighted 
	%       average marginal cost from the baseline    
	ref_chang_r    = (annual_mean_mc_r./ref_mc_r)-1;
	ex_ref_chang_r = [id_econr, ref_chang_r];
	
    %		calcuate percentage change in electricity generation weighted 
	%       average marginal cost from the previous iteration    
	prev_chang_r    = (annual_mean_mc_r./prev_mc_r)-1;
	ex_prev_chang_r = [id_econr, prev_chang_r];
	
    % --------------------------------------------------------------------
    % add generation weighted average marginal cost and percentage change
    % to main results
    % --------------------------------------------------------------------
    
    %       reference generation weighted average marginal cost
    z(:,end+1)  = NaN; 
    z(:,end)    = ref_mc_r;
    
    %       generation weighted average marginal cost
    z(:,end+1)  = NaN; 
    z(:,end)    = annual_mean_mc_r; 
    
    %       percentage change from baseline in generation weighted average
    %       marginal cost 
    [zpos,zidx] = ismember(z(:,1),ex_ref_chang_r(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = ex_ref_chang_r(zidx(zpos),end);
    clear zpos zidx
    
    %       percentage change from previous in generation weighted average
    %       marginal cost 
    [zpos,zidx] = ismember(z(:,1),ex_prev_chang_r(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = ex_prev_chang_r(zidx(zpos),end);
    clear zpos zidx
	
	% --------------------------------------------------------------------
    % calcuate non-electricity (Sullivan) productivity impacts from NSE
	% for REM
    % --------------------------------------------------------------------

	% 		convert import demand and NSE from cell to matrix format
	annual_dem_rd = cell2mat(dem_rd); % (econr*hour x 1)
	annual_nse_rd = cell2mat(nse_rd); 

	%		call script: reshape function hour x region matrix to single
	%		             vector --?
	annual_dem_dr = annual_dem_rd';
	annual_nse_dr = annual_nse_rd';

	%		combine demand and NSE into a single matrix
	annual_nse_data = [annual_dem_dr(2:end,:), annual_nse_dr(2:end,:)];
	
	% 		calculate total annual NSE
	annual_nse_r = [annual_nse_rd(:,1), sum(annual_nse_rd(:,2:end),2)];   

	%		call script: to lable hourly date by hour, date, time of day, 
	%                    and seasonally. script creates a csv file--?
	output = Datetimestamp(annual_nse_data,IterationCount);

    % 		call script: non-electricity (Sullivan) productivity impacts 
	%                    from NSE
 	[sector_losses_ALL] = CalcSullivanLosses(IterationCount, annual_nse_data);
	
	% --------------------------------------------------------------------
    % create gams gdx files for REM: MtoCGE.gdx and MtoCGE2.gdx
	%
	% (MtoCGE == electricity productivity impacts == ; 
	%  MtoCGE2 == non-electricity productivity impacts)                          
    % --------------------------------------------------------------------	

	%		change working directory to the REM directory
	cd modelcge
    
    % define number of economic regions present in DREM
    econr_drem = size(sector_losses_ALL,2);
	
    %		define structure of parameter 
	econr.name    = 'econr';
	econr.uels    = {1:econr_drem};
	nele.name     = 'nele';
	nele.uels     = {'osa', 'grn', 'vna', 'oca', 'pfb', 'cba', 'apa', 'frs', 'con', 'fin', 'fbm', 'bom', 'wpm', 'per', 'cpm', 'cem', 'pfm', 'tec', 'trm', 'cru', 'coa', 'min', 'pub', 'hlt', 'bos', 'ngd', 'tel', 'rtl', 'trn'};
	% nele.uels = {'AGR', 'MIN', 'CONST', 'MANUF', 'TEL_UTL', 'TRD_RTL', 'FIN', 'SRV', 'PUB'}; % for static REM (SREM)

	% define parameter for non-electricity productivity impacts
	sull_impacts.name = 'sull_impacts';
	sull_impacts.type = 'parameter';
    sull_impacts.val  = sector_losses_ALL;
	sull_impacts.form = 'full';
	sull_impacts.uels = {nele.uels,econr.uels};
	
	% export non-electricity productivity impacts gdx
	wgdx ('MtoCGE', sull_impacts);
	
	ele_impacts.name = 'ele_impacts';
	ele_impacts.type = 'parameter';
	ele_impacts.val  = ex_ref_chang_r(1:econr_drem,:);
	ele_impacts.form = 'sparse';
	ele_impacts.uels = {econr.uels} ;
	
	% export electricity productivity impacts gdx
	wgdx('MtoCGE2', ele_impacts);

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
	switch_coupled = 5;       % 0 == NONE , 
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
    % import REM results from gdx (CCEtoM.gdx): demand, electricity price
    % and armington electricity price 
    % --------------------------------------------------------------------

	% crosswalk between economic regions in DREM and zones      
    tmpstruct.name = 'econz_map';
    tmpstruct.form = 'sparse';                                       
	tmpmat         = rgdx('CGEtoM',tmpstruct);
    econrz = tmpmat.val; % economic region id and zones id
    
    y = econrz;          % define temp crosswalk array
    
    clear tmpstruct tmpmat;
    
    % percentage change from the baseline: economic region
	%	sectoral production for the domestic market 
	%	region CA for sector electricity at time period 2050
    % electricity generation (MWH) - economic region r - week
    tmpstruct.name = 'd_ele_econr_pchg';
    tmpstruct.form = 'sparse';                                        
	tmpmat         = rgdx('CGEtoM',tmpstruct);     
    [zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
    
    demand_chg_econr = tmpmat.val'; % define current change in economic region demand from the REM
    
	clear tmpstruct tmpmat zpos zidx;

    % percentage change from the baseline: economic region
	%	sectoral production for the domestic market 
	%	region RO for sector electricity at time period 2050
    tmpstruct.name = 'p_ele_econr_pchg';
    tmpstruct.form = 'sparse';                                        
	tmpmat         = rgdx('CGEtoM',tmpstruct);     
    [zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;
    
    
    % percentage change from the baseline: economic region
	%	armington aggregate price 
	%	region CA for sector electricity at time period 2050
    tmpstruct.name = 'pa_ele_econr_pchg';
    tmpstruct.form = 'sparse';                                        
	tmpmat         = rgdx('CGEtoM',tmpstruct);     
    [zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;
    
	% percentage change from the baseline: psm zones
	%	sectoral production for the domestic market 
	%	region CA for sector electricity at time period 2050
    % electricity generation (MWH) - economic region r - week
    tmpstruct.name = 'd_ele_zones_pchg';
    tmpstruct.form = 'sparse';                                        
	tmpmat         = rgdx('CGEtoM',tmpstruct);     
	[ypos,yidx] = ismember(y(:,2),tmpmat.val(:,1));  % combine zone data with economic region crosswalk    
	y(:,end+1)  = NaN;                               % add NaN column at end of temp arrya
	y(ypos,end) = tmpmat.val(yidx(ypos),end);        % merge imported parameter value column to temp array
    
    demand_chg_zones = tmpmat.val';                  % define current change in demand from the REM
    
	clear tmpstruct tmpmat ypos yidx;  
    
	% percentage change from the baseline: psm zones
	%	sectoral production for the domestic market 
	%	region RO for sector electricity at time period 2050
    tmpstruct.name = 'p_ele_zones_pchg';
    tmpstruct.form = 'sparse';                                        
	tmpmat         = rgdx('CGEtoM',tmpstruct);     
	[ypos,yidx] = ismember(y(:,2),tmpmat.val(:,1));   
	y(:,end+1)  = NaN;                              
	y(ypos,end) = tmpmat.val(yidx(ypos),end);   
	clear tmpstruct tmpmat ypos yidx;  
	
    
	% percentage change from the baseline: psm zones
	%	armington aggregate price 
	%	region CA for sector electricity at time period 2050
	tmpstruct.name = 'pa_ele_zones_pchg';
    tmpstruct.form = 'sparse';                                        
	tmpmat         = rgdx('CGEtoM',tmpstruct);     
	[ypos,yidx] = ismember(y(:,2),tmpmat.val(:,1));   
	y(:,end+1)  = NaN;                              
	y(ypos,end) = tmpmat.val(yidx(ypos),end);      
	clear tmpstruct tmpmat ypos yidx;  

	% --------------------------------------------------------------------
    % change working directory to the psm directory
    % --------------------------------------------------------------------
	
	% 		move to main coupled model directory: ./codedir
	cd ..

	% --------------------------------------------------------------------
    % define stepsize logic for convergence
    % --------------------------------------------------------------------

	%		display current percentage change in electricity demand (from rem)
	disp('demand change for zones from DREM');
	disp(demand_chg_zones);
    
    disp('demand change for economic regions from DREM');
	disp(demand_chg_econr);
	
    %       call demand change adjustment function
    demand_chg_adj_zones = DemandAdj(IterationCount,demand_chg_zones,'zones');
    demand_chg_adj_econr = DemandAdj(IterationCount,demand_chg_econr,'econr');
    
    % --------------------------------------------------------------------
    % change working directory to the psm directory
    % --------------------------------------------------------------------
    
    % 		move to main coupled model directory: ./codedir/Gams_data
	cd Gams_data

	% --------------------------------------------------------------------
    % export adjusted percentage change in electricity demand for next
	% iteration
    % --------------------------------------------------------------------

    %       export: demand adjustment change for zones
	% time headers
    cell_econr = cellstr(string(demand_chg_zones(1,:)));                        
    EconHeader = matlab.lang.makeValidName(cell_econr,'Prefix','DemandChang_'); % write demand header names per zone
    commaHeader    = [EconHeader;repmat({','},1,numel(EconHeader))];            % define header name with commas
	commaHeader    = commaHeader(:)';                                           % something
	EcontextHeader = cell2mat(commaHeader);                                     % convert header names to matrix
    
    % write header to file in main results 
	filename = ['Input_data_UC', '.csv'];           %deifine file name
	fid = fopen(filename,'w');                   %begin writing
	fprintf(fid,'%s\n',EcontextHeader');         %write header name to file
	fclose(fid);                                 %close file writing
    
    % write data to end of file (under header)
	dlmwrite(filename,demand_chg_adj_zones,'-append');   %append data to file
    
    %       export: demand adjustment change for economic region
    % define header title ids
    cell_econr = cellstr(string(demand_chg_econr(1,:)));
    EconHeader = matlab.lang.makeValidName(cell_econr,'Prefix','DemandChang_'); % write demand header names per region
    commaHeader    = [EconHeader;repmat({','},1,numel(EconHeader))]; 
	commaHeader    = commaHeader(:)';                                 
	EcontextHeader = cell2mat(commaHeader);                         

	% write header to file in main results 
	filename = ['Input_data_DC', '.csv'];                               
	fid = fopen(filename,'w');                                       
	fprintf(fid,'%s\n',EcontextHeader');                             
	fclose(fid);                                                     

	% write data to end of file (under header)
	dlmwrite(filename,demand_chg_adj_econr,'-append');                      
	
    % 		change working directory to main coupled model directory
	cd ..	%  move to main coupled model directory: ./codedir

    % save adjusted percentage again in psm model directory)
	savefilename = sprintf('results/Iteration%d/Input_data_UC.csv', IterationCount); %display iteration main results directory
	copyfile('Gams_data/Input_data_UC.csv', savefilename);                           %save psm model directory
    
    savefilename = sprintf('results/Iteration%d/Input_data_DC.csv', IterationCount); %display iteration main results directory
	copyfile('Gams_data/Input_data_DC.csv', savefilename);                           %save psm model directory

    % --------------------------------------------------------------------
    % change working directory to current iteration in main coupled model
    % results directory 
    % --------------------------------------------------------------------
    
	cd results		%  move to main coupled model directory: ./codedir/results
    
	% --------------------------------------------------------------------
    % export main coupled model results to a csv file:
    % --------------------------------------------------------------------

	%		define header names
	EconHeader     = {'econr', ...
                      'nse_r','toc_r','gen_r','coal_gen_r','ngas_gen_r', ...
                      'steam_gen_r','nuclear_gen_r','biogas_gen_r','geo_gen_r','petrol_gen_r','intertiel_gen_r','wind_gen_r', ...
                      'hydro_gen_r','pshydro_gen_r','motorload_gen_r','solar_gen_r', ...
                      'toc_coal_r','toc_ngas_r', ...
                      'annual_mc_r', ...
                      'ref_mc_r','annual_mean_mc_r','ex_ref_chang_r','ex_prev_chang_r', ...
                      'd_ele_econr_pchg','p_ele_econr_pchg','pa_ele_econr_pchg'};
	commaHeader    = [EconHeader;repmat({','},1,numel(EconHeader))];
	commaHeader    = commaHeader(:)';
	EcontextHeader = cell2mat(commaHeader);
	
	filename = ['Iteration', int2str(IterationCount)];
	cd(filename)    %  move to main coupled model directory: ./codedir/results/iterationX
	
	%		export: write header to file -- main results
	filename = ['Shares_Price', int2str(IterationCount), '.csv'];
	fid      = fopen(filename,'w');
	fprintf(fid,'%s\n',EcontextHeader);
	fclose(fid);

	%		export: write data to end of file (under header)  -- main results
	dlmwrite(filename,z,'-append');
    
    % --------------------------------------------------------------------
    % export main coupled model results to a csv file: 
    % --------------------------------------------------------------------
  
    EconHeader     = {'econr', 'zones','d_ele_zones_pchg','p_ele_zones_pchg','pa_ele_zones_pchg'};
	commaHeader    = [EconHeader;repmat({','},1,numel(EconHeader))];
	commaHeader    = commaHeader(:)';
	EcontextHeader = cell2mat(commaHeader);
    
    filename = ['Shares_Price_z', int2str(IterationCount), '.csv'];
	fid      = fopen(filename,'w');
	fprintf(fid,'%s\n',EcontextHeader);
	fclose(fid);
    
    dlmwrite(filename,y,'-append');
    
    % --------------------------------------------------------------------
    % export main coupled model results to a csv file
    % --------------------------------------------------------------------
    
    EconHeader     = {'econr', 'annual_nse_r'};
	commaHeader    = [EconHeader;repmat({','},1,numel(EconHeader))];
	commaHeader    = commaHeader(:)';
	EcontextHeader = cell2mat(commaHeader);

	%		export: write header to file and data -- NSE
	filename = ['TotalUnserved', int2str(IterationCount), '.csv'];
	fid = fopen(filename,'w');
	fprintf(fid,'%s\n',EcontextHeader);
	fclose(fid);
    
    dlmwrite(filename,annual_nse_r,'-append');

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
	EconHeader     = {'econr', 'ex_ref_chang_r', 'ex_prev_chang_r', 'annual_mean_mc_r'};
	commaHeader    = [EconHeader;repmat({','},1,numel(EconHeader))];
	commaHeader    = commaHeader(:)';
	EcontextHeader = cell2mat(commaHeader);

	%		export: write header to file and data -- NSE
	filename = ['Change', int2str(IterationCount), '.csv'];
	fid = fopen(filename,'w');
	fprintf(fid,'%s\n',EcontextHeader);
	fclose(fid);
    
    Temp = [ex_ref_chang_r ex_prev_chang_r(:,end) annual_mean_mc_r];
    dlmwrite(filename,Temp,'-append');
	
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
    % 
    % --------------------------------------------------------------------
	
    %       define convergence criterion
    CC = 0.0003;
    
    %       create vector of convergence criterion to match economic
    %       regions
    C = repmat(CC,length(id_econr),1);
    
    %       define percentage change from baseline
    change_r = ex_prev_chang_r(:,end);
    
    %       write a csv file which will then be used to Check for convergence
	if (abs(change_r)<C)                              % if converges, exit status of true -- 
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
