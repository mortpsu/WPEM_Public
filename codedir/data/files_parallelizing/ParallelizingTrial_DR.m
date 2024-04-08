%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Parallelizing PSM
% I.  This program imports the necessary data - demand in 2010, power 
%     outages from the WBM, generator data (including availability
% II. Export data to .gdx 
% III. Run Unit Commitment and DCOPF models
% 
%     i - the bus ID number: 1 - 312 
%     g - the generator number: 1 - 3569
%  zone - "states": 1 - 13
% econr - economic regions between PSM and DREM 
%         Not including CAN and MEX    4, 6, 14
%             Including CAN and MEX2   2, 4, 11 or 10 (OR and WA combined)
%  line - transmission lines between buses - 
% lines - transmission lines between zones - 36
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z=ParallelizingTrial(week, workingdir)

	% --------------------------------------------------------------------
    % define working directory for parallelizing the PSM
    % --------------------------------------------------------------------
	
    newFolder = strtrim(workingdir);
    fprintf('Running ParallelizingTrial for week %2d in directory %40s\n', week, newFolder)
    oldFolder = cd(newFolder);
    fprintf('Switching from %80s to %80s\n', oldFolder, newFolder)

    % --------------------------------------------------------------------
    % remove when on ACI: this is the test a single week
    % --------------------------------------------------------------------
	
% 	  cd Gams_data
%     week = 28;
%     week = 35;
	
	% --------------------------------------------------------------------
    % begin timer
    % --------------------------------------------------------------------
	
    tic
	
	% --------------------------------------------------------------------
    % load percentage change in demand from the economic model (initially
    % set to zeros
    % --------------------------------------------------------------------
    
    %		define demand multipler for the unit commitment program
	fid = fopen('Input_data_UC.csv'); 
	Input_data = textscan(fid, '%f','headerLines',1,'delimiter',',','collectoutput',1);
	fclose(fid);
	DemandMult_UC      = (Input_data{1,1})'; 
    
    %       define demand multipler for the DCOPF program
    fid = fopen('Input_data_DC.csv'); 
	Input_data = textscan(fid, '%f','headerLines',1,'delimiter',',','collectoutput',1);
	fclose(fid);
	DemandMult_DC       = (Input_data{1,1})'; 

	% --------------------------------------------------------------------
    % load WBM power outage data
    % --------------------------------------------------------------------
    
	load('PowerOutages.mat', 'PowerOutagesbyGen', 'GeneratorsOff')
    % load('PowerOutages.mat', 'PowerOutagesbyGenCap')      % DATA NOT USED

	% --------------------------------------------------------------------
    % load demand by bus for the whole year
    % --------------------------------------------------------------------
    
	fid=fopen('Demandfor2010.csv');
    demand2012CA=textscan(fid, ['%s',repmat('%f',[1,8784])],'headerlines',1,'delimiter',',','collectoutput',1);
    fclose(fid);
    demand = demand2012CA{1,2};

	% --------------------------------------------------------------------
    % load demand response from the demand response emulator by bus for the whole year
    % --------------------------------------------------------------------
	
	% delete extra day from base WECC demand
	demand(:,8761:end) = [];

	% Load the demand response factors for this case
	load('load_scale_factors.mat','busdata');
	
	% modify demand(bus,hour) by DR factor
	demand = demand .* busdata;
	
	% --------------------------------------------------------------------
    % define the start and final hour of the week
    % --------------------------------------------------------------------
    
	FirstHour  = 1+(week-1)*168;
    LastHour   = week*168;
	
    % 		define demand by bus for the specific week	
    Demandweek = demand(:,FirstHour:LastHour);

	% --------------------------------------------------------------------
    % alter week 48: maintain feasibility in this
    % particular week to meeting demand - this can be changed but will
    % require more fine tuning of the models than is currently present
    % --------------------------------------------------------------------
    
    if week==48
        
        for ii=2:1:168
            if Demandweek(199,ii)>9000
                Demandweek(199,ii) = Demandweek(199,ii-1);
            end
        end
        
        Demandweek(202,:)            = Demandweek(202,:)-200;
        columns=[152 345; 153 150; 154 27; 161 94; 155 52; 156 110; 162 182; 163 289];
        
        Demandweek(98,columns(:,1))  = Demandweek(98,columns(:,1))-columns(:,2)';
        columns=[154 51 161 60];
        
        Demandweek(272,columns(:,1)) = Demandweek(272,columns(:,1))-columns(:,2)';
        
        Demandweek(145,:) = 0;
        
    end

	% --------------------------------------------------------------------
    % load generator and renewable availability
    % --------------------------------------------------------------------
	
    load('Gendata2010.mat', 'Geninmonth', 'AvailabilityFinal');
    % load('Gendata2010.mat', 'm');

	% --------------------------------------------------------------------
    % define information for creating GDX file
    % --------------------------------------------------------------------
    
    num_gen        = size(Geninmonth,1);   % define the number of generators - 3569
    num_bus        = 312;                  % define number of buses - 312
    num_hours_week = 168;                  % define the number of hours in a week
    num_zones      = max(Geninmonth(:,9)); % define number of zones - 13
    num_econr      = max(Geninmonth(:,8)); % define number of economic regions - 2, 6, 14
    num_ftype      = max(Geninmonth(:,5)); % define number of fuel types - 20
    
	% --------------------------------------------------------------------
    % define function to create identifications for GAMS programs
    % --------------------------------------------------------------------
    
	guel = (@(s,v) strcat(s,strsplit(num2str(v))));

	% --------------------------------------------------------------------
    % define main sets for GAMS program
    % --------------------------------------------------------------------
	
    %       define structure for generator id
    g.name  = 'g';                 % GAMS name
    g.uels  = guel('g',1:num_gen); % GAMS UELs 

    %       define structure for bus id
    i.name  = 'i';
    i.uels  = guel('i',1:num_bus);

    %       define structure for demand hours in a week
    d.name  = 'd'; 
    d.uels  = guel('d',1:num_hours_week);
	
	%       define structure for generator fuel type
    ft.name = 'ft';
    ft.uels = guel('',1:num_ftype);

    %       define structure for zones in WECC
    zones.name = 'zones';
    zones.uels = guel('',1:num_zones);
    
    %       define structure for economic regions
    econr.name = 'econr';
    econr.uels = guel('',1:num_econr);

    %       define structure for economics regions (replace and remove later)
    GenLocation.name = 'GenLocation';
    GenLocation.uels = econr.uels;
                  
    % define structure for generator affected by WBM (g \in WBM_{d})
    gOff.name = 'gOff';
    gOff.uels = guel('g',GeneratorsOff');

    % --------------------------------------------------------------------
    % define crosswalk sets for GAMS program
    % --------------------------------------------------------------------
    
    %       define structure for generator to bus crosswalk (g \in i)
    ig.name     = 'ig';
    ig.val      = [Geninmonth(:,6) Geninmonth(:,1)];        
    ig.uels     = {i.uels, g.uels};         
                % bus id, gen i

    %       define structure for generator to zone crosswalk (g \in zone) (replace and remove later)
    Gentozone.name = 'Gentozone';
    Gentozone.val  = [Geninmonth(:,1) Geninmonth(:,9)]; 
    Gentozone.uels = {g.uels,zones.uels};    
                % g id, zones id
                     
    %       define structure for generator to zone crosswalk (g \in zone)
    gzones.name = 'gzones';
    gzones.val  = [Geninmonth(:,1) Geninmonth(:,9)]; 
    gzones.uels = {g.uels,zones.uels};    
                % g id, zones id   
                     
    % define structure for generator to zone crosswalk (g \in zone)
    geconr.name = 'geconr';
    geconr.val  = [Geninmonth(:,1) Geninmonth(:,8)]; 
    geconr.uels = {g.uels, econr.uels};    
                % g id, econ region id                     
	
    % define structure for generator to bus crosswalk (g \in i)
    gi.name     = 'gi';
    gi.val      = [Geninmonth(:,1) Geninmonth(:,6)];        
    gi.uels     = {g.uels, i.uels};         
                % gen id, bus id - reverse map of ig     
              
    % define structure for generator to bus crosswalk (g \in i)
    gftype.name = 'gftype';
    gftype.val  = [Geninmonth(:,1) Geninmonth(:,5)];        
    gftype.uels = {g.uels, ft.uels};         
                % gen id, bus id                                               

    % --------------------------------------------------------------------
    % define demand and renewable parameters for GAMS program
    % --------------------------------------------------------------------          
              
    % 		define structure for demand at bus i in hour d of week w (demand in 2010)
    p1Demand.name = 'p1Demand';
    p1Demand.type = 'parameter';
    p1Demand.form = 'full';
    p1Demand.val  = Demandweek;
    p1Demand.uels = {guel('i',1:num_bus),guel('d',1:num_hours_week)}; 
                    % bus id, hour id

    % 		define structure for renewable availability in hour d of week w
    Avail.name    = 'Avail';
    Avail.type    = 'parameter';
    Avail.form    = 'full';
    Avail.val     = AvailabilityFinal(:,FirstHour:LastHour);
    Avail.uels    = {guel('g',1:num_gen),guel('d',1:num_hours_week)};
                    % g id, hour id
                    
    % 		define structure for generation outage of generator g from WBM 
    GenLevel.name = 'GenLevel';
    GenLevel.type = 'parameter';
    GenLevel.form = 'full';
    GenLevel.val  = PowerOutagesbyGen(:,FirstHour:LastHour);
    GenLevel.uels = {guel('g',GeneratorsOff'),guel('d',1:num_hours_week)}; 
                       %g id, d id                
	
	% --------------------------------------------------------------------
    % define generator parameters for GAMS program
    % --------------------------------------------------------------------   
	
    % 		define structure for minimum power generation of generator g
    Pmin.name        = 'Pmin';
    Pmin.type        = 'parameter';
    Pmin.form        = 'full';
    Pmin.val         = Geninmonth(:,2);
    Pmin.uels        = {guel('g',1:num_gen)};

    % 		define structure for maximum power generation of generator g
    Pmax.name        = 'Pmax';
    Pmax.type        = 'parameter';
    Pmax.form        = 'full';
    Pmax.val         = Geninmonth(:,3);
    Pmax.uels        = {guel('g',1:num_gen)};

    % 		define structure for production cost (variable) for generator g
    ProdCost.name    = 'a';
    ProdCost.type    = 'parameter';
    ProdCost.form    = 'full';
    ProdCost.val     = Geninmonth(:,4);
    ProdCost.uels    = {guel('g',1:num_gen)};

    % 		define structure for fuel type (1-20) of generator g
    FuelType.name    = 'FuelType';
    FuelType.type    = 'parameter';
    FuelType.form    = 'full';
    FuelType.val     = Geninmonth(:,5); 
    FuelType.uels    = {guel('g',1:num_gen)}; 

    % 		define structure for startup cost of generator g
    Startupcost.name = 'Startupcost';
    Startupcost.type = 'parameter';
    Startupcost.form = 'full';
    Startupcost.val  = Geninmonth(:,7); 
    Startupcost.uels = {guel('g',1:num_gen)};

    % 		define structure for generator by economic region (g \in R)
    GenArea.name     = 'GenArea';
    GenArea.type     = 'parameter';
    GenArea.form     = 'full';
    GenArea.val      = Geninmonth(:,8);
    GenArea.uels     = {guel('g',1:num_gen)};		 

    % 		define structure for heat rate of generator g
    HeatRate.name    = 'HeatRate';
    HeatRate.type    = 'parameter';
    HeatRate.form    = 'full';
    HeatRate.val     = Geninmonth(:,10);
    HeatRate.uels    = {guel('g',1:num_gen)};

    % 		define structure for fuel cost of production (variable) for generator g
    FuelCost.name    = 'FuelCost';
    FuelCost.type    = 'parameter';
    FuelCost.form    = 'full';
    FuelCost.val     = Geninmonth(:,11);
    FuelCost.uels    = {guel('g',1:num_gen)};

    % 		define structure for variable operating and maintenance cost for generator g 
    VOM.name         = 'VOM';
    VOM.type         = 'parameter';
    VOM.form         = 'full';
    VOM.val          = Geninmonth(:,12);
    VOM.uels         = {guel('g',1:num_gen)};

    % 		define structure for ramp rate of generator g
    Ramp.name        = 'Ramp';
    Ramp.type        = 'parameter';
    Ramp.form        = 'full';
    Ramp.val         = Geninmonth(:,13);
    Ramp.uels        = {guel('g',1:num_gen)};

    % 		define structure for uptime requirement of generator g
    Uptime.name      = 'Uptime';
    Uptime.type      = 'parameter';
    Uptime.form      = 'full';
    Uptime.val       = Geninmonth(:,14); 
    Uptime.uels      = {guel('g',1:num_gen)};

    % 		define structure for downtime requirement of generator g
    Downtime.name    = 'Downtime';
    Downtime.type    = 'parameter';
    Downtime.form    = 'full';
    Downtime.val     = Geninmonth(:,15); 
    Downtime.uels    = {guel('g',1:num_gen)};

    % --------------------------------------------------------------------
    % define coupled model parameters for GAMS program
    % -------------------------------------------------------------------- 
    
    % 		define structure for percentage change in demand from REM
    DChange.name    = 'DChange';
    DChange.type    = 'parameter';
    DChange.form    = 'full';
    DChange.val     = DemandMult_UC; 
    DChange.uels    = {guel('',[1,2,3,4,6,7,8,9,12,13])};
    
    % 		define structure for percentage change in demand from REM
    DChange_DC.name = 'DChange_DC';
    DChange_DC.type = 'parameter';
    DChange_DC.form = 'full';
    DChange_DC.val  = DemandMult_DC; 
    DChange_DC.uels = {guel('',1:num_econr-3)};
                   
	% -------------------------------------------------------------------- 
    % post summer code -- WHY? Mort, can we REMOVE?  
    % --------------------------------------------------------------------       
    
    % %if (week >= 26) && (week <= 35)
    %    GenLevel.val=PowerOutagesbyGen(:,FirstHour:LastHour); copied from here
    % %else
    % %   [r,c] = size(PowerOutagesbyGen(:,FirstHour:LastHour));
    % %   GenLevel.val = zeros(r,c);
    % %end
    % 
    % %Genoff.val=GenonoffperHour(:,FirstHour:LastHour);
    % %mkdir(sprintf('week%d',week))

	% --------------------------------------------------------------------
    % define a pause in the paralellization: required for ACI
    % (crashes happen when multiple licenses of GAMS are called at once)
    % --------------------------------------------------------------------
    
	rng('shuffle')
    sleepy_time=rand()*60;
    fprintf('Sleeping for %f seconds\n', sleepy_time)
    pause(sleepy_time);

	% --------------------------------------------------------------------
    % export structures to a single GAMS gdx file for TrialImport.gms
    % --------------------------------------------------------------------   
	
    filename = ['MtoG', '.gdx'];
    wgdx(filename,d,i,g,ft,ig,econr,gzones,geconr,gi,gftype,gOff,GenLevel,GenLocation,Pmin,Pmax,ProdCost,DChange,DChange_DC,FuelType,Avail,p1Demand,Startupcost,GenArea,zones,Gentozone,HeatRate,FuelCost,VOM,Ramp,Uptime,Downtime);

    % --------------------------------------------------------------------
	% post summer code -- WHY? Mort, can we REMOVE?
    % --------------------------------------------------------------------
    % %        if (week >= 26) && (week <= 35)
    % %           copyfile 'TrialImportData_onpeak.gms'  'TrialImportData.gms'
    % %           fprintf('Running Summer Peak Version\n') 
    % %        else
    % %           copyfile 'TrialImportData_offpeak.gms'  'TrialImportData.gms'
    % %           fprintf('Running Non-Summer Offpeak Version\n') 
    % %        end
    % %    
    
    % --------------------------------------------------------------------
    % Call the Unit Commitment program (GAMS)
    % --------------------------------------------------------------------
    myrc = system('gams "UnitCommitment" lo=3'); % GAMS call
    
    % quit Matlab if the GAMS program files
    if (myrc ~= 0)&&(myrc ~= 3)
        disp("GAMS job failure in UnitCommitment, aborting MATLAB job")
        quit(myrc)
    end
    fprintf('GAMS UnitCommitment return code is %d\n', myrc)

    %       import Unit Commitment model outputs
    % UC model status code
    dummygen.name = 'ModelStatus';
    dummygen.form = 'sparse';
    Temp          = rgdx('ModelStatusCodes',dummygen);
    clear dummygen;
    
    z2UC = Temp.val;
    fprintf('GAMS UnitCommitment Model Status code is %d\n', z2UC)

    % UC model solver status
    dummygen.name = 'SolverStatus';
    dummygen.form = 'sparse';
    Temp2         = rgdx('ModelStatusCodes',dummygen);
    clear dummygen;
    
    z1UC = Temp2.val;
    fprintf('GAMS UnitCommitment Solver Status code is %d\n', z1UC)

    % UC model MIP gap constraint
    dummygen.name = 'pMIPGAP';
    dummygen.form = 'sparse';
    Temp2         = rgdx('ModelStatusCodes',dummygen);
    clear dummygen;
    
    mipgap = Temp2.val;
    fprintf('GAMS UnitCommitment MIPGAP is %d\n', mipgap)
    
    % UC model MIP gap percentage
    dummygen.name = 'pGAPPERC';
    dummygen.form = 'sparse';
    Temp2         = rgdx('ModelStatusCodes',dummygen);
    clear dummygen;
    
    gapperc = Temp2.val;
    fprintf('GAMS UnitCommitment GAP percentage is %d\n', gapperc)

    % --------------------------------------------------------------------
    % Call DCOPF model (GAMS)
    % --------------------------------------------------------------------    
    myrc = system('gams "DCOPF" lo=3'); % GAMS call
    
    % quit Matlab if the GAMS program files
    if (myrc ~= 0)&&(myrc ~= 3)
        disp("GAMS job failure in DCOPF, aborting MATLAB job")
        quit(myrc)
    end
    fprintf('GAMS DCOPF return code is %d\n', myrc)

    %       import DCOPF model outputs
    % DCOPF model status code
    dummygen.name = 'ModelStatus';
    dummygen.form = 'sparse';
    Temp          = rgdx('ModelStatusCodesOPF',dummygen);
    clear dummygen;
    
    z2OPF=Temp.val;
    fprintf('GAMS DCOPF Model Status code is %d\n', z2OPF)
    
    % DCOPF model solver status
    dummygen.name = 'SolverStatus';
    dummygen.form = 'sparse';
    Temp2         = rgdx('ModelStatusCodesOPF',dummygen);
    clear dummygen;
    
    z1OPF=Temp2.val;
    fprintf('GAMS DCOPF Solver Status code is %d\n', z1OPF)

	% --------------------------------------------------------------------
    % Repeat both Unit Commitment and DCOPF if model status code of 
	% DCOPF is not equal to 1
    % --------------------------------------------------------------------
   
    if (z2OPF ~= 1)

       % 		Call the Unit Commitment model (GAMS) 
        myrc = system('gams "UnitCommitment" lo=3'); % call GAMS UC model
        
        %       quit Matlab if the GAMS program files
        if (myrc ~= 0)&&(myrc ~= 3)                  
            disp("GAMS job failure in UnitCommitment, aborting MATLAB job")
            quit(myrc)
        end
        fprintf('GAMS UnitCommitment return code is %d\n', myrc)
        
        % 		import Unit Commitment model outputs
        % UC model status
        dummygen.name = 'ModelStatus';
        dummygen.form = 'sparse';
        Temp          = rgdx('ModelStatusCodes',dummygen);
        clear dummygen;
        z2UC          = Temp.val;
        fprintf('GAMS UnitCommitment Model Status code is %d\n', z2UC)

        % UC model solver status
        dummygen.name = 'SolverStatus';
        dummygen.form = 'sparse';
        Temp2         = rgdx('ModelStatusCodes',dummygen);
        clear dummygen;
        z1UC          = Temp2.val;
        fprintf('GAMS UnitCommitment Solver Status code is %d\n', z1UC)

        % UC model MIP gap status
        dummygen.name = 'pMIPGAP';
        dummygen.form = 'sparse';
        Temp2         = rgdx('ModelStatusCodes',dummygen);
        clear dummygen;
        mipgap        = Temp2.val;
        fprintf('GAMS UnitCommitment MIPGAP is %d\n', mipgap)

        % UC model MIP gap percentage
        dummygen.name = 'pGAPPERC';
        dummygen.form = 'sparse';
        Temp2         = rgdx('ModelStatusCodes',dummygen);
        clear dummygen;
        gapperc       = Temp2.val;
        fprintf('GAMS UnitCommitment GAP percentage is %d\n', gapperc)

        % 		Call DCOPF model (GAMS)		
        myrc = system('gams "DCOPF" lo=3'); % call GAMS program
        
        %       quit Matlab if the GAMS program files
        if (myrc ~= 0)&&(myrc ~= 3)
            disp("GAMS job failure in DCOPF, aborting MATLAB job")
            quit(myrc)
        end
        fprintf('GAMS DCOPF return code is %d\n', myrc)

        % 		import DCOPF model outputs		
        % DCOPF model status
        dummygen.name = 'ModelStatus';
        dummygen.form = 'sparse';
        Temp          = rgdx('ModelStatusCodesOPF',dummygen);
        clear dummygen;
        z2OPF         = Temp.val;
        fprintf('GAMS DCOPF Model Status code is %d\n', z2OPF)

        % DCOPF model solver status
        dummygen.name = 'SolverStatus';
        dummygen.form = 'sparse';
        Temp2         = rgdx('ModelStatusCodesOPF',dummygen);
        clear dummygen;  
        z1OPF         = Temp2.val;
        fprintf('GAMS DCOPF Solver Status code is %d\n', z1OPF)
    end

    % --------------------------------------------------------------------
    % define main result output matrix array
    % --------------------------------------------------------------------
    
    z = (1:num_econr)'; % first columns is the economic region id

	% --------------------------------------------------------------------
    % import results.gdx from the DCOPF model
    % --------------------------------------------------------------------
	
	%       non-served electricity (MWH)- economic region econr - week - added 5
	tmpstruct.name = 'nse_r';                       % define parameter name to import
	tmpstruct.form = 'sparse';                      % define form as sparese to get uels
    tmpstruct.uels = {econr.uels};                  % define the uels labels (if present)
	tmpmat         = rgdx('results',tmpstruct);     % import gams parameter
    [zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); % find uels index to merge results
    z(:,end+1)  = NaN;                              % define extra column of NaN at end of the results matrix
    z(zpos,end) = tmpmat.val(zidx(zpos),end)-5;     % merge parameter values to end of the results matrix
    clear tmpstruct tmpmat zpos zidx;               % clear tmp variables
	
	%       total cost of electricity generation ($) - economic region r - week
	tmpstruct.name = 'toc_r';
	tmpstruct.form = 'sparse';                      
    tmpstruct.uels = {econr.uels};                  
	tmpmat         = rgdx('results',tmpstruct);     
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;   
	
    %       electricity generation (MWH) - economic region r - week
    tmpstruct.name = 'gen_r';
    tmpstruct.form = 'sparse';                      
    tmpstruct.uels = {econr.uels};                  
	tmpmat         = rgdx('results',tmpstruct);     
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;  
	
    %       coal electricity generation (MWH) - economic region r - week
    tmpstruct.name = 'coal_gen_r';
    tmpstruct.form = 'sparse';                      
    tmpstruct.uels = {econr.uels};                  
	tmpmat         = rgdx('results',tmpstruct);     
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;  

    %       natural gas (CT, CCGT, and steam turbine) electricity generation (MWH) - economic region r - week
    tmpstruct.name = 'ngas_gen_r';
    tmpstruct.form = 'sparse';                      
    tmpstruct.uels = {econr.uels};                  
	tmpmat         = rgdx('results',tmpstruct);     
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;  
	
	%       electricity demand (MWH) - economic region r - hour d  - added 5
	tmpstruct.name = 'dem_rd';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels,d.uels};
	tmpmat         = rgdx('results',tmpstruct);
	ex_dem_rd        = tmpmat.val;                         % non-main results
	ex_dem_rd(:,end) = tmpmat.val(:,end) - 5;              % subtrate 5 
    wd_dem_rd        = ReshapeMatrix(ex_dem_rd);      % reshape long (r*d x 1) to wide (r x d) 
    clear tmpstruct tmpmat;

	%       shed electricity (MWH) - bus i - hour d  - added 5
	tmpstruct.name = 'shed_ird';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {i.uels,econr.uels,d.uels};
	tmpmat         = rgdx('results',tmpstruct);
	ex_shed_ird        = tmpmat.val;
	ex_shed_ird(:,end) = tmpmat.val(:,end) - 5;    
    wd_shed_ird        = ReshapeMatrix(ex_shed_ird);  % reshape long (i*d x 1) to wide (i \in r x d) 
    clear tmpstruct tmpmat;	

	% --------------------------------------------------------------------
    % import results2.gdx from the DCOPF model
    % --------------------------------------------------------------------
	%The below are inputs are for analysis purpose and have to not be included later on
	
	% electricity generation (MWH) - steam turbine - economic region econr - week
	tmpstruct.name = 'steam_gen_r';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels};
	tmpmat         = rgdx('results2',tmpstruct); 
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;   

	% electricity generation (MWH) - nuclear - economic region econr - week
	tmpstruct.name = 'nuclear_gen_r';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels};
	tmpmat         = rgdx('results2',tmpstruct); 
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;  

	% electricity generation (MWH) - bio-gass - economic region econr - week
	tmpstruct.name = 'biogas_gen_r';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels};
	tmpmat         = rgdx('results2',tmpstruct); 
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;  

	% electricity generation (MWH) - geothermal - economic region econr - week
	tmpstruct.name = 'geo_gen_r';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels};
	tmpmat         = rgdx('results2',tmpstruct); 
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;    

	% electricity generation (MWH) - petroleum - economic region econr - week
	tmpstruct.name = 'petrol_gen_r';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels};
	tmpmat         = rgdx('results2',tmpstruct); 
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;   

	% electricity generation (MWH) - itertie - economic region econr - week
	tmpstruct.name = 'intertiel_gen_r';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels};
	tmpmat         = rgdx('results2',tmpstruct); 
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;   

	% electricity generation (MWH) - wind (onshore) - economic region econr - week
	tmpstruct.name = 'wind_gen_r';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels};
	tmpmat         = rgdx('results2',tmpstruct); 
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx;  

	% 'non-served electricity (MWH) - economic region econr - hour d  - added 5'
	tmpstruct.name = 'nse_rd';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels,d.uels};
	tmpmat         = rgdx('results2',tmpstruct);   	
	ex_nse_rd        = tmpmat.val;
	ex_nse_rd(:,end) = tmpmat.val(:,end) - 5;
    wd_nse_rd        = ReshapeMatrix(ex_nse_rd);      % reshape long (r*d x 1) to wide (r x d)
	clear tmpstruct tmpmat; 
	
	% --------------------------------------------------------------------
    % import results3.gdx from the DCOPF model
    % --------------------------------------------------------------------
	
	%       electricity generation (MWH) - hydro - economic region econr - week
	tmpstruct.name = 'hydro_gen_r';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels};
	tmpmat         = rgdx('results3',tmpstruct); 
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx; 

	%       electricity generation (MWH) - ps hydro - economic region econr - week
	tmpstruct.name = 'pshydro_gen_r';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels};
	tmpmat         = rgdx('results3',tmpstruct); 
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx; 

	%       electricity generation (MWH) - motorload - economic region econr - week
	tmpstruct.name = 'motorload_gen_r';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels};
	tmpmat         = rgdx('results3',tmpstruct); 
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx; 

	%       electricity generation (MWH) - solar (all) - economic region econr - week
	tmpstruct.name = 'solar_gen_r';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels};
	tmpmat         = rgdx('results3',tmpstruct); 
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx; 
	
	% --------------------------------------------------------------------
    % import MarginalElectricityPrice.gdx from the DCOPF model
    % --------------------------------------------------------------------
	
	%       total cost of electricity generation ($) - coal - economic region r - week
	tmpstruct.name  = 'toc_coal_r';
	tmpstruct.form  = 'sparse';
	tmpstruct.uels  = {econr.uels};
	tmpmat          = rgdx('MarginalElectricityPrice',tmpstruct);
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx; 

	%       total cost of electricity generation ($) - natural gas (CT, CCGT) - economic region r - week
	tmpstruct.name  = 'toc_ngas_r';
	tmpstruct.form  = 'sparse';
	tmpstruct.uels  = {econr.uels};
	tmpmat          = rgdx('MarginalElectricityPrice',tmpstruct);
	[zpos,zidx] = ismember(z(:,1),tmpmat.val(:,1)); 
	z(:,end+1)  = NaN;                              
	z(zpos,end) = tmpmat.val(zidx(zpos),end);     
	clear tmpstruct tmpmat zpos zidx; 
		
	%       electricity generation (MWH) - economic region r - hour d   - added 5
	tmpstruct.name = 'mc_ird';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {i.uels,econr.uels,d.uels};
	tmpmat         = rgdx('MarginalElectricityPrice',tmpstruct);
	ex_mc_ird        = tmpmat.val;
	ex_mc_ird(:,end) = tmpmat.val(:,end) - 5;
    wd_mc_ird        = ReshapeMatrix(ex_mc_ird);      % reshape long (i*d x 1) to wide (i \in r x d)
	clear tmpstruct tmpmat;  

	%       marginal cost of electricity generation ($\MWH) - bus i in economic region r - hour d  - added 5
	tmpstruct.name = 'mdem_ird';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {i.uels,econr.uels,d.uels};
	tmpmat         = rgdx('MarginalElectricityPrice',tmpstruct);
	clear tmpstruct;	 
	ex_mdem_ird        = tmpmat.val;
	ex_mdem_ird(:,end) = tmpmat.val(:,end) - 5;
    wd_mdem_ird        = ReshapeMatrix(ex_mdem_ird);  % reshape long (i*d x 1) to wide (r x d)
	clear tmpstruct tmpmat;	

	%       non-served electricity (MWH) - bus i - hour d  - added 5
	tmpstruct.name = 'nse_ird';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {i.uels,econr.uels,d.uels};
	tmpmat = rgdx('results',tmpstruct);
	clear tmpstruct;	
	ex_nse_ird        = tmpmat.val;
	ex_nse_ird(:,end) = tmpmat.val(:,end) - 5;
    wd_nse_ird        = ReshapeMatrix(ex_nse_ird);    % reshape long (i*d x 1) to wide (r x d)
	clear tmpstruct tmpmat;	

	%       electricity generation (MWH) - economic region r - hour d   - added 5
	tmpstruct.name = 'gen_rd';
	tmpstruct.form = 'sparse';
	tmpstruct.uels = {econr.uels,d.uels};
	tmpmat         = rgdx('MarginalElectricityPrice',tmpstruct);    
	ex_gen_rd = tmpmat.val;
    wd_gen_rd = ReshapeMatrix(ex_gen_rd);             % reshape long (r*d x 1) to wide (r x d)
	clear tmpstruct tmpmat;

	% --------------------------------------------------------------------
    % slow program for process purposes
    % --------------------------------------------------------------------
	
	disp("Calculations complete, now about to write files")
	fprintf('Sleeping for %f seconds\n', sleepy_time)
	pause(sleepy_time);

	% --------------------------------------------------------------------
    % export main results from results,2,3 to week w csv file
    % --------------------------------------------------------------------
	
	filename = ['week', int2str(week), '.csv'];
	csvwrite(filename,z);

	% --------------------------------------------------------------------
    % export results from  week w to csv files
    % --------------------------------------------------------------------

	% 		csv file for marginal cost  - bus i in economic region r - hour d
	filename = ['MarginalCost', int2str(week), '.csv'];
	csvwrite(filename,wd_mc_ird);

	
	% 		csv file for electricity demand - bus i in economic region r - hour d
	filename = ['MultipliedDemand', int2str(week), '.csv'];
	csvwrite(filename,wd_mdem_ird);

	% 		csv file for non-served electricity at bus i for hour d
	filename = ['UnservedEnergywhole', int2str(week), '.csv'];
	csvwrite(filename,wd_nse_ird);

	% 		csv file for electricity generation - economic region r  - hour d
	filename = ['Genperhour', int2str(week), '.csv'];
	csvwrite(filename,wd_gen_rd);

	% 		csv file for electricity demand - economic region r - hour d
	filename = ['NewDemandwhole', int2str(week), '.csv'];
	csvwrite(filename,wd_dem_rd);

	% 		csv file for non-served electricity - economic region r - hour d
	filename = ['NewEnergywholeunserved', int2str(week), '.csv'];
	csvwrite(filename,wd_nse_rd);

	% 		csv file for shed energy - bus i - hour d
	filename = ['ShedEnergywhole', int2str(week), '.csv'];
	csvwrite(filename,wd_shed_ird);

	% --------------------------------------------------------------------
    % export unit commitment and DCOPF model status results
    % --------------------------------------------------------------------

	% 		output ACI job id -- WHY?
	myjobid  = getenv('PBS_JOBID');
	filename = ['CheckJob', int2str(week), '.csv'];
	csvwrite(filename,myjobid);

	% 		csv file for all unit commitment and DCOPF model and solver status information
	filename   = ['SolverStatus', int2str(week), '.csv'];
	modelstats = [z2UC, z1UC, mipgap, z2OPF, z1OPF, gapperc]; 
	csvwrite(filename,modelstats);

	% --------------------------------------------------------------------
    % export unit commitment and DCOPF model status results
    % --------------------------------------------------------------------

	disp("File write operations are complete, returning to control program")
	
	% --------------------------------------------------------------------
    % end timer
    % --------------------------------------------------------------------
	
	toc  
	
end
