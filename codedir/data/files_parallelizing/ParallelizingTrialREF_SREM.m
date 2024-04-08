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
% state - states 1 - 14 -- we might change this later
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
    oldFolder = cd(newFolder)
    fprintf('Switching from %80s to %80s\n', oldFolder,newFolder)

    % --------------------------------------------------------------------
    % remove when on ACI this is the test a single week
    % --------------------------------------------------------------------
	
%	 cd Gams_data
%    week = 28;
	
	% --------------------------------------------------------------------
    % begin timer
    % --------------------------------------------------------------------
	
    tic
	
	% --------------------------------------------------------------------
    % load percentage change in demand from the economic model (intially
    % set to zeros
    % --------------------------------------------------------------------
    
	fid = fopen('Input_data.csv');
	Input_data = textscan(fid, '%f%f','headerLines',1,'delimiter',',','collectoutput',1);
	fclose(fid);
	
	%		define demand multipler from REM
	DemandMult       = Input_data{1,1};

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
	
%	% delete extra day from base WECC demand
%	demand(:,8761:end) = [];
%
%	% Load the demand response factors for this case
%	load('load_scale_factors.mat','busdata');
%	
%	% modify demand(bus,hour) by DR factor
%	demand = demand .* busdata;
	
	% --------------------------------------------------------------------
    % define the start and final hour of the week
    % --------------------------------------------------------------------
    
	FirstHour  = 1+(week-1)*168;
    LastHour   = week*168;
	
    % 		define demand by bus for the specific week
	
    Demandweek = demand(:,FirstHour:LastHour);

	% --------------------------------------------------------------------
    % alter week 48 -- WHY?
    % --------------------------------------------------------------------
    if week==48
        for ii=2:1:168
        if Demandweek(199,ii)>9000
        Demandweek(199,ii) = Demandweek(199,ii-1);
        end
        end
        Demandweek(202,:) = Demandweek(202,:)-200;
        columns=[152 345; 153 150; 154 27; 161 94; 155 52; 156 110; 162 182; 163 289];
        Demandweek(98,columns(:,1)) = Demandweek(98,columns(:,1))-columns(:,2)';
        columns=[154 51 161 60];
        Demandweek(272,columns(:,1)) = Demandweek(272,columns(:,1))-columns(:,2)';
        Demandweek(145,:) = 0;
    end

	% --------------------------------------------------------------------
    % SECTION IS NOT USED -- REMOVE?
    % --------------------------------------------------------------------
    %This is the availability Data for the renewable generators.
    % fid=fopen('Availshapes.csv');
    % Availabilitytemp=textscan(fid, [repmat('%f',[1,169])],'headerlines',1,'delimiter',',','collectoutput',1);
    % fclose(fid);
    % Availability = Availabilitytemp{1,1};
    % [m1,n1]=size(Availability);
    % Availability=max(Availability,0);

    % fid=fopen('AvailTrial2.csv');
    % Availabilitytemp=textscan(fid, [repmat('%f',[1,8761])],'headerlines',1,'TreatAsEmpty',{'#N/A','na'},'delimiter',',','collectoutput',1);
    % fclose(fid);
    % Availability = Availabilitytemp{1,1};
    % Availability(isnan(Availability)) = 10000;
    % [m1,n1]=size(Availability);
    % Availability=max(Availability,0);
    % Availability=min(Availability,10000);
    % 
    % 
    % filename = 'Gendatawithyears.csv';
    % delimiter = ',';
    % formatSpec = '%f%s%s%f%f%f%f%f%f%f%f%f%f%f%f%f%f';
    % fileID = fopen(filename,'r');
    % dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines',1, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    % fclose(fileID);
    % 
    % GenRetire=dataArray{:,3};
    % GenCommission=dataArray{:,2};
    % GenNumber=dataArray{:,1};
    % MinPower=dataArray{:,4};
    % MaxPower=dataArray{:,5};
    % ProdCostvalues=dataArray{:,6};
    % FuelTypevalues=dataArray{:,7};
    % BusID=dataArray{:,8};
    % StartUpC=dataArray{:,9};
    % AreaCodeCAorRO=dataArray{:,10};
    % StateZone=dataArray{:,11};
    % HeatRateAll=dataArray{:,12};
    % FuelCostAll=dataArray{:,13};
    % VariableCostAll=dataArray{:,14};
    % RampRateAll=dataArray{:,15};
    % UptimeAll=dataArray{:,16};
    % DownTimeAll=dataArray{:,17};
    % 
    % %This is the particular year to analyze. In this case we have taken it for
    % %2014 year.
    % DateStringInitial = '31-Dec-2009';
    % DateTimeInitial=datetime(DateStringInitial);
    % formatIn = 'dd-mmm-yyyy';
    % row1=datenum(DateStringInitial,formatIn);
    % 
    % DateStringFinal = '31-Dec-2010';
    % formatIn = 'dd-mmm-yyyy';
    % row2=datenum(DateStringFinal,formatIn);
    % 
    % 
    % 
    % l=1;
    % [m,n]=size(GenRetire);
    %     %Here we check for he condition of the Retirements and Commisioning.
    %     for ii=1:1:m
    %         condition1=(datenum(GenRetire{ii})>datenum(DateStringFinal));
    %         condition2=(datenum(GenCommission{ii})<datenum(DateStringFinal));
    %         if (condition1 && condition2)
    %             Geninmonth(l,1) = GenNumber(ii,1);
    %             Geninmonth(l,2) = MinPower(ii,1);
    %             Geninmonth(l,3) = MaxPower(ii,1);
    %             Geninmonth(l,4) = ProdCostvalues(ii,1);
    %             Geninmonth(l,5) = FuelTypevalues(ii,1);
    %             Geninmonth(l,6) = BusID(ii,1);
    %             Geninmonth(l,7) = StartUpC(ii,1);
    %             Geninmonth(l,8) = AreaCodeCAorRO(ii,1);
    %             Geninmonth(l,9) = StateZone(ii,1);
    %             Geninmonth(l,10) = HeatRateAll(ii,1);
    %             Geninmonth(l,11) = FuelCostAll(ii,1);
    %             Geninmonth(l,12) = VariableCostAll(ii,1);
    %             Geninmonth(l,13) = RampRateAll(ii,1);
    %             Geninmonth(l,14) = UptimeAll(ii,1); 
    %             Geninmonth(l,15) = DownTimeAll(ii,1); 
    %             l=l+1;
    %         end
    %     end
    %     
    %     GeninmonthTemp=Geninmonth;
    %     [m,n]=size(Geninmonth);
    %    %Below is the availability data for that particular year.
    %     k=1;
    %     for ii=1:1:m1
    %         for j=1:1:m
    %         if Availability(ii,1) = =GeninmonthTemp(j,1)
    %           AvailabilityFinal(k,:) = Availability(ii,:);
    %           k=k+1;
    %         end
    %         end
    %     end
    %     
    %     AvailabilityFinal(:,1) = [];
    %     Geninmonth(:,1) = linspace(1,m,m); %This is to rename the generators to help facilitate it in GAMS
    %     
    %     save('Gendata2010.mat', 'AvailabilityFinal','Geninmonth','m'); 

	% --------------------------------------------------------------------
    % load generator and renewable availability
    % --------------------------------------------------------------------
	
    load('Gendata2010.mat', 'Geninmonth', 'AvailabilityFinal');
    % load('Gendata2010.mat', 'm');

	% --------------------------------------------------------------------
    % define information for creating GDX file
    % --------------------------------------------------------------------
	
    [num_gen,~]=size(Geninmonth);      % define the number of generators - 3569
    num_bus = 312;                     % define number of buses - 312
    num_hours_week = 168;              % define the number of hours in a week
    num_zone = max(Geninmonth(:,9));   % define number of zones - 13
    num_econ = max(Geninmonth(:,8));   % define number of economic regions - 14 (might change to 11)

    Geninmonth(:,1) = linspace(1,num_gen,num_gen); % --WHY? -- REMOVE?

	% --------------------------------------------------------------------
    % define function to create identifications for GAMS programs
    % --------------------------------------------------------------------
    
	guel = (@(s,v) strcat(s,strsplit(num2str(v))));

	% --------------------------------------------------------------------
    % define sets for GAMS program
    % --------------------------------------------------------------------
	
    % define structure for generator id
    g.name = 'g';                 % GAMS name
    g.uels = guel('g',1:num_gen); % GAMS UELs 

    % define structure for bus id
    i.name = 'i';
    i.uels = guel('i',1:num_bus);

    % define structure for demand hours in a week
    d.name = 'd'; 
    d.uels = guel('d',1:num_hours_week);

    % define structure for zones in WECC
    zones.name = 'zones';
    zones.uels = guel('',1:num_zone);

    % define structure for economics regions (this is based on DREM)
    GenLocation.name = 'GenLocation';
    GenLocation.uels = guel('',1:num_econ);

    % define structure for generator to bus crosswalk (g \in i)
    ig.name = 'ig';
    ig.val  = [Geninmonth(:,6) Geninmonth(:,1)];        
    ig.uels = {i.uels, g.uels};         
              % bus id, gen i

    % define structure for generator to zone crosswalk (g \in zone)
    Gentozone.name = 'Gentozone';
    Gentozone.val  = [Geninmonth(:,1) Geninmonth(:,9)]; 
    Gentozone.uels = {g.uels,zones.uels};    
                     % g id, zones id
					 
    % define structure for generator affected by WBM (g \in WBM_{d})
    gOff.name = 'gOff';
    gOff.uels = guel('g',GeneratorsOff');

    % define structure for demand at bus i in hour d of week w (demand in 2010)
    p1Demand.name = 'p1Demand';
    p1Demand.type = 'parameter';
    p1Demand.form = 'full';
    p1Demand.val  = Demandweek;
    p1Demand.uels = {guel('i',1:num_bus),guel('d',1:num_hours_week)}; 
                    % bus id, hour id

    % define structure for renewable availability in hour d of week w
    Avail.name = 'Avail';
    Avail.type = 'parameter';
    Avail.form = 'full';
    Avail.val  = AvailabilityFinal(:,FirstHour:LastHour);
    Avail.uels = {guel('g',1:num_gen),guel('d',1:num_hours_week)};
                 % g id, hour id
				 
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

    % 		define structure for percentage change in demand from REM
    DChange.name     = 'DChange';
    DChange.type     = 'parameter';
    DChange.form     = 'full';
    DChange.val      = DemandMult; 
    DChange.uels     = {guel('',1:num_econ)};

    % 		define structure for generation outage of generator g from WBM 
    GenLevel.name    = 'GenLevel';
    GenLevel.type    = 'parameter';
    GenLevel.form    = 'full';
    GenLevel.val     = PowerOutagesbyGen(:,FirstHour:LastHour);
    GenLevel.uels    = {guel('g',GeneratorsOff'),guel('d',1:num_hours_week)}; 
                       %g id, d id
    
	% -------------------------------------------------------------------- 
    % REMOVE?  
    % --------------------------------------------------------------------       
    % p1Demand.val = Demandweek;                         %copied from here
    % Avail.val=AvailabilityFinal(:,FirstHour:LastHour); %copied from here
    % 
    % %if (week >= 26) && (week <= 35)
    %    GenLevel.val=PowerOutagesbyGen(:,FirstHour:LastHour);
    % %else
    % %   [r,c] = size(PowerOutagesbyGen(:,FirstHour:LastHour));
    % %   GenLevel.val = zeros(r,c);
    % %end
    % 
    % %Genoff.val=GenonoffperHour(:,FirstHour:LastHour);
    % %mkdir(sprintf('week%d',week))

	% --------------------------------------------------------------------
    % define a pause in the paralellization -- WHY?
    % --------------------------------------------------------------------
    
	rng('shuffle')
    sleepy_time=rand()*60;
    fprintf('Sleeping for %f seconds\n', sleepy_time)
    pause(sleepy_time);

	% --------------------------------------------------------------------
    % export structures to a single GAMS gdx file
    % --------------------------------------------------------------------   
	
    filename = ['MtoG', '.gdx'];
    wgdx(filename,d,i,g,ig,gOff,GenLevel,GenLocation,Pmin,Pmax,ProdCost,DChange,FuelType,Avail,p1Demand,Startupcost,GenArea,zones,Gentozone,HeatRate,FuelCost,VOM,Ramp,Uptime,Downtime);

    % --------------------------------------------------------------------
	% REMOVE?
    % --------------------------------------------------------------------
    % %        if (week >= 26) && (week <= 35)
    % %           copyfile 'TrialImportData_onpeak.gms'  'TrialImportData.gms'
    % %           fprintf('Running Summer Peak Version\n') 
    % %        else
    % %           copyfile 'TrialImportData_offpeak.gms'  'TrialImportData.gms'
    % %           fprintf('Running Non-Summer Offpeak Version\n') 
    % %        end
    % %        
    % %z2=10;
    % %while  z2 > 8
    % %         system 'gams "UnitCommitment" lo=3';
    % --------------------------------------------------------------------
	
    % --------------------------------------------------------------------
    % Call the Unit Commitment model (GAMS)
    % --------------------------------------------------------------------
    myrc = system('gams "UnitCommitment" lo=3'); % GAMS call
    
    % quit Matlab if the GAMS program files
    if (myrc ~= 0)&&(myrc ~= 3)
        disp("GAMS job failure in UnitCommitment, aborting MATLAB job")
        quit(myrc)
    end
    fprintf('GAMS UnitCommitment return code is %d\n', myrc)

    % import Unit Commitment model outputs
    % --------------------------------------------------------------------
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
    %end 

    %system 'gams "UnitCommitment" lo=3';
    %z2=10;
    %while  z2 > 8
    %         system 'gams "DCOPF" lo=3';
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

    % import DCOPF model outputs
    % --------------------------------------------------------------------
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

    %end 
    
 
	% --------------------------------------------------------------------
    % Repeat both Unit Commitment and DCOPF if model status code of 
	% DCOPF is not equal to 1
    % --------------------------------------------------------------------
    
%    z2OPF = 0;
    if (z2OPF ~= 1)

       % 		Call the Unit Commitment model (GAMS)
	   
        myrc = system('gams "UnitCommitment" lo=3'); % call GAMS UC model
        
        % quit Matlab if the GAMS program files
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
        
        % quit Matlab if the GAMS program files
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
    % import results.gdx from the DCOPF model
    % --------------------------------------------------------------------
	
	% 		non-served electricity - total - annual
	dummygen.name  = 'UnservedEnergy';
	dummygen.form  = 'sparse';
	UnservedEnergy = rgdx('results',dummygen);
	clear dummygen;	
	z(1) = UnservedEnergy.val;
	% fprintf('Unserved Energy is %f\n', z(1));   

	% 		total cost of electricity generation in CA - total - annual
	dummygen.name = 'TotalCostCali';
	dummygen.form = 'sparse';
	TotalcostCali = rgdx('results',dummygen);
	clear dummygen;	 
	z(2) = TotalcostCali.val;

	% 		total cost of electricity generation in ROWECC - total - annual
	dummygen.name = 'TotalCostRO';
	dummygen.form = 'sparse';
	TotalcostRO   = rgdx('results',dummygen);
	clear dummygen;	
	z(3) = TotalcostRO.val;

	% 		electricity generation in CA - total - annual
	dummygen.name ='GenCali';
	dummygen.form ='sparse';
	GenCali       = rgdx('results',dummygen);
	clear dummygen;	
	z(4) = GenCali.val;

	% 		electricity generation in ROWECC - total - annual
	dummygen.name = 'GenRO';
	dummygen.form = 'sparse';
	GenRO         = rgdx('results',dummygen);
	clear dummygen;
	z(5) = GenRO.val;

	% 		coal electricity generation in CA - total - annual
	dummygen.name = 'CaliCoal';
	dummygen.form = 'sparse';
	CaliCoal      = rgdx('results',dummygen);
	clear dummygen;
	z(6) = CaliCoal.val;

	% 		natural gas CT, CCGT, steam, turbine electricity generation in CA  - total - annual
	dummygen.name = 'CaliGas';
	dummygen.form = 'sparse';
	CaliGas       = rgdx('results',dummygen);
	clear dummygen;
	z(7) = CaliGas.val;

	% 		coal electricity generation in ROWECC - total - annual
	dummygen.name = 'RoweccCoal';
	dummygen.form = 'sparse';
	RoweccCoal    = rgdx('results',dummygen);
	clear dummygen;
	z(8) = RoweccCoal.val;
	
	% 		natural gas CT, CCGT, steam, turbine electricity generation in ROWECC  - total - annual
	dummygen.name = 'RoweccGas';
	dummygen.form = 'sparse';
	RoweccGas     = rgdx('results',dummygen);
	clear dummygen;
	z(9) = RoweccGas.val;
	
	% 		electricity demand in CA - total for hour d
	dummygen.name = 'DemandCaliperHour';
	dummygen.form = 'sparse';
	dummygen.uels = {guel('d',1:168)};
	DemandCaliNew = rgdx('results',dummygen);
	clear dummygen;
	DemandwholeCali      = DemandCaliNew.val;
	DemandwholeCali(:,2) = DemandwholeCali(:,2)-5;
	NewDemandwholeCali   = DemandwholeCali(:,2);

	% 		electricity demand in ROWECC - total for hour d 
	dummygen.name = 'DemandROperHour';
	dummygen.form = 'sparse';
	dummygen.uels = {guel('d',1:168)};
	DemandRONew   = rgdx('results',dummygen);
	clear dummygen;
	DemandwholeRO      = DemandRONew.val;
	DemandwholeRO(:,2) = DemandwholeRO(:,2)-5;
	NewDemandwholeRO   = DemandwholeRO(:,2);

	% 		shed energy at bus i for hour d
	dummygen.name = 'ShedEnergyValue';
	dummygen.form = 'sparse';
	dummygen.uels = {guel('i',1:312),guel('d',1:168)};
	ShedEnergywhole = rgdx('results',dummygen);
	clear dummygen;
	EnergywholeShed      = ShedEnergywhole.val;
	EnergywholeShed(:,3) = EnergywholeShed(:,3)-5;
	NewEnergywholeShed   = ReshapeMatrix(EnergywholeShed);

	% --------------------------------------------------------------------
    % import results2.gdx from the DCOPF model
    % --------------------------------------------------------------------
	%The below are inputs are for analysis purpose and have to not be included later on
	
	% 		steam turbine electricity generation in CA - total - annual
	dummygen.name = 'CaliSteam';
	dummygen.form = 'sparse';
	CaliSteam     = rgdx('results2',dummygen);
	clear dummygen;
	z(10) = CaliSteam.val;

	% 		steam turbine electricity generation in ROWECC - total - annual
	dummygen.name = 'RoweccSteam';
	dummygen.form = 'sparse';
	RoweccSteam   = rgdx('results2',dummygen);
	clear dummygen;
	z(11) = RoweccSteam.val;

	% 		nuclear electricity generation in CA - total - annual
	dummygen.name = 'CaliNuclear';
	dummygen.form = 'sparse';
	CaliNuclear   = rgdx('results2',dummygen);
	clear dummygen;
	z(12) = CaliNuclear.val;

	% 		nuclear electricity generation in ROWECC - total - annual
	dummygen.name = 'RoweccNuclear';
	dummygen.form = 'sparse';
	RoweccNuclear = rgdx('results2',dummygen);
	clear dummygen;
	z(13) = RoweccNuclear.val;

	% 		biogas-other electricity generation in CA - total - annual
	dummygen.name = 'CaliBiogas';
	dummygen.form = 'sparse';
	CaliBiogas    = rgdx('results2',dummygen);
	clear dummygen;
	z(14) = CaliBiogas.val;

	% 		biogas-other electricity generation in ROWECC - total - annual
	dummygen.name = 'RoweccBiogas';
	dummygen.form = 'sparse';
	RoweccBiogas  = rgdx('results2',dummygen);
	clear dummygen;
	z(15) = RoweccBiogas.val;

	% 		geothermal electricity generation in CA - total - annual
	dummygen.name = 'CaliGeo';
	dummygen.form = 'sparse';
	CaliGeo       = rgdx('results2',dummygen);
	clear dummygen;
	z(16) = CaliGeo.val;

	% 		geothermal electricity generation in ROWECC - total - annual
	dummygen.name = 'RoweccGeo';
	dummygen.form = 'sparse';
	RoweccGeo     = rgdx('results2',dummygen);
	clear dummygen;
	z(17) = RoweccGeo.val;

	% 		petroleum electricity generation in CA - total - annual
	dummygen.name = 'CaliPetrol';
	dummygen.form = 'sparse';
	CaliPetrol    = rgdx('results2',dummygen);
	clear dummygen;
	z(18) = CaliPetrol.val;

	% 		petroleum electricity generation in ROWECC - total - annual
	dummygen.name = 'RoweccPetrol';
	dummygen.form = 'sparse';
	RoweccPetrol  = rgdx('results2',dummygen);
	clear dummygen;
	z(19) = RoweccPetrol.val;

	% 		itertie electricity generation in CA - total - annual
	dummygen.name = 'CaliIntertie';
	dummygen.form = 'sparse';
	CaliIntertie  = rgdx('results2',dummygen);
	clear dummygen;
	z(20) = CaliIntertie.val;


	% 		itertie electricity generation in ROWECC - total - annual
	dummygen.name = 'RoweccIntertie';
	dummygen.form = 'sparse';
	RoweccIntertie = rgdx('results2',dummygen);
	clear dummygen;
	z(21) = RoweccIntertie.val;


	% 		wind-onshore electricity generation in CA - total - annual
	dummygen.name = 'CaliWind';
	dummygen.form = 'sparse';
	CaliWind      = rgdx('results2',dummygen);
	clear dummygen;
	z(22) = CaliWind.val;

	% 		wind-onshore electricity generation in ROWECC - total - annual
	dummygen.name = 'RoweccWind';
	dummygen.form = 'sparse';
	RoweccWind    = rgdx('results2',dummygen);
	clear dummygen;
	z(23) = RoweccWind.val;

	% 		non-served electricity in CA for hour d 
	dummygen.name = 'UnservedEnergyCaliperHour';
	dummygen.form = 'sparse';
	dummygen.uels = {guel('d',1:168)};
	UnservedEnergyCaliNew  = rgdx('results2',dummygen);
	clear dummygen;
	EnergywholeunservedCali      = UnservedEnergyCaliNew.val;
	EnergywholeunservedCali(:,2) = EnergywholeunservedCali(:,2)-5;
	NewEnergywholeunservedCali   = EnergywholeunservedCali(:,2);

	% 		non-served electricity in ROWECC for hour d
	dummygen.name = 'UnservedEnergyROperHour';
	dummygen.form = 'sparse';
	dummygen.uels = {guel('d',1:168)};
	UnservedEnergyRONew  = rgdx('results2',dummygen);
	clear dummygen;
	EnergywholeunservedRO      = UnservedEnergyRONew.val;
	EnergywholeunservedRO(:,2) = EnergywholeunservedRO(:,2)-5;
	NewEnergywholeunservedRO   = EnergywholeunservedRO(:,2);
	
	% --------------------------------------------------------------------
    % import results3.gdx from the DCOPF model
    % --------------------------------------------------------------------
	
	% 		hydro electricity generation in CA - total - annual
	dummygen.name = 'CaliHydro';
	dummygen.form = 'sparse';
	CaliHydro     = rgdx('results3',dummygen);
	clear dummygen;
	z(24) = CaliHydro.val;

	% 		hydro electricity generation in ROWECC - total - annual
	dummygen.name = 'RoweccHydro';
	dummygen.form = 'sparse';
	RoweccHydro   = rgdx('results3',dummygen);
	clear dummygen;
	z(25) = RoweccHydro.val;

	% 		pump\storage hydro electricity generation in CA - total - annual
	dummygen.name = 'CaliPSHydro';
	dummygen.form = 'sparse';
	CaliPSHydro   = rgdx('results3',dummygen);
	clear dummygen;
	z(26) = CaliPSHydro.val;

	% 		pump\storage hydro electricity generation in ROWECC - total - annual
	dummygen.name = 'RoweccPSHydro';
	dummygen.form = 'sparse';
	RoweccPSHydro = rgdx('results3',dummygen);
	clear dummygen;
	z(27) = RoweccPSHydro.val;

	% 		motor load electricity generation in CA - total - annual
	dummygen.name = 'CaliMotorLoad';
	dummygen.form = 'sparse';
	CaliMotorLoad = rgdx('results3',dummygen);
	clear dummygen;
	z(28) = CaliMotorLoad.val;

	% 		motor load electricity generation in ROWECC - total - annual
	dummygen.name   = 'RoweccMotorLoad';
	dummygen.form   = 'sparse';
	RoweccMotorLoad = rgdx('results3',dummygen);
	clear dummygen;
	z(29) = RoweccMotorLoad.val;

	% 		solar (three types) electricity generation in CA - total - annual
	dummygen.name = 'CaliSolar';
	dummygen.form = 'sparse';
	CaliSolar     = rgdx('results3',dummygen);
	clear dummygen;
	z(30) = CaliSolar.val;


	% 		solar (three types) electricity generation in ROWECC - total - annual
	dummygen.name = 'RoweccSolar';
	dummygen.form = 'sparse';
	RoweccSolar   = rgdx('results3',dummygen);
	clear dummygen;
	z(31) = RoweccSolar.val;

	% ------------------- REMOVE?
	%           dummygen.name = 'GenCoalCali';
	%           dummygen.form = 'sparse';
	%           %VOM.type='parameter';
	%           dummygen.uels = {guel('g',1:m),guel('d',1:168)};
	%           GenCoalCali  = rgdx('results3',dummygen);
	%           CoalCAGen=GenCoalCali.val;
	% -------------------
	
	% 		coal electricity generation in CA for generator g	
	dummygen.name     = 'CoalGenCaliPerGen';
	dummygen.form     = 'sparse';
	dummygen.uels     = {guel('g',1:num_gen)};
	Califroniacoalgen = rgdx('results3',dummygen);
	clear dummygen;
	CoalGenCaliPerGen = (Califroniacoalgen.val);
	

	% 		natural gas CT, CCGT, steam, turbine electricity generation in CA for generator g
	dummygen.name = 'GasGenCaliPerGen';
	dummygen.form = 'sparse';
	dummygen.uels = {guel('g',1:num_gen)};
	CalifroniaGasgen  = rgdx('results3',dummygen);
	clear dummygen;
	GasGenCaliPerGen=(CalifroniaGasgen.val);
	

	% 		coal electricity generation in RO for generator g
	dummygen.name = 'CoalGenROPerGen';
	dummygen.form = 'sparse';
	dummygen.uels = {guel('g',1:num_gen)};
	ROcoalgen     = rgdx('results3',dummygen);
	clear dummygen;
	CoalGenROPerGen=(ROcoalgen.val);
	

	% 		natural gas CT, CCGT, steam, turbine electricity generation in CA for generator g
	dummygen.name = 'GasGenROPerGen';
	dummygen.form = 'sparse';
	dummygen.uels = {guel('g',1:num_gen)};
	ROGasgen      = rgdx('results3',dummygen);
	clear dummygen;
	GasGenROPerGen=(ROGasgen.val);


	% 		non-served electricity in CA - total - annual  - added 5
	dummygen.name      = 'UnservedEnergyCali';
	dummygen.form      = 'sparse';
	CaliUnservedEnergy = rgdx('results3',dummygen);
	clear dummygen;
	z(36) = CaliUnservedEnergy.val-5;


	% 		non-served electricity in RO - total - annual  - added 5
	dummygen.name    = 'UnservedEnergyRO';
	dummygen.form    = 'sparse';
	ROUnservedEnergy = rgdx('results3',dummygen);
	clear dummygen;
	z(37) = ROUnservedEnergy.val-5;

	% --------------------------------------------------------------------
    % import MarginalElectricityPrice.gdx from the DCOPF model
    % --------------------------------------------------------------------
    
	% 		marginal cost in CA of bus i for hour d  - added 5
	dummygen.name    = 'MarginalCostCali';
	dummygen.form    = 'sparse';
	dummygen.uels    = {guel('i',1:312),guel('d',1:168)};
	MarginalCostCali = rgdx('MarginalElectricityPrice',dummygen);
	clear dummygen;
	MargCostCali      = MarginalCostCali.val;
	MargCostCali(:,3) = MargCostCali(:,3)-5;
	NewMargCostCali   = ReshapeMatrix(MargCostCali);
	
	
	% 		marginal cost in RO of bus i for hour d  - added 5   
	dummygen.name  = 'MarginalCostRO';
	dummygen.form  = 'sparse';
	dummygen.uels  = {guel('i',1:312),guel('d',1:168)};
	MarginalCostRO = rgdx('MarginalElectricityPrice',dummygen);
	clear dummygen;
	MargCostRO      = MarginalCostRO.val;
	MargCostRO(:,3) = MargCostRO(:,3)-5;
	NewMargCostRO   = ReshapeMatrix(MargCostRO);
	
	% 		 electricity demand in CA of bus i for hour d  - added 5
	dummygen.name        = 'MultipliedDemandCali';
	dummygen.form        = 'sparse';
	dummygen.uels        = {guel('i',1:312),guel('d',1:168)};
	MultipliedDemandCali = rgdx('MarginalElectricityPrice',dummygen);
	clear dummygen;
	MDemandCali      = MultipliedDemandCali.val;
	MDemandCali(:,3) = MDemandCali(:,3)-5;
	NewMDemandCali   = ReshapeMatrix(MDemandCali);
	
	% 		 electricity demand in RO of bus i for hour d  - added 5
	dummygen.name      = 'MultipliedDemandRO';
	dummygen.form      = 'sparse';
	dummygen.uels      = {guel('i',1:312),guel('d',1:168)};
	MultipliedDemandRO = rgdx('MarginalElectricityPrice',dummygen);
	clear dummygen;
	MDemandRO      = MultipliedDemandRO.val;
	MDemandRO(:,3) = MDemandRO(:,3)-5;
	NewMDemandRO   = ReshapeMatrix(MDemandRO);

	% 		non-served electricity at bus i for hour d
	dummygen.name       = 'UnservedEnergyValue';
	dummygen.form       = 'sparse';
	dummygen.uels       = {guel('i',1:312),guel('d',1:168)};
	UnservedEnergywhole = rgdx('MarginalElectricityPrice',dummygen);
	clear dummygen;
	EnergywholeUnserved      = UnservedEnergywhole.val;
	EnergywholeUnserved(:,3) = EnergywholeUnserved(:,3)-5;
	NewEnergywholeUnserved   = ReshapeMatrix(EnergywholeUnserved);

	% 		electricity generation in CA for hour d 
	dummygen.name  = 'GenerationCali';
	dummygen.form  = 'sparse';
	dummygen.uels  = {guel('d',1:168)};
	GenerationCali = rgdx('MarginalElectricityPrice',dummygen);
	clear dummygen;
	GenerationCaliperHour = (GenerationCali.val(:,2))';
	
	% 		electricity generation in RO for hour d
	dummygen.name = 'GenerationRO';
	dummygen.form = 'sparse';
	dummygen.uels = {guel('d',1:168)};
	GenerationRO  = rgdx('MarginalElectricityPrice',dummygen);
	clear dummygen;
	GenerationROperHour=(GenerationRO.val(:,2))';	

	% 		electricity total variable cost of coal in CA for generator g
	dummygen.name    = 'HeatRateCaliCoal';
	dummygen.form    = 'sparse';
	dummygen.uels    = {guel('g',1:num_gen)};
	HeatRateCaliCoal = rgdx('MarginalElectricityPrice',dummygen);
	clear dummygen;
	CaliforniaHeatRateCoal=(HeatRateCaliCoal.val);	

	% 		electricity total variable cost of natural gas (CT and CCGT) in CA for generator g
	dummygen.name   = 'HeatRateCaliGas';
	dummygen.form   = 'sparse';
	dummygen.uels   = {guel('g',1:num_gen)};
	HeatRateCaliGas = rgdx('MarginalElectricityPrice',dummygen);
	clear dummygen;
	CaliforniaHeatRateGas=(HeatRateCaliGas.val);	
	
	% 		electricity total variable cost of coal in CA for generator g
	dummygen.name  = 'HeatRateROCoal';
	dummygen.form  = 'sparse';
	dummygen.uels  = {guel('g',1:num_gen)};
	HeatRateROCoal = rgdx('MarginalElectricityPrice',dummygen);
	clear dummygen;
	ROHeatRateCoal=(HeatRateROCoal.val);

	% 		electricity total variable cost of natural gas (CT and CCGT) in ROWECC for generator g
	dummygen.name = 'HeatRateROGas';
	dummygen.form = 'sparse';
	dummygen.uels = {guel('g',1:num_gen)};
	HeatRateROGas = rgdx('MarginalElectricityPrice',dummygen);
	clear dummygen;
	ROHeatRateGas=(HeatRateROGas.val);
 
	% --------------------------------------------------------------------
    % calcuate fuel usage for non-renewable generators 
    % --------------------------------------------------------------------
	
	% 		coal fuel usage CA - total - annual
	indices           = ismember(CaliforniaHeatRateCoal(:,1),CoalGenCaliPerGen(:,1)); % determin which generator g has heat rate
	ReqHeatrate       = CaliforniaHeatRateCoal(indices,2);                            % create new varaible with only generators with heat rate
	FuelUsageCoalCali = sum(ReqHeatrate.*CoalGenCaliPerGen(:,2));                     % calculate fuel usage for the year
	z(32) = FuelUsageCoalCali;
	
	% 		natural gas (CT and CCGT) fuel usage CA - total - annual
	indices           = ismember(CaliforniaHeatRateGas(:,1),GasGenCaliPerGen(:,1));
	ReqHeatrate       = CaliforniaHeatRateGas(indices,2);
	FuelUsageGasCali  = sum(ReqHeatrate.*GasGenCaliPerGen(:,2));
	z(33) = FuelUsageGasCali;
	
	% 		coal fuel usage ROWEE - total - annual
	indices           = ismember(ROHeatRateCoal(:,1),CoalGenROPerGen(:,1));
	ReqHeatrate       = ROHeatRateCoal(indices,2);
	FuelUsageCoalRO   = sum(ReqHeatrate.*CoalGenROPerGen(:,2));
	z(34) = FuelUsageCoalRO;
	
	% 		natural gas (CT and CCGT) fuel usage ROWECC - total - annual
	indices           = ismember(ROHeatRateGas(:,1),GasGenROPerGen(:,1));
	ReqHeatrate       = ROHeatRateGas(indices,2);
	FuelUsageGasRO    = sum(ReqHeatrate.*GasGenROPerGen(:,2));
	z(35) = FuelUsageGasRO;

	% --------------------------------------------------------------------
    % slow program for process purposes -- WHY?
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
	
	% ---------------------- REMOVE?
	% 		csv file for  parameters.
	%           filename = ['CaliforniaCoalGen', int2str(week), '.csv'];
	%           csvwrite(filename,CoalCAGen);
	% ----------------------

	% 		csv file for marginal cost in CA of bus i for hour d
	filename = ['MarginalCostCali', int2str(week), '.csv'];
	csvwrite(filename,NewMargCostCali);

	% 		csv file for marginal cost in ROWECC of bus i for hour d
	filename = ['MarginalCostRO', int2str(week), '.csv'];
	csvwrite(filename,NewMargCostRO);
	
	% 		csv file for electricity demand in CA of bus i for hour d
	filename = ['MultipliedDemandCali', int2str(week), '.csv'];
	csvwrite(filename,NewMDemandCali);

	% 		csv file for electricity demand in ROWECC of bus i for hour d
	filename = ['MultipliedDemandRO', int2str(week), '.csv'];
	csvwrite(filename,NewMDemandRO);

	% 		csv file for non-served electricity at bus i for hour d
	filename = ['UnservedEnergywhole', int2str(week), '.csv'];
	csvwrite(filename,NewEnergywholeUnserved);

	% 		csv file for electricity generation in CA for hour d
	filename = ['CaliGenperhour', int2str(week), '.csv'];
	csvwrite(filename,GenerationCaliperHour);

	% 		csv file for electricity generation in CA for hour d
	filename = ['ROGenperhour', int2str(week), '.csv'];
	csvwrite(filename,GenerationROperHour);

	% ---------------------- REMOVE?
	%           %output a csv filw which has the cali heat rate per generation of
	%           %coal
	%           filename = ['CaliHeatrateCoal', int2str(week), '.csv'];
	%           csvwrite(filename,CaliforniaHeatRateCoal);
	%           
	%           %output a csv filw which has the cali heat rate per generation of
	%           %coal
	%           filename = ['CaliHeatrateGas', int2str(week), '.csv'];
	%           csvwrite(filename,CaliforniaHeatRateGas);
	%           
	%           %output a csv filw which has the cali heat rate per generation of
	%           %coal
	%           filename = ['ROHeatrateCoal', int2str(week), '.csv'];
	%           csvwrite(filename,ROHeatRateCoal);
	%           
	%           %output a csv filw which has the cali heat rate per generation of
	%           %coal
	%           filename = ['ROHeatrateGas', int2str(week), '.csv'];
	%           csvwrite(filename,ROHeatRateGas);
	% ----------------------           

	% 		csv file for electricity demand in CA - total for hour d
	filename = ['NewDemandwholeCali', int2str(week), '.csv'];
	csvwrite(filename,NewDemandwholeCali);

	% 		csv file for electricity demand in ROWECC - total for hour d
	filename = ['NewDemandwholeRO', int2str(week), '.csv'];
	csvwrite(filename,NewDemandwholeRO);

	% 		csv file for non-served electricity in CA for hour d
	filename = ['NewEnergywholeunservedCali', int2str(week), '.csv'];
	csvwrite(filename,NewEnergywholeunservedCali);

	% 		csv file for non-served electricity in ROWECC for hour d
	filename = ['NewEnergywholeunservedRO', int2str(week), '.csv'];
	csvwrite(filename,NewEnergywholeunservedRO);

	% 		csv file for shed energy at bus i for hour d
	filename = ['ShedEnergywhole', int2str(week), '.csv'];
	csvwrite(filename,NewEnergywholeShed);

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
