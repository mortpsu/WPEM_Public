function output=Datetimestamp(z1,IterationCount)
t1 = datetime(2010,1,1,1,0,0);
t2 = datetime(2010,12,30,24,0,0);
timeValue = t1:hours(1):t2;

for iTime = 1:length(timeValue)
   Time{iTime} = datestr(timeValue(iTime), 'yyyy-mm-dd-HH-MM') ; 
end

%For seasons use inbetween
%The below is for the hours of the day
%Extract rows for morning
   morning=timeValue(hour(timeValue)>=6 & hour(timeValue)<=11);
%Extract rows for afternoon
 afternoon=timeValue(hour(timeValue)>=12 & hour(timeValue)<=17);
 %Extract rows for evening
   evening=timeValue(hour(timeValue)>=18 & hour(timeValue)<=23);
%Extract rows for afternoon
 night=timeValue(hour(timeValue)>=0 & hour(timeValue)<=5) ;
 


Timeofday=cell(8736,1);
for ii=1:1:8736
    count = ismember(timeValue,morning);
    if count(ii)==1
        Timeofday{ii}='morning';
    end
   
    
end

for ii=1:1:8736
    count = ismember(timeValue,afternoon);
    if count(ii)==1
        Timeofday{ii}='afternoon';
    end
   
    
end

for ii=1:1:8736
    count = ismember(timeValue,evening);
    if count(ii)==1
        Timeofday{ii}='evening';
    end
   
    
end

for ii=1:1:8736
    count = ismember(timeValue,night);
    if count(ii)==1
        Timeofday{ii}='night';
    end
   
    
end

%Below is for day of the week
[DayNumber,DayName] = weekday(timeValue);
Temp=cellstr(DayName);
    for ii=1:1:8736
        if Temp{ii}=='Mon'
            NameDay{ii}='weekday';
            
        elseif Temp{ii}=='Tue'
            NameDay{ii}='weekday';
        elseif Temp{ii}=='Wed'
            NameDay{ii}='weekday';
        elseif Temp{ii}=='Thu'
            NameDay{ii}='weekday';
        elseif Temp{ii}=='Fri'
            NameDay{ii}='weekday';    
        elseif Temp{ii}=='Sat'
            NameDay{ii}='weekend';
        elseif Temp{ii}=='Sun'
            NameDay{ii}='weekend';    
        end
    end
    
    %Now for seasons
    
    %First is spring
    tlower = datetime(2010,03,21,1,0,0);
    tupper = datetime(2010,06,21,24,0,0);
    
    spring = isbetween(timeValue,tlower,tupper);
    
for ii=1:1:8736
    if spring(ii)==1
        seasonofyear{ii}='spring';
    end
   
    
end
    
    tlower = datetime(2010,06,22,1,0,0);
    tupper = datetime(2010,09,22,24,0,0);
    
    summer = isbetween(timeValue,tlower,tupper);
    
for ii=1:1:8736
    if summer(ii)==1
        seasonofyear{ii}='summer';
    end
   
    
end 
    tlower = datetime(2010,09,23,1,0,0);
    tupper = datetime(2010,12,20,24,0,0);
    
    fall = isbetween(timeValue,tlower,tupper);

for ii=1:1:8736
    if fall(ii)==1
        seasonofyear{ii}='fall';
    end
   
    
end
    tlower1 = datetime(2010,12,21,1,0,0);
    tupper1 = datetime(2010,12,30,24,0,0);
    
    tlower2 = datetime(2010,01,01,1,0,0);
    tupper2 = datetime(2010,03,20,24,0,0);
    
    winter1 = isbetween(timeValue,tlower1,tupper1); 
    winter2 = isbetween(timeValue,tlower2,tupper2);
    
    winter = winter1 | winter2;

for ii=1:1:8736
    if winter(ii)==1
        seasonofyear{ii}='winter';
    end
   
    
end

%z1=xlsread('Econ_data_Larger_Shock.xlsx');

%       create head for csv output
% identifiy the number of regions
num_econr = size(z1,2)/2;

% create header titles
header1 = {'Time','Season','Dayofweek','TimeofDay'};                    % time headers
cell_econr = cellstr(string(1:1:num_econr));                            % write econr as cell
header2    = matlab.lang.makeValidName(cell_econr,'Prefix','mdemand_'); % write demand header names per region
header3    = matlab.lang.makeValidName(cell_econr,'Prefix','nse_');     % write nse header names per region

EconHeader2 = [header1, header2 header3];
% EconHeader2 = {'Time','Season','Dayofweek','TimeofDay','DemandCali','DemandRO','UnservedEnergyCali','UnservedEnergyRO'}; 
commaHeader2 = [EconHeader2;repmat({','},1,numel(EconHeader2))]; 
commaHeader2 = commaHeader2(:)';
EcontextHeader2 = cell2mat(commaHeader2); 

z1_col = size(z1,2);
fmt_spec = [repmat('%12.3f, ', 1, z1_col-1),'%12.3f\n'];

cd results
filename = ['Iteration', int2str(IterationCount)];
cd(filename)
filename = ['UnservedEnergyData', int2str(IterationCount), '.csv'];
fid = fopen(filename, 'w') ;
for iLine = 1:size(z1, 1)+1 % Loop through each time/value row
    if iLine ==1
     fprintf(fid,'%s\n',EcontextHeader2);
    else
   fprintf(fid, '%s,', Time{iLine-1}) ;            % Print the time string
   fprintf(fid, '%s,', seasonofyear{iLine-1}) ;    % Print the season string
   fprintf(fid, '%s,', NameDay{iLine-1}) ;         % Print the day of week string
   fprintf(fid, '%s,', Timeofday{iLine-1}) ;       % Print the timeofday string 
   fprintf(fid, fmt_spec, z1(iLine-1, 1:z1_col)) ; % Print the data values
%     fprintf(fid, '%12.3f, %12.3f,%12.3f,%12.3f\n', z1(iLine-1, 1:size(z1,2))) ; 
    % Print the data values

    end    
end
fclose(fid) ;
cd ..
cd ..
output=1;
end
