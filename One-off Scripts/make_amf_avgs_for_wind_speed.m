% Make gridded AMF averages with ever-increasing cloud fractions
% allowed

work_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed';
cd(fullfile(work_path,'SE US BEHR Hybrid - No ghost'));
str = sprintf('Now in %s',pwd);
accent = repmat('*',1,length(str));
fprintf('\n\t%s\n\t%s\n\t%s\n',accent,str,accent);
[~,NO2_GRID,LON_GRID,LAT_GRID,COUNT] = no2_column_map_2014('2013-06-01','2013-08-30',...
    [-87.1 -81.9], [31.9, 35.6], 'mapfield','BEHRAMFTrop', 'behrdir','.','cloudfraccrit',0.4,'makefig',false);
save('Gridded_avg_AMF_20130601-20130830_Cld04_NormalRow.mat','NO2_GRID','LON_GRID','LAT_GRID','COUNT');
fprintf('Saved hybrid 0.4 cloud\n')
clear NO2_GRID LON_GRID LAT_GRID COUNT

[~,NO2_GRID,LON_GRID,LAT_GRID,COUNT] = no2_column_map_2014('2013-06-01','2013-08-30',...
    [-87.1 -81.9], [31.9, 35.6], 'mapfield','BEHRAMFTrop', 'behrdir','.','cloudfraccrit',0.6,'makefig',false);
save('Gridded_avg_AMF_20130601-20130830_Cld06_NormalRow.mat','NO2_GRID','LON_GRID','LAT_GRID','COUNT');
fprintf('Saved hybrid 0.6 cloud\n')
clear NO2_GRID LON_GRID LAT_GRID COUNT

[~,NO2_GRID,LON_GRID,LAT_GRID,COUNT] = no2_column_map_2014('2013-06-01','2013-08-30',...
    [-87.1 -81.9], [31.9, 35.6], 'mapfield','BEHRAMFTrop', 'behrdir','.','cloudfraccrit',0.8,'makefig',false);
save('Gridded_avg_AMF_20130601-20130830_Cld08_NormalRow.mat','NO2_GRID','LON_GRID','LAT_GRID','COUNT');
fprintf('Saved hybrid 0.8 cloud\n')
clear NO2_GRID LON_GRID LAT_GRID COUNT

cd(fullfile(work_path,'SE US BEHR Monthly - No ghost'));
str = sprintf('Now in %s',pwd);
accent = repmat('*',1,length(str));
fprintf('\n\t%s\n\t%s\n\t%s\n',accent,str,accent);
[~,NO2_GRID,LON_GRID,LAT_GRID,COUNT] = no2_column_map_2014('2013-06-01','2013-08-30',...
    [-87.1 -81.9], [31.9, 35.6], 'mapfield','BEHRAMFTrop', 'behrdir','.','cloudfraccrit',0.4,'makefig',false);
save('Gridded_avg_AMF_20130601-20130830_Cld04_NormalRow.mat','NO2_GRID','LON_GRID','LAT_GRID','COUNT');
fprintf('Saved monthly 0.4 cloud\n')
clear NO2_GRID LON_GRID LAT_GRID COUNT

[~,NO2_GRID,LON_GRID,LAT_GRID,COUNT] = no2_column_map_2014('2013-06-01','2013-08-30',...
    [-87.1 -81.9], [31.9, 35.6], 'mapfield','BEHRAMFTrop', 'behrdir','.','cloudfraccrit',0.6,'makefig',false);
save('Gridded_avg_AMF_20130601-20130830_Cld06_NormalRow.mat','NO2_GRID','LON_GRID','LAT_GRID','COUNT');
fprintf('Saved monthly 0.6 cloud\n')
clear NO2_GRID LON_GRID LAT_GRID COUNT

[~,NO2_GRID,LON_GRID,LAT_GRID,COUNT] = no2_column_map_2014('2013-06-01','2013-08-30',...
    [-87.1 -81.9], [31.9, 35.6], 'mapfield','BEHRAMFTrop', 'behrdir','.','cloudfraccrit',0.8,'makefig',false);
save('Gridded_avg_AMF_20130601-20130830_Cld08_NormalRow.mat','NO2_GRID','LON_GRID','LAT_GRID','COUNT');
fprintf('Saved monthly 0.8 cloud\n')
clear NO2_GRID LON_GRID LAT_GRID COUNT

cd(fullfile(work_path,'SE US BEHR Hourly - No ghost'));
str = sprintf('Now in %s',pwd);
accent = repmat('*',1,length(str));
fprintf('\n\t%s\n\t%s\n\t%s\n',accent,str,accent);
[~,NO2_GRID,LON_GRID,LAT_GRID,COUNT] = no2_column_map_2014('2013-06-01','2013-08-30',...
    [-87.1 -81.9], [31.9, 35.6], 'mapfield','BEHRAMFTrop', 'behrdir','.','cloudfraccrit',0.4,'makefig',false);
save('Gridded_avg_AMF_20130601-20130830_Cld04_NormalRow.mat','NO2_GRID','LON_GRID','LAT_GRID','COUNT');
fprintf('Saved hourly 0.4 cloud\n')
clear NO2_GRID LON_GRID LAT_GRID COUNT

[~,NO2_GRID,LON_GRID,LAT_GRID,COUNT] = no2_column_map_2014('2013-06-01','2013-08-30',...
    [-87.1 -81.9], [31.9, 35.6], 'mapfield','BEHRAMFTrop', 'behrdir','.','cloudfraccrit',0.6,'makefig',false);
save('Gridded_avg_AMF_20130601-20130830_Cld06_NormalRow.mat','NO2_GRID','LON_GRID','LAT_GRID','COUNT');
fprintf('Saved hourly 0.4 cloud\n')
clear NO2_GRID LON_GRID LAT_GRID COUNT

[~,NO2_GRID,LON_GRID,LAT_GRID,COUNT] = no2_column_map_2014('2013-06-01','2013-08-30',...
    [-87.1 -81.9], [31.9, 35.6], 'mapfield','BEHRAMFTrop', 'behrdir','.','cloudfraccrit',0.8,'makefig',false);
save('Gridded_avg_AMF_20130601-20130830_Cld08_NormalRow.mat','NO2_GRID','LON_GRID','LAT_GRID','COUNT');
fprintf('Saved hourly 0.4 cloud\n')
clear NO2_GRID LON_GRID LAT_GRID COUNT