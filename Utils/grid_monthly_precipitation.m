% Grid_monthly_temperature: Interpolates NOAA temperature data to a lat/lon
% grid for plotting temperature to compare against BEHR maps

start_date = '04/01/2011';
end_date = '09/30/2011';

lat_bdy = [25 50];
lon_bdy = [-125 -65];
resolution = 0.05;

date_mode = 'precip'; % set to monthly or daily

temp_data_dir = '/Users/Josh/Documents/MATLAB/Data for BEHR/Temperature Data';


warning off all

data_ind = 1; % Use this to keep track of how many temperature data entries we've read in

for d = 1:numel(datenums)
    curr_day = datestr(datenums(d),29);
    year = curr_day(1:4);
    month = curr_day(6:7);
    day = curr_day(9:10);
    
    % Open and read in the station list
    full_station_dir = fullfile(temp_data_dir,year,month,sprintf('%s%sstation.txt',year,month));
    fid_station = fopen(full_station_dir);
    fgetl(fid_station); % Skip the header row
    stations = textscan(fid_station,'%d %d %s %d %d %d %s %s %s %f %f %d %d %d %d','Delimiter','|');
    fclose(fid_station);
    
    % Save the relevant fields as more useful variables
    wban_id = stations{1};
    station_lats = stations{10};
    station_lons = stations{11};
    
    clear stations
    
    % The first time through, use the number of stations and datenums vector to
    % guess how many sites we will have. If extra need to be added, they will
    % just get tacked on at the end, and any unused will be removed later in
    % the code.
    if d==1;
        temp_vec = NaN(1,numel(wban_id)*numel(datenums)*24);
        lat_vec = NaN(1,numel(wban_id)*numel(datenums)*24);
        lon_vec = NaN(1,numel(wban_id)*numel(datenums)*24);
    end
    
    % Open the desired dataset (monthly or daily)
    full_data_dir = fullfile(temp_data_dir,year,month,sprintf('%s%s%s.txt',year,month,date_mode));
    fid_data = fopen(full_data_dir);
    fgetl(fid_data); % Grab the head line
    
    data = textscan(fid_data,'%d %d %d %d');
    x=1;
    
end

% Remove any unfilled temperature values
nans = find(isnan(temp_vec));
temp_vec(nans) = [];
lat_vec(nans) = [];
lon_vec(nans) = [];

% Make our lat/lon grid and interpolate the data to it
lons = min(lon_bdy):resolution:max(lon_bdy);
lats = min(lat_bdy):resolution:max(lat_bdy);
[longrid, latgrid] = meshgrid(lons,lats);
gridded_temp = griddata(lon_vec, lat_vec, temp_vec, longrid, latgrid);
