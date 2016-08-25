function [ winds ] = calc_all_city_winds( cities, start_date, end_date )
%CALC_ALL_CITY_WINDS Computes daily wind speed and direction for cities requested
%   CALC_ALL_CITY_WINDS( CITIES, START_DATE, END_DATE ) will compute the
%   wind speed and direction for all cities specified in the structure
%   CITIES, which must contain fields Name, Longitude, and Latitude, i.e.
%   this is the structure returned by CITY_INFO(). It returns a structure
%   with those same fields but with, additionally, the fields dnums, utchr
%   windvel, and winddir. These contain the vectors with the date number,
%   UTC hour that the winds are taken from, wind speed (in m/s) and wind
%   direction in degrees (CCW from E is positive).  Wind vectors are read
%   from WRF-Chem output files and are averaged over a 3x3 grid cell unit
%   centered on the city. In the vertical dimension, it averages over 500 m
%   above ground level.


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

global onCluster
global numThreads
global wrf_path
if isempty(onCluster)
    onCluster = false;
end
if isempty(numThreads) || ~onCluster
    numThreads = 0;
end
if isempty(wrf_path)
    wrf_path = fullfile('/Volumes','share2','USERS','LaughnerJ','WRF','US_BEHR_FULLEMIS','hourly');
end

if ~exist(wrf_path,'dir')
    E.badinput('WRF_PATH %s is not a valid directory', wrf_path)
end

req_fields = {'Longitude', 'Latitude'};
if ~isstruct(cities) || any(~isfield(cities, req_fields))
    E.badinput('CITIES must be a structure with fields %s', strjoin(req_fields, ', '));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DATA LOADING AND PREP %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F = dir(fullfile(wrf_path,'*.nc'));
dnums = nan(numel(F),1);
for b=1:numel(F)
    [s,e] = regexp(F(b).name, '\d\d\d\d-\d\d-\d\d');
    dnums(b) = datenum(F(b).name(s:e),'yyyy-mm-dd');
end
% Remove any days before June 1st - allow WRF spinup
F = F(dnums >= datenum(start_date) & dnums <= datenum(end_date));
dnums = dnums(dnums >= datenum(start_date) & dnums <= datenum(end_date));
[XLON, XLAT, COSALPHA, SINALPHA, zlev, utchr] = read_wrf_vars(wrf_path, F(1), {'XLONG', 'XLAT', 'COSALPHA', 'SINALPHA', 'zlev','utchr'},false,0);

utchr = double(utchr); % imported as int
n_hrs = numel(unique(utchr(:,1)));
unique_hrs = unique(utchr(:,1));

winds = cities;
n_cities = numel(cities);
windvec = nan(numel(dnums), n_hrs);

for b=1:numel(winds)
    winds(b).dnums = dnums;
    winds(b).utchr = unique_hrs'; % transpose so that the dimensions of dnums and utchrs match the corresponding dims of windvel and winddir
    winds(b).windvel = windvec;
    winds(b).winddir = windvec;
end

% We can take just one 2D slice of lon, lat, cos, and sin because these do
% not change in time. U and V we will average for each hour; which hour to
% use can be decided later based on the OMI overpass times.
for b=1:size(utchr,1)
    if any(utchr(b,:) ~= utchr(b,1))
        E.badvar('utchr','Not all files have the same set of utc hours')
    end
end
XLON = XLON(:,:,1,1);
XLAT = XLAT(:,:,1,1);
COSALPHA = COSALPHA(:,:,1,1);
SINALPHA = SINALPHA(:,:,1,1);


% Now calculate the height above ground level for the top of each box, and
% use only Ue and Ve for boxes below 500 m.
zlev = cumsum(zlev,3);
too_high = zlev > 500;


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

local_wrf_path = wrf_path;
city_lons = [winds.Longitude];
city_lats = [winds.Latitude];
windvels = nan(numel(dnums), n_hrs, n_cities);
winddirs = nan(numel(dnums), n_hrs, n_cities);

if numel(F) ~= numel(dnums)
    E.sizeMismatch('F','dnums')
end

% Do not run in parallel if not on the cluster, if numThreads is 0, it will
% run in serial.
%parfor(a=1:numel(F), numThreads)
for a=1:numel(F)
    
    [Ue, Ve] = read_wrf_vars(local_wrf_path, F(a), {'U','V'});
    [Ue, Ve] = wrf_winds_transform(Ue, Ve, COSALPHA, SINALPHA);
    
    % Remove any levels entirely too high, this will cut down on memory
    % usage somewhat, though peak memory usage may not decrease that much.
    Ue(too_high) = nan;
    Ve(too_high) = nan;
    keep_level = true(size(too_high,3),1);
    for i=1:size(too_high,3)
        this_slice = squeeze(too_high(:,:,i,:,:));
        if all(this_slice(:))
            keep_level(i) = false;
        end
    end
    Ue = Ue(:,:,keep_level,:);
    Ve = Ve(:,:,keep_level,:);
    
    for b=1:n_cities
        city_lon = city_lons(b); %#ok<PFBNS>
        city_lat = city_lats(b); %#ok<PFBNS>
        
        
        for h=1:n_hrs
            [this_windvel, this_winddir] = calc_avg_wind(XLON, XLAT, Ue(:,:,:,h), Ve(:,:,:,h), city_lon, city_lat, 'timedim', 4);
            windvels(a,h,b) = this_windvel;
            winddirs(a,h,b) = this_winddir;
        end
    end
end

for b=1:n_cities
    winds(b).windvel = windvels(:,:,b);
    winds(b).winddir = winddirs(:,:,b);
end
end

