function [ windvel, theta, avail_dnums ] = read_ecmwf( files, center_lon, center_lat )
%READ_ECMWF Read and average wind data from an ECMWF netCDF file
%   WINDS = READ_ECMWF( FILES, CENTER_LON, CENTER_LAT ) Read in the U and V
%   wind data from the netCDF file given by FILES and find the
%   average wind at OMI overpass time for a 3x3 set of wind data grid cells
%   around the CENTER_LON and CENTER_LAT. File must be a single file name
%   as a string

windvel_cell = cell(size(files));
theta_cell = cell(size(files));
avail_dnums_cell = cell(size(files));

for f=1:numel(files)
    check_file(files{f});
    [windvel_cell{f}, theta_cell{f}, avail_dnums_cell{f}] = read_single_ecmwf(files{f}, center_lon, center_lat);
end

windvel = veccat(windvel_cell{:});
theta = veccat(theta_cell{:});
avail_dnums= veccat(avail_dnums_cell{:});

% Ensure outputs are in order by date
[avail_dnums, sortvec] = sort(avail_dnums);
windvel = windvel(sortvec);
theta = theta(sortvec);

end

function check_file(filename)
E = JLLErrors;
if ~exist(filename,'file')
    E.filenotfound(filename);
end

try
    ni = ncinfo(filename);
catch err
    if strcmp(err.identifier,'MATLAB:imagesci:netcdf:libraryFailure')
        E.callError('not_ecmwf','%s does not appear to be a netcdf file',filename);
    else
        rethrow(err)
    end
end

file_vars = {ni.Variables.Name};
req_vars = {'longitude','latitude','level','time','u','v'};
xx = ~ismember(req_vars, file_vars);
if any(xx)
    E.callError('not_ecmwf','%s does not contain the expected variables: %s', filename, strjoin(req_vars(xx), ', '));
end
if ~strcmpi(ncreadatt(ni.Filename, 'level', 'long_name'), 'pressure_level')
    E.callError('ecmwf_level', '%s does not have pressure as the vertical coordinate', filename)
end
end

function [ windvel, theta, avail_dnums ] = read_single_ecmwf( file, center_lon, center_lat )


%%%%%%%%%%%%%%%%%%%%%
%%%%% CONSTANTS %%%%%
%%%%%%%%%%%%%%%%%%%%%

omi_overpass = 13.75;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
if ~ischar(file)
    E.badinput('FILES must be a string')
elseif ~exist(file, 'file')
    E.filenotfound(file)
end

if ~isnumeric(center_lon) || ~isscalar(center_lon) || center_lon < -180 || center_lon > 180
    E.badinput('CENTER_LON must be a scalar number between -180 and 180.')
end

if ~isnumeric(center_lat) || ~isscalar(center_lat) || center_lat < -90 || center_lat > 90
    E.badinput('CENTER_LAT must be a scalar number between -90 and 90.')
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Pressure that corresponds to 500 m altitude. Use to subset vertical
% levels
p500m = 1013 * exp(-0.5/7.4);

% # of grid cells to either side of the center one
box_sz = 1;

ni = ncinfo(file);
plev = double(ncread(ni.Filename, 'level'));
zz = plev >= p500m;
zi = find(zz, 1);
ze = find(zz, 1, 'last');

nclon = double(ncread(ni.Filename, 'longitude'));
nclon(nclon>180) = nclon(nclon>180)-360;
nclat = double(ncread(ni.Filename, 'latitude'));
% ECMWF time is given in hours since 1 Jan 1900. This adjusts it to be
% consistent with Matlab's datenums
nctime = double(ncread(ni.Filename, 'time'));
tt = nctime > 0;
if any(~tt)
    warning('Times are 0 in %s, these will be ignored', file)
end
ti = find(tt,1);
te = find(tt,1,'last');

ncdnum = nctime(tt)/24 + datenum('1900-01-01');
avail_dnums = unique(floor(ncdnum));

% Find the point closest in lat/lon
[~,xi] = min(abs(nclon - center_lon));
[~,yi] = min(abs(nclat - center_lat));

start = [xi - box_sz, yi - box_sz, zi, ti];
count = [2*box_sz+1, 2*box_sz+1, ze-zi+1, te-ti+1];

% NCREAD handles offset, scaling, and fill values automatically
ncU = ncread(ni.Filename, 'u', start, count);
ncV = ncread(ni.Filename, 'v', start, count);

% For each time, compute the average wind field from the grid cells around
% the center one.
windvel_all = nan(size(ncdnum));
theta_all = nan(size(ncdnum));
for a=1:numel(ncdnum)
    u_field = nanmean(ncU(:,:,:,a),3);
    u_bar = nanmean(u_field(:));
    v_field = nanmean(ncV(:,:,:,a),3);
    v_bar = nanmean(v_field(:));
    
    windvel_all(a) = sqrt(u_bar.^2 + v_bar.^2);
    theta_all(a) = atan2d(v_bar, u_bar);
end

% Finally compute the approximate OMI overpass time, assuming an overpass
% of 13:45 local standard. Then interpolate the ECMWF wind fields to that
% time.

utc_offset = round(center_lon/15);
utc_overpass = omi_overpass - utc_offset;
avail_overpass = avail_dnums + utc_overpass/24;

windvel = interp1(ncdnum, windvel_all, avail_overpass);
theta = interp1(ncdnum, theta_all, avail_overpass);

end

