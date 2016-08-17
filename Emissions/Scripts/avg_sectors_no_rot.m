function [ Sectors ] = avg_sectors_no_rot( start_date, end_date, city, varargin )
%AVG_SECTORS_NO_ROT Averages days when the wind is in one of 8 sectors

DEBUG_LEVEL = 1;
HOME = getenv('HOME');

p = inputParser;
p.addParameter('behrdir',fullfile(HOME,'Documents', 'MATLAB', 'BEHR', 'Workspaces','Wind speed','SE US BEHR Hourly - No ghost'),@(x) exist(x,'dir'))
p.addParameter('fileprefix','OMI_BEHR_',@isstr);
p.addParameter('clouds','omi',@isstr);
p.addParameter('cloudfraccrit',-1,@isscalar)
p.addParameter('rowanomaly','XTrackFlags',@(x) any(strcmpi(x,{'AlwaysByRow','RowsByTime','XTrackFlags','XTrackFlagsLight'}))) %Ensure that the rowanomaly value is one of the allowed 4
p.addParameter('rows',[],@(x) (isnumeric(x) && (numel(x) == 0 || numel(x) == 2)));

p.parse(varargin{:});
pout = p.Results;

behrdir = pout.behrdir;
behrprefix = pout.fileprefix;
cloud_type = pout.clouds;
if ~ismember(lower(cloud_type), {'omi', 'modis'})
    E.badinput('%s is not a valid cloud type', cloud_type);
end
cloud_frac_crit = pout.cloudfraccrit;
if cloud_frac_crit < 0 
    switch lower(cloud_type)
        case 'omi'
            cloud_frac_crit = 0.2;
        case 'modis'
            cloud_frac_crit = 0;
    end
end
row_anomaly_type = pout.rowanomaly;
rows_reject = pout.rows;



% Load the proper wind file for the given city
[theta, theta_dnums] = load_windfile(city);

% Set the bounds for the city to be within ~200 km
[lonlim, latlim] = city_lims(city);

dnums = datenum(start_date):datenum(end_date);

init_bool = true;

for d=1:numel(dnums);
    if DEBUG_LEVEL > 0; fprintf('Adding VCDs from %s to ', datestr(dnums(d))); end
    % Load the file for the given day, figure out its wind direction
    fname = sprintf('%s%04d%02d%02d.mat',behrprefix,year(dnums(d)),month(dnums(d)),day(dnums(d)));
    F = load(fullfile(behrdir, fname), 'OMI');
    OMI = F.OMI;
    if init_bool
        init_bool = false;
        xx = OMI(1).Longitude(1,:) >= lonlim(1) & OMI(1).Longitude(1,:) <= lonlim(2);
        yy = OMI(1).Latitude(:,1) >= latlim(1) & OMI(1).Latitude(:,1) <= latlim(2);
        
        blank_mat = nan(sum(yy), sum(xx));
        lon_mat = OMI(1).Longitude(yy,xx);
        lat_mat = OMI(1).Latitude(yy,xx);
        SectSubstruct = struct('VCDs', blank_mat, 'sum_weighted_col', blank_mat,...
            'sum_weight', blank_mat, 'Longitude', lon_mat, 'Latitude', lat_mat);
        Sectors = struct('E', SectSubstruct, 'NE', SectSubstruct, 'N', SectSubstruct, 'NW',SectSubstruct,...
            'W',SectSubstruct,'SW',SectSubstruct,'S',SectSubstruct,'SE',SectSubstruct);
    else
        lontest = OMI(1).Longitude(yy,xx);
        lattest = OMI(1).Latitude(yy,xx);
        if any(lontest(:) ~= lon_mat(:)) || any(lattest(:) ~= lat_mat(:))
            E.callError('inconsistent_lon-lat','OMI longitude/latitude grid for %s is inconsistent with previous days',datestr(dnums(d)));
        end
    end
    
    dd = theta_dnums == dnums(d);
    if sum(dd) < 1
        E.callError('bad_wind','Wind direction for %s not found',datestr(dnums(d)));
    elseif sum(dd) > 1
        E.callError('bad_wind','Multiple entries for wind direction on %s found', datestr(dnums(d)));
    end
    todays_sect = theta_to_sect(theta(dd));
    if DEBUG_LEVEL > 0; fprintf('%s\n',todays_sect); end
    [this_WeightedColumn, this_Weight] = BEHR_day_no2(OMI,'mapfield', 'BEHRColumnAmountNO2Trop',...
        'cloud_prod', cloud_type, 'cloud_frac_max', cloud_frac_crit,...
        'row_anomaly', row_anomaly_type, 'rows', rows_reject);
    
    Sectors.(todays_sect).sum_weighted_col = nansum2(cat(3,Sectors.(todays_sect).sum_weighted_col, this_WeightedColumn(yy,xx)),3);
    Sectors.(todays_sect).sum_weight = nansum2(cat(3,Sectors.(todays_sect).sum_weight, this_Weight(yy,xx)),3);
end

% Go through each sector and get the final VCDs
fns = fieldnames(Sectors);
for f=1:numel(fns)
    Sectors.(fns{f}).VCDs = Sectors.(fns{f}).sum_weighted_col ./ Sectors.(fns{f}).sum_weight;
end

end

function [theta, dnums] = load_windfile(city)
HOME = getenv('HOME');
root_path = fullfile(HOME, 'Documents', 'MATLAB', 'BEHR', 'Workspaces', 'Wind speed');
% Ensure city starts with a capital letter
city = strcat(upper(city(1)), lower(city(2:end)));
fname = sprintf('%s-Wind-Conditions-1900UTC-5layers-earthrel.mat', city);
W = load(fullfile(root_path, fname), 'theta', 'dnums');
theta = W.theta;
dnums = W.dnums;
end

function [lonlim, latlim] = city_lims(city)
E = JLLErrors;
switch lower(city)
    case 'atlanta'
        lonlim = [-86.4 -82.4];
        latlim = [31.80 35.80];
    case 'birmingham'
        lonlim = [-88.80 -84.80];
        latlim = [31.5 35.5];
    otherwise
        E.badinput('%s not a recognized city', city)
end
end

function sect = theta_to_sect(theta)
E = JLLErrors;
if theta < -180 || theta > 180
    E.badinput('theta should be between -180 and +180 (theta = %g)',theta)
elseif theta >= -22.5 && theta < 22.5
    sect = 'E';
elseif theta >= 22.5 && theta < 67.5
    sect = 'NE';
elseif theta >= 67.5 && theta < 112.5
    sect = 'N';
elseif theta >= 112.5 && theta < 157.5
    sect = 'NW';
elseif theta >= -67.5 && theta < -22.5
    sect = 'SE';
elseif theta >= -112.5 && theta < -67.5
    sect = 'S';
elseif theta >= -157.5 && theta < -112.5
    sect = 'SW';
elseif theta < -157.5 || theta >= 157.5
    sect = 'W';
end
end