function [ data ] = read_modis_cloud( modis_directory, date_in, data, omi_swath_start, omi_next_swath_start, lonlim, latlim, varargin )
%READ_MODIS_CLOUD Reads MODIS MYD06_L2 files and averages the cloud fraction to OMI pixels
%   DATA = READ_MODIS_CLOUD( MODIS_DIRECTORY, DATE_IN, DATA,
%   OMI_SWATH_START, OMI_NEXT_SWATH_START, LONLIM, LATLIM) Reads MYD06_L2
%   files' Cloud_Fraction dataset and averages those values to OMI pixels.
%   The result is stored in the MODISCloud field of DATA. The MODIS files
%   used are stored in a cell array in the MODISCloudFiles field of DATA.
%   Description of all inputs follows:
%
%       MODIS_DIRECTORY: the directory containing the MYD06_L2 files, which
%       should be organized by year in subdirectories.
%
%       DATE_IN: a representation of the date to average for as either a
%       Matlab datenum or a string that can be parsed by DATENUM() without
%       a format specified.
%
%       DATA: a structure to which the MODIS data will be added. It must
%       contain fields MODISCloud and MODISCloudFiles to accept the data,
%       plus the field Longitude, Latitude, and the pixel corner fields
%       specified by the parameters 'LoncornField' and 'LatcornField'.
%
%       OMI_SWATH_START: the starting UTC time of the OMI swath that is
%       being operated on (the number following the first "t" in the NASA
%       SP filename) given as a number (not a string). Running str2double()
%       on the time given in the filename will return the appropriate
%       value, which should be HHMM or HMM (depending if the hour has two
%       digits or not). This along with OMI_NEXT_SWATH_START are used to
%       determine which MODIS cloud files are temporally coincident with
%       the current swath.
%
%       OMI_NEXT_SWATH_START: the starting UTC time of the next OMI swath
%       after the one being operated on currently, again given as a number.
%
%       LONLIM: a two element vector giving the longitude limits of the
%       retrieval domain. The more negative element should go first.
%
%       LATLIM: same as LONLIM but for the latitude limits.
%
%   Parameter values:
%
%       'DEBUG_LEVEL': takes a scalar number, sets the verbosity of the
%       function. Defaults to 0 (i.e. no debugging messages printed).
%       Higher numbers increase the verbosity.
%
%       'AllowNoFile': boolean (default is false); if true, will not error
%       if no MODIS cloud files are found.
%
%       'LoncornField': indicates which field describes the corner points
%       of the pixels. This field must be 3D and have size == 
%       [4, size(DATA.Latitude)]. Defaults to 'FoV75CornerLongitude'.
%
%       'LatcornField': similar to 'LoncornField', but for latitude.
%       Defaults to 'FoV75CornerLatitude'.

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
p = inputParser;
p.addParameter('DEBUG_LEVEL', 0);
p.addParameter('AllowNoFile', false);
p.addParameter('LoncornField', 'FoV75CornerLongitude');
p.addParameter('LatcornField', 'FoV75CornerLatitude');
p.parse(varargin{:});
pout = p.Results;

DEBUG_LEVEL = pout.DEBUG_LEVEL;
allow_no_file = pout.AllowNoFile;
loncorn_field = pout.LoncornField;
latcorn_field = pout.LatcornField;

% Validation - positional inputs %

if ~ischar(modis_directory)
    E.badinput('MODIS_DIRECTORY must be a string')
elseif ~exist(modis_directory, 'dir')
    E.badinput('MODIS_DIRECTORY is not a directory')
end


if isnumeric(date_in)
    if ~isscalar(date_in)
        E.badinput('If given as a number, DATE_IN must be scalar')
    end
elseif ischar(date_in)
    try
        date_in = datenum(date_in);
    catch err
        if strcmp(err.identifier, 'MATLAB:datenum:ConvertDateString')
            E.badinput('DATE_IN could not be recognized as a valid format for a date string')
        else
            rethrow(err)
        end
    end
else
    E.badinput('DATE_IN must be a string or number')
end


if ~isnumeric(omi_swath_start) || ~isscalar(omi_swath_start) || omi_swath_start < 0 || omi_swath_start > 2359
    E.badinput('OMI_SWATH_START must be a scalar number between 0 and 2359')
end
if ~isnumeric(omi_next_swath_start) || ~isscalar(omi_next_swath_start) || omi_next_swath_start < 0 || omi_next_swath_start > 2359
    E.badinput('OMI_NEXT_SWATH_START must be a scalar number between 0 and 2359')
end

if ~isnumeric(lonlim) || numel(lonlim) ~= 2 || lonlim(1) > lonlim(2)
    E.badinput('LONLIM must be a two element numeric vector with the second element greater than the first')
end
if ~isnumeric(latlim) || numel(latlim) ~= 2 || latlim(1) > latlim(2)
    E.badinput('LATLIM must be a two element numeric vector with the second element greater than the first')
end

% Validation - Parameters %

if ~isnumeric(DEBUG_LEVEL) || ~isscalar(DEBUG_LEVEL)
    E.badinput('The parameter DEBUG_LEVEL must be a scalar number')
end

if ~islogical(allow_no_file) || ~isscalar(allow_no_file)
    E.badinput('The parameter AllowNoFile must be a scalar logical value');
end

if ~ischar(loncorn_field) 
    E.badinput('The parameter LonconrField must be a string')
elseif ~isfield(data, loncorn_field)
    E.badinput('The LoncornField "%s" is not in DATA', loncorn_field)
end

if ~ischar(latcorn_field) 
    E.badinput('The parameter LatcornField must be a string')
elseif ~isfield(data, latcorn_field)
    E.badinput('The LatcornField "%s" is not in DATA', latcorn_field)
end

% If you ever read more data fields from the Data structure, or save more
% fields to the Data structre, add them to this list.
req_fields = {'Longitude', 'Latitude', 'MODISCloud', 'MODISCloudFiles', loncorn_field, latcorn_field};
if ~isstruct(data)
    E.badinput('DATA must be a structure')
end
xx = ~ismember(req_fields, fieldnames(data));
if any(xx)
    E.badinput('DATA is missing the required field(s) %s. The lat/lon corner fields (%s, %s) can be changed by the corresponding parameters (see function help)',...
        strjoin(req_fields(xx), ', '), loncorn_field, latcorn_field)
end

% The latcorn and loncorn fields should have the corners as the field
% dimension and, then the same size as the Longitude fields
sz_loncorn = size(data.(loncorn_field));
sz_latcorn = size(data.(latcorn_field));
sz_check = [4, size(data.Latitude)];

if ~isequal(sz_loncorn, sz_check)
    E.badinput('The shape of the longitude corner field (%s) is not %s (i.e. corners then same size as Longitude)', loncorn_field, mat2str(sz_loncorn));
end
if ~isequal(sz_latcorn, sz_check)
    E.badinput('The shape of the latitude corner field (%s) is not %s (i.e. corners then same size as Longitude)', latcorn_field, mat2str(sz_latcorn));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

julian_day = modis_date_to_day(date_in);

%Find all MODIS files that occur after the current OMI file
%but before the next OMI file.
modis_filepattern = sprintf('MYD06_L2.A%04d%03d*.hdf',year(date_in),julian_day);
year_str = sprintf('%04d', year(date_in));
modis_files = dir(fullfile(modis_directory, year_str, modis_filepattern));

% Error if no MODIS cloud files found: this will prevent publishing data
% without MODIS cloud information.
if isempty(modis_files)
    if allow_no_file
        % If no MODIS cloud files found, then fill the MODISCloud field
        % with NaNs (representing no data)
        data.MODISCloud = nan(size(data.Latitude));
        return
    else
        % In most cases though, if there's no files, we shouldn't be
        % retrieving this day yet (not until the cloud files are available)
        % so error.
        E.filenotfound(sprintf('MODIS cloud for %s',datestr(date_in)));
    end
end

% Where in the filename the time for that granule is given. Assumes file
% name like MYD06_L2.A2013042.1805.006.2014112100744.hdf, where "1805" is
% the time in question
mod_time_indices = 19:22;

all_lon = [];
all_lat = [];
all_cldfrac = [];

files_used = false(size(modis_files));

for a=1:length(modis_files);
    mod_filename = modis_files(a).name;
    mod_granule_time = str2double(mod_filename(mod_time_indices));
    % Skip any modis files that do not occur during the
    % time period of the current swath
    if mod_granule_time < omi_swath_start;
        continue
    elseif mod_granule_time > omi_next_swath_start;
        continue
    end
    
    %For each file that fits the criteria mentioned
    %above, import its latitude, longitude, and cloud
    %fraction.
    if DEBUG_LEVEL > 0; fprintf('  Averaging MODIS cloud file %s\n',mod_filename); end
    mod_filename = fullfile(modis_directory, year_str, modis_files(a).name);
    mod_fileinfo = hdfinfo(mod_filename);
     
    % At the moment, since we're only reading three variables, there's no
    % reason to do something like read_omi_sp that can read an arbitrary
    % list of variables. This is also a different application because we
    % aren't just copying variables from one to another.
    latitude = hdfreadmodis(mod_filename, hdfdsetname(mod_fileinfo, 1, 1, 'Latitude'));
    longitude = hdfreadmodis(mod_filename, hdfdsetname(mod_fileinfo, 1, 1, 'Longitude'));
    cloud_fraction = hdfreadmodis(mod_filename, hdfdsetname(mod_fileinfo, 1, 2, 'Cloud_Fraction'));
    
    xx = longitude > lonlim(1) - 0.25 & longitude < lonlim(2) + 0.25;
    yy = latitude > latlim(1) - 0.25 & latitude < latlim(2) + 0.25;
    
    longitude = longitude(xx & yy);
    latitude = latitude(xx & yy);
    cloud_fraction = cloud_fraction(xx & yy);
    
    if ~isempty(cloud_fraction)
        all_lon = cat(1, all_lon, longitude);
        all_lat = cat(1, all_lat, latitude);
        all_cldfrac = cat(1, all_cldfrac, cloud_fraction);
        files_used(a) = true;
    end
end

%Find all the MODIS cloud pixels in each OMI pixel and average them
%together. Pixels that have no MODIS clouds available (towards the edge of
%the swath) should have a NaN.
data.MODISCloud = nan(size(data.Latitude));
if ~isempty(all_cldfrac)
    loncorn = data.(loncorn_field);
    latcorn = data.(latcorn_field);
    
    for a = 1:numel(data.Latitude);
        xall=[loncorn(: ,a); loncorn(1, a)];
        yall=[latcorn(:, a); latcorn(1, a)];
        
        % If there is an invalid corner coordinate, skip because we cannot
        % be sure the correct polygon will be used.
        if any(isnan(xall)) || any(isnan(yall))
            continue
        end
        
        xx_cld = inpolygon(all_lat, all_lon, yall, xall);
        
        cld_vals=all_cldfrac(xx_cld);
        data.MODISCloud(a)=nanmean(cld_vals);
    end
    
    data.MODISCloudFiles={modis_files(files_used).name};
    for a = 1:numel(data.MODISCloudFiles)
        data.MODISCloudFiles{a} = fullfile(modis_directory, data.MODISCloudFiles{a});
    end
end

end

