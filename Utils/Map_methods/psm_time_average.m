function [ values, lon_grid, lat_grid ] = psm_time_average( start_date, end_date, varargin  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARSING AND VALIDATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addOptional('lon_bdy', GlobeGrid(0.05, 'domain', 'us'));
p.addOptional('lat_bdy', []);
p.addParameter('behr_dir', behr_paths.behr_mat_dir);
p.addParameter('filepattern', sprintf('OMI_BEHR_%s_*.mat', BEHR_version));
p.addParameter('dayofweek', 'UMTWRFS');
p.addParameter('holidays', false);
p.addParameter('filterpsm', false);
p.addParameter('avgfield', 'BEHRColumnAmountNO2Trop');
p.addParameter('rejectmode', 'detailed');
p.addParameter('clouds', 'omi');
p.addParameter('cloudfraccrit', 0.2);
p.addParameter('rowanomaly', 'XTrackFlags');
p.addParameter('DEBUG_LEVEL', 2);

p.parse(varargin{:});
pout = p.Results;

if ~iscell(start_date)
    start_date = {start_date};
end
if ~iscell(end_date)
    end_date = {end_date};
end

start_date = cellfun(@validate_date,start_date);
end_date = cellfun(@validate_date,end_date);
if numel(start_date) ~= numel(end_date)
    E.badinput('START_DATE and END_DATE must have the same number of elements');
end

if isa(pout.lon_bdy, 'GlobeGrid')
    input_grid = pout.lon_bdy;
    lon_bdy = input_grid.DomainLon;
    lat_bdy = input_grid.DomainLat;
elseif isnumeric(pout.lon_bdy) && numel(pout.lon_bdy) == 2
    lon_bdy = pout.lon_bdy;
    lat_bdy = pout.lat_bdy;
    if isempty(lat_bdy)
        E.badinput('If giving "lon_bdy" as a 2-element vector, "lat_bdy" must also be given the same way')
    elseif ~isnumeric(lat_bdy) || numel(lat_bdy) ~= 2
        E.badinput('"lat_bdy" must be a two-element numeric vector, if given');
    end
else
    E.badinput('A lat/lon domain must be specified either by giving an instance of GlobeGrid as the third input or by giving two two-element vectors as the third and fourth inputs');
end


if ~ischar(pout.behr_dir)
    E.badinput('"behr_dir" must be a string')
elseif ~exist(pout.behr_dir, 'dir');
    E.badinput('"behr_dir" %s does not exist', pout.behr_dir)
else
    behr_dir = pout.behr_dir;
end

if ~ischar(pout.filepattern)
    E.badinput('"filepattern" must be a string');
else
    file_pattern = pout.filepattern;
end

if ~ischar(pout.dayofweek)
    E.badinput('"dayofweek" must be a string');
else
    days_of_week = pout.dayofweek;
end

if ~ischar(pout.avgfield)
    E.badinput('"avgfield" must be a string')
else
    avg_field = pout.avgfield;
end

if ~ischar(pout.rejectmode)
    % Validation of what values it may have should be done in
    % omi_pixel_reject.
    E.badinput('"rejectmode" must be a string')
else
    reject_mode = pout.rejectmode;
end

% Likewise, these get checked in omi_pixel_reject
reject_details.cloud_type = pout.clouds;
reject_details.cloud_frac = pout.cloudfraccrit;
reject_details.row_anom_mode = pout.rowanomaly;
reject_details.check_behr_amf = true;


if ~isnumeric(pout.DEBUG_LEVEL) || ~isscalar(pout.DEBUG_LEVEL)
    E.badinput('"DEBUG_LEVEL" must be a scalar number')
else
    DEBUG_LEVEL = pout.DEBUG_LEVEL;
end

% Make sure the lesser boundary is first
lon_bdy = sort(lon_bdy);
lat_bdy = sort(lat_bdy);


% Convert 1 letter day abbreviations to a vector that isbusday can
% understand
day_abbrevs = {'U', 'M', 'T', 'W', 'R', 'F', 'S'};
weekend = true(size(day_abbrevs));
for a=1:numel(day_abbrevs)
    weekend(a) = isempty(strfind(upper(days_of_week), day_abbrevs{a}));
end

if (~islogical(pout.holidays) && ~isnumeric(pout.holidays)) || ~isscalar(pout.holidays)
    E.badinput('"holidays" must be a scalar logical or numeric value');
elseif pout.holidays
    holidays = [];
else
    holidays = 0;
end

if ~isscalar(pout.filterpsm) || (~islogical(pout.filterpsm) && ~isnumeric(pout.filterpsm))
    E.badinput('"filterpsm" must be a scalar logical or numeric value');
else
    filter_psm = pout.filterpsm;
end


%%%%%%%%%%%%%%%%%
% MAIN FUNCTION %
%%%%%%%%%%%%%%%%%

grid_var = 'OMI';

F = dir(fullfile(behr_dir, file_pattern));
if isempty(F)
    E.filenotfound('No files matching %s were found in %s', file_pattern, behr_dir);
end

first_file = true;

datevec = [];
for a=1:numel(start_date)
    datevec = cat(2, datevec, start_date(a):end_date(a));
end

for a=1:numel(F)
    file_date = get_file_date(F(a).name);
    if ~ismember(file_date, datevec)
        if DEBUG_LEVEL > 1
            fprintf('File %s outside of date range specified, skipping\n', F(a).name);
        end
        continue
    elseif ~isbusday(file_date, holidays, weekend)
        if DEBUG_LEVEL > 1
            fprintf('Skipping %s due to day of week specification (%s)\n', F(a).name, days_of_week);
        end
        continue
    end
    
    if DEBUG_LEVEL > 0
        fprintf('Adding data from %s\n', datestr(file_date));
    end
    
    D = load(fullfile(behr_dir, F(a).name), grid_var);
    if first_file
        if ~isfield(D.(grid_var), avg_field)
            E.badinput('%s is not a field in %s', avg_field, fullfile(behr_dir, F(a).name));
        end
        
        % The PSM gridded fields each have their own weights. All the CVM
        % gridded fields use the 'Areaweight' field. The PSM weights fields
        % are named as the field they are weights for + 'Weights', so we
        % check if that field is available and if so, use it as the
        % weighting field.
        weights_field = sprintf('%sWeights', avg_field);
        is_psm = true;
        if ~isfield(D.(grid_var), weights_field)
            weights_field = 'Areaweight';
            is_psm = false;
        end
        
        lon_grid = D.(grid_var)(1).Longitude;
        lat_grid = D.(grid_var)(1).Latitude;
        
        xx = lon_grid(1,:) >= lon_bdy(1) & lon_grid(1,:) <= lon_bdy(2);
        yy = lat_grid(:,1) >= lat_bdy(1) & lat_grid(:,1) <= lat_bdy(2);
        lon_grid = lon_grid(yy,xx);
        lat_grid = lat_grid(yy,xx);
        
        % Initialize the value and weights matrices
        sz = size(lon_grid);
        values = zeros(sz);
        weights = zeros(sz);
        
        first_file = false;
    end
    
    for b=1:numel(D.(grid_var))
        this_swath = D.(grid_var)(b);
        % PSM gridded fields have already been filtered for bad pixels,
        % including clouds, row anomaly, and VCD quality. CVM fields need to be
        % filtered, not because of the difference in the gridding algorithm,
        % but because I use the generic preprocessing function in PSM_Main for
        % all the CVM fields.
        if ~is_psm || filter_psm
            this_swath = omi_pixel_reject(this_swath, reject_mode, reject_details, 'weight_field', weights_field);
        end
        
        this_swath_values = this_swath.(avg_field)(yy,xx);
        this_swath_weights = this_swath.(weights_field)(yy,xx);
        values = nansum(cat(3, values, this_swath_values .* this_swath_weights), 3);
        weights = nansum(cat(3, weights, this_swath_weights),3);
    end
end

values = values ./ weights;

end

function file_date = get_file_date(file_name)
fdate_str = regexp(file_name, '\d\d\d\d\d\d\d\d', 'match', 'once');
file_date = datenum(fdate_str, 'yyyymmdd');
end

