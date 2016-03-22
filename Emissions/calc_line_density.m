function [ no2_x, no2_linedens, lon, lat, no2_mean ] = calc_line_density( fpath, fnames, center_lon, center_lat, theta, varargin )
%CALC_LINE_DENSITY Calculate a wind-aligned line density for a given time period
%   Calculates a line density of NO2 up and downwind of a city by aligning
%   each day's plume to the x-axis as described in Valin 2013.  By fitting
%   an exponentially modified Gaussian function to the line density,
%   certain features of the NOx emissions and chemistry can be derived.
%
%   c.f.    Beirle et al., Science, 2015, pp. 1737-1739
%           Valin et al., Geophys. Res. Lett., 2013, pp. 1856-1860
%           de Foy et al., Atmos. Environ., 2014, pp. 66-77
%           Lu et al., ACP, 2015, pp. 10367-10383
%
%   Required inputs:
%
%       fpath - the path to where the files containing a Data structure are
%       located. This is the structure in files output from BEHR_main.m and
%       must be acceptable as the first input to rotate_plume.m
%
%       fnames - any one of several methods of specifying the files to
%       read. Could be a structure output from the dir() command, a
%       string that when passed as dir(fullfile(fpath,fnames)) will return
%       the file(s) desired, or a cell array of file names.
%
%       center_lon, center_lat - the coordinates considered the center of
%       the domain, usually the coordinate of a city or other point
%       emission source.
%
%       theta - a vector of wind directions given as degrees CCW from east
%       between -180 and 180. Must match the number of files to be loaded.
%
%   Parameter inputs:
%
%
%       data_ind - which top-level indices of the Data structure loaded from
%       the files to use. Can be a single index, a vector of indicies, or
%       (as is default) an empty vector which will use all the swaths present
%       in each file.
%
%       'windvel' - a vector of wind velocities to be used in filtering out
%       days that do not meet desired criteria. If not given, all days that
%       have sufficient coverage (pixels not removed for clouds or row
%       anomaly)
%
%       'windcrit' - the number to compare windvel values to, if given, the
%       parameter 'windop' must also be specified.
%
%       'windop' - can be '<', '>', '<=', or '>=' and will be used with
%       'windcrit' to evalute windvel values in the expression windvel
%       <windop> windcrit. If that evaluates to false, the day will be
%       skipped.
%
%       'crit_logical' - an alternative to the wind inputs, this should be a 
%       logical vector that is true for days that shold be included. Must
%       have the same number of elements as the number of files to use.
%
%       'rel_box_corners' - a four element vector to be passed to rotate
%       plume describing how large a box to use to circumscribe the plumes.
%       See rotate_plume for more information.
%
%       'force_calc' - if true will override the criteria that rejects days
%       with too many unfilled pixels and use all days that meet the
%       windvel criterion. Defaults to false.
%
%       'interp' - boolean, defaults to true. If true, individual days are
%       chosen if they have a sufficient number of viable observations to
%       represent the entire domain well. Any missing elements are filled
%       my interpolation.  If false, bad observations (cloud fraction or 
%       row anomaly) are removed, but all days are averaged together. This
%       mode is closer to how Lu et al. 2015 did their analysis, I believe.
%       In the future, the default may change to false.
%
%       'DEBUG_LEVEL' - level of output to console. Defaults to 1, 0 shuts
%       off all output.
%
%   Josh Laughner <joshlaugh5@gmail.com> 5 Feb 2016

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

p=inputParser;
p.addOptional('nox_or_no2','no2',@(x) ismember(lower(x),{'nox','no2'}));
p.addParameter('data_ind',[]);
p.addParameter('windvel',[]);
p.addParameter('windcrit',[]);
p.addParameter('windop','');
p.addParameter('crit_logical',[]);
p.addParameter('rel_box_corners',[]);
p.addParameter('force_calc',false);
p.addParameter('interp',true);
p.addParameter('DEBUG_LEVEL',1);

p.parse(varargin{:});

pout=p.Results;
nox_or_no2 = pout.nox_or_no2;
data_ind = pout.data_ind;
windvel = pout.windvel;
windcrit = pout.windcrit;
windop = pout.windop;
crit_logical = pout.crit_logical;
rel_box_corners = pout.rel_box_corners;
force_calc = pout.force_calc;
interp_bool = pout.interp;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

if ~ischar(fpath)
    E.badinput('fpath must be a string')
elseif ~exist(fpath,'dir')
    E.badinput('fpath (%s) does not exist',fpath);
end

if ischar(fnames)
    fnames_struct = dir(fullfile(fpath,fnames));
elseif iscellstr(fnames);
    fnames_struct = struct('name', repmat({''},size(fnames)));
    for f=1:numel(fnames)
        fnames_struct(f).name = fnames{f};
    end
elseif isstruct(fnames) && isfield(fnames,'name')
    fnames_struct = fnames;
else
    E.badinput('Input fnames is not a valid format; see documentation')
end

use_data_ind = true;
if ~isempty(data_ind) && ~isnumeric(data_ind)
    E.badinput('data_ind must be an empty matrix or numeric vector')
elseif isempty(data_ind)
    use_data_ind = false;
end
if ~isscalar(center_lon) || ~isnumeric(center_lon)
    E.badinput('center_lon must be a scalar number')
end
if ~isscalar(center_lat) || ~isnumeric(center_lat)
    E.badinput('center_lat must be a scalar number')
end

if ~isnumeric(theta) || any(theta < -180 | theta > 180) || numel(theta) ~= numel(fnames_struct)
    E.badinput('theta must be a numeric vector with values between -180 and +180 that has the same number of elements as the number of files to be loaded')
end

if ~isempty(crit_logical) && (numel(crit_logical) ~= numel(fnames_struct) || ~islogical(crit_logical))
    E.badinput('crit_logical must be a logical vector with the same number of elements as the number of files to be loaded');
end

% Check that if any one of windvel, windop, or windcrit are passed,
% all of them are. Alternatively, if the crit_logical matrix is
% available, none of those should be given.
E.addCustomError('windvelcrit','If any of windvel, windop, or windcrit are given, all must be given');
E.addCustomError('crit_conflict','crit_logical and the windvel/windop/windcrit parameters are mutually exclusive. You may only set one of them.');
windvel_set = false;
if ~isempty(windvel)
    if ~isnumeric(windvel) || any(windvel < 0) || numel(windvel) ~= numel(fnames_struct)
        E.badinput('windvel (if given) must be a numeric vector with values >= 0 that has the same number of elements as the number of files to be loaded')
    end
    windvel_set = true;
    if ~isempty(crit_logical)
        E.callCustomError('crit_conflict');
    end
end
if ~isempty(windcrit)
    if ~isnumeric(windcrit) || ~isscalar(windcrit) || windcrit < 0
        E.badinput('windcrit (if given) must be a numeric vector with values >= 0 that has the same number of elements as the number of files to be loaded')
    elseif ~windvel_set
        E.callCustomError('windvelcrit');
    end
    if ~isempty(crit_logical)
        E.callCustomError('crit_conflict');
    end
elseif windvel_set
    E.callCustomError('windvelcrit');
end
if ~isempty(windop)
    if ~ischar(windop) || ~ismember(windop,{'<','>','<=','>='})
        E.badinput('windop (if given) must be one of ''<'', ''>'', ''<='', ''>=''')
    elseif ~windvel_set
        E.callCustomError('windvelcrit');
    end
    if ~isempty(crit_logical)
        E.callCustomError('crit_conflict');
    end
elseif windvel_set
    E.callCustomError('windvelcrit')
end
% Finished with the wind crit input checking

if ~isempty(rel_box_corners) && ( ~isvector(rel_box_corners) || numel(rel_box_corners) ~= 4 )
    E.badinput('rel_box_corners must be a 4 element vector');
end
if ~isscalar(force_calc) || ~islogical(force_calc)
    E.badinput('force_calc must be a scalar logical');
end
if ~isscalar(DEBUG_LEVEL) || ~isnumeric(DEBUG_LEVEL) || DEBUG_LEVEL < 0
    E.badinput('DEBUG_LEVEL must be a scalar number >= 0.')
end

if ~isscalar(interp_bool) || ~islogical(interp_bool)
    E.badinput('interp must be a scalar logical')
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Constants
if strcmpi(nox_or_no2,'nox')
    nox_no2_scale = 1.32; % scales NO2 to NOx columns, c.f. supporting info for Beirle et al. 2011 (Science)
else
    nox_no2_scale = 1;
end

if isempty(crit_logical) && windvel_set
    crit_logical = eval(sprintf('windvel %s windcrit',windop));
end

create_array = true;
i = 0;
for d=1:numel(fnames_struct)
    D = load(fullfile(fpath,fnames_struct(d).name),'Data');
    
    if ~crit_logical(d)
        if DEBUG_LEVEL > 0; fprintf('Condition criteria not met, skipping %s\n',fnames_struct(d).name); end
        continue
    end
    
    if ~use_data_ind
        n_swath = numel(D.Data);
    else 
        n_swath = numel(data_ind);
    end

    for s=1:n_swath
        i = i+1;
        if ~use_data_ind
            e = s;
        else
            e = data_ind(s);
        end
        if DEBUG_LEVEL > 0; disp('Rotating plume'); end
        if ~isempty(rel_box_corners)
            OMI = rotate_plume(D.Data(e), center_lon, center_lat, theta(d), rel_box_corners);
        else
            OMI = rotate_plume(D.Data(e), center_lon, center_lat, theta(d));
        end
        if isempty(OMI.Longitude)
            if DEBUG_LEVEL > 0; fprintf('No grid cells in %s\n',fnames_struct(d).name); end 
            continue
        end
        OMI = omi_pixel_reject(OMI,'omi',0.2,'XTrackFlags');
        xx = OMI.Areaweight > 0;
        if create_array
            nox = nan(size(OMI.Longitude,1), size(OMI.Longitude,2), numel(fnames_struct)*n_swath);
            aw = nan(size(OMI.Longitude,1), size(OMI.Longitude,2), numel(fnames_struct)*n_swath);
            lon = OMI.Longitude;
            lat = OMI.Latitude;
            create_array = false;
        end

        if interp_bool
            % This criterion accounts for how many neighbors are empty, giving more
            % weight to large clumps of NaNs (due to row anomaly or clouds) and
            % less weight to scattered pixels. It still looks like 50% is a good
            % cutoff.
            mfrac = badpix_metric(~xx);
            if mfrac > .5 && ~force_calc
                if DEBUG_LEVEL > 0; fprintf('Too many clumped missing pixels, skipping %s\n',fnames_struct(d).name); end
                continue
            end
            
            % Fill in empty grid boxes by interpolation to prevent discontinuity in the
            % line density due to an uneven number of missing values for different
            % distances from the city.
            if DEBUG_LEVEL > 0; disp('Interpolating to fill in gaps in NO2 matrix'); end
            F = scatteredInterpolant(OMI.Longitude(xx), OMI.Latitude(xx), OMI.BEHRColumnAmountNO2Trop(xx));
            Faw = scatteredInterpolant(OMI.Longitude(xx), OMI.Latitude(xx), OMI.Areaweight(xx));
            nox(:,:,i) = F(OMI.Longitude, OMI.Latitude)*nox_no2_scale;
            aw(:,:,i) = F(OMI.Longitude, OMI.Latitude);
        else
            
            OMI.BEHRColumnAmountNO2Trop(~xx) = nan;
            OMI.Areaweight(~xx) = nan;
            nox(:,:,i) = OMI.BEHRColumnAmountNO2Trop*nox_no2_scale;
            aw(:,:,i) = OMI.Areaweight;
        end
    end
end

%no2_mean = nanmean(nox,3);
no2_mean = nansum2(nox .* aw, 3) ./ nansum2(aw, 3);

% Calculate the line density. See de Foy et al., Atmos. Environ. (2014) p.
% 66. Basically an integration along the line perpendicular to the plume.

if DEBUG_LEVEL > 0; disp('Calculating line density'); end
no2_linedens = zeros(1,size(lon,2));
no2_x = nan(1,size(lon,2));
d_cm = nan(size(lon));
for a=1:numel(no2_linedens)
    for b=1:size(lon,1)-1 % because we need a latitudinal difference to average over, we can't do the last row.
        d_cm(b,a) = m_lldist(lon(b:b+1,a),lat(b:b+1,a)) * 1e5;
        if ~isnan(no2_mean(b,a))
            no2_linedens(a) = no2_linedens(a) + no2_mean(b,a) * d_cm(b,a);
        end
    end
    % Calculate x in km distant from the center lon/lat. OMI is a gridded
    % representation so OMI.Longitude(:,a) are all the same, as are
    % OMI.Latitude(b,:).
    no2_x(a) = m_lldist([lon(1,a), center_lon], [center_lat, center_lat]) * sign(lon(1,a) - center_lon);
end

% Remove any values that never got anything added to them b/c there were no
% non-nan values for that transect.
rr = ~all(isnan(no2_mean),1);
no2_x = no2_x(rr);
no2_linedens = no2_linedens(rr);

% Convert line density from molec/cm to moles/km
no2_linedens = no2_linedens * 1e5 / 6.022e23;

end

function [mfrac, msum] = badpix_metric(xx)
% Calculates a metric for the badness of missing pixels. A missing grid
% cell doesn't contribute to this metric unless it has at least one
% neighbor that is also missing
sz = size(xx);
m = zeros(sz);
for a=1:sz(1)
    for b=1:sz(2)
        if a > 1 && b > 1
            m(a,b) = m(a,b) + xx(a-1, b-1);
        end
        if a > 1
            m(a,b) = m(a,b) + xx(a-1,b);
        end
        if a > 1 && b < sz(2)
            m(a,b) = m(a,b) + xx(a-1,b+1);
        end
        if b < sz(2)
            m(a,b) = m(a,b) + xx(a,b+1);
        end
        if a < sz(1) && b < sz(2)
            m(a,b) = m(a,b) + xx(a+1,b+1);
        end
        if a < sz(1)
            m(a,b) = m(a,b) + xx(a+1,b);
        end
        if a < sz(1) && b > 1
            m(a,b) = m(a,b) + xx(a+1,b-1);
        end
        if b > 1
            m(a,b) = m(a,b) + xx(a,b-1);
        end
    end
end

% How many neighbors there are, i.e. the maximum value m could take on.
% Corner cells have 3 neighbors, non-corner edge cells have 5 and non-edge
% cells have 8.
n_neighbors = prod(sz-2)*8 + (prod(sz) - prod(sz-2) - 4) * 5 + 12;
msum = sum(m(:));
mfrac = msum / n_neighbors;
end
