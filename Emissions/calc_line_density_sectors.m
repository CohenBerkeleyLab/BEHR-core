function [ no2_x, no2_linedens, no2_linedens_std, lon, lat, no2_mean, no2_std, num_valid_obs, nox, debug_cell ] = calc_line_density_sectors( fdates, center_lon, center_lat, theta, windvel, varargin )
%[ NO2_X, NO2_LINEDENS, NO2_LINEDENS_STD, LON, LAT, NO2_MEAN, NO2_STD, NUM_VALID_OBS] = CALC_LINE_DENSITY( FPATH, FNAMES, CENTER_LON, CENTER_LAT, THETA )
%   Calculate a wind-aligned line density for a given time period.
%
%   Calculates a line density of NO2 up and downwind of a city for 8
%   different wind sectors a la Beirle et al. 2011.
%
%   c.f.    Beirle et al., Science, 2015, pp. 1737-1739
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
%       'windvel' - a vector of wind velocities to be used in separated
%       days into slow and fast wind conditions. The separation value can
%       be altered with the parameter 'windsepcrit'.
%
%   Outputs:
%
%       no2_x - the x-coordinates of the line density (in km from center
%       lon/lat) as a structure for each direction.
%
%       no2_linedensity - the line density in mol/km.
%
%       no2_linedens_std - the standard deviation of the line density,
%       calculated by first calculated a standard deviation of the column
%       densities across the time average weighted by the grid cell
%       areaweight, then adding up those in quadrature along the line of
%       integration, multiplied by the width of the grid cells.
%
%       lon, lat - the longitude and latitude coordinates of the
%       wind-aligned column densities.
%
%       no2_mean - the average wind-aligned column densities.
%
%       no2_std - the area-weighted standard deviation of the wind aligned
%       column densities.
%
%       num_valid_obs - an array the same size as no2_mean that counts the
%       number of observations that went into each column density. NaNs and
%       cases where areaweight = 0 do not count.
%
%       nox - 3D matrix of individual rotated swaths, mainly for debugging
%       purposes.
%
%       debug_cell - a cell array with the file name and swath number
%       corresponding to each 2D slice of the output NOX.
%
%   Parameter inputs:
%
%
%       data_ind - which top-level indices of the Data structure loaded from
%       the files to use. Can be a single index, a vector of indicies, or
%       (as is default) an empty vector which will use all the swaths present
%       in each file.
%
%       'windsepcrit' - the speed used to separate slow and fast winds.
%       Defaults to 2 (m/s).
%
%       'rel_box_corners' - a four element vector to be passed to rotate
%       plume describing how large a box to use to circumscribe the plumes.
%       See rotate_plume for more information.
%
%       'force_calc' - if true will override the criteria that rejects days
%       with too many unfilled pixels and use all days that meet the
%       windvel criterion. Defaults to false.
%
%       'interp' - boolean, defaults to false. If true, individual days are
%       chosen if they have a sufficient number of viable observations to
%       represent the entire domain well. Any missing elements are filled
%       my interpolation.  If false, bad observations (cloud fraction or
%       row anomaly) are removed, but all days are averaged together. This
%       mode is closer to how Lu et al. 2015 did their analysis, I believe.
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
p.addParameter('windsepcrit',2);
p.addParameter('rel_box_corners',[2 4 2 2]);
p.addParameter('datatype','domino');
p.addParameter('griddingmethod','interp');
p.addParameter('n_rows_excl',0);
p.addParameter('DEBUG_LEVEL',1);

p.parse(varargin{:});

pout=p.Results;
nox_or_no2 = pout.nox_or_no2;
data_ind = pout.data_ind;
windsepcrit = pout.windsepcrit;
rel_box_corners = pout.rel_box_corners;
datatype = pout.datatype;
grid_method = pout.griddingmethod;
n_rows_excl = pout.n_rows_excl;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

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

if ~isnumeric(theta) || any(theta < -180 | theta > 180) || numel(theta) ~= numel(fdates)
    E.badinput('theta must be a numeric vector with values between -180 and +180 that has the same number of elements as the number of files to be loaded')
end

if ~isscalar(windsepcrit) || ~isnumeric(windsepcrit)
    E.badinput('windsepcrit must be a scalar number')
elseif windsepcrit <= 0
    warning('windsepcrit should be >0')
end

if ~isempty(rel_box_corners) && ( ~isvector(rel_box_corners) || numel(rel_box_corners) ~= 4 )
    E.badinput('rel_box_corners must be a 4 element vector');
end
if ~isscalar(DEBUG_LEVEL) || ~isnumeric(DEBUG_LEVEL) || DEBUG_LEVEL < 0
    E.badinput('DEBUG_LEVEL must be a scalar number >= 0.')
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

switch lower(datatype)
    case 'domino'
        Grid = design_dom_grid(0.18, rel_box_corners(2), center_lon, center_lat);
end

nox_fast = [];
nox_slow = [];

create_array = true;
i = 0;
for d=1:numel(fdates)
    D = load_files_for_day(fdates(d), datatype);
    
    if ~use_data_ind
        xx_swath = find_useful_swaths(D.Data, center_lon, center_lat, n_rows_excl);
        data_ind = find(xx_swath);
        n_swath = sum(xx_swath);
    else
        n_swath = numel(data_ind);
    end
    
    if d==1
        % Estimate how many swaths will be included in the average.
        estimated_ntimes = numel(fdates)*n_swath;
    end
    
    % We'll still "rotate" each day but divide it by sectors. This isn't so
    % much trying to reproduce Beirle 11 as give me a way to find out which
    % directions contribute to certain features of the line density.
    for s=1:n_swath
        i = i+1;
        e = data_ind(s);
        
        if strcmpi(grid_method,'rotate')
            E.notimplemented('grid_method == rotate (not updated to use NO2WindSector class');
            if DEBUG_LEVEL > 0; disp('Rotating plume'); end
            if ~isempty(rel_box_corners)
                OMI = rotate_plume(D.Data(e), center_lon, center_lat, theta(d), rel_box_corners);
            else
                OMI = rotate_plume(D.Data(e), center_lon, center_lat, theta(d));
            end
            if isempty(OMI.Longitude)
                if DEBUG_LEVEL > 0; fprintf('No grid cells in %s\n',datestr(fdates(d))); end
                continue
            end
            OMI = omi_pixel_reject(OMI,'omi',0.2,'XTrackFlags');
            xx = OMI.Areaweight > 0;
            if create_array
                directions = {'W','SW','S','SE','E','NE','N','NW'};
                theta_bin_edges = [-180, -157.5, -112.5, -67.5, -22.5, 22.5, 67.5, 112.5, 157.5];
                nox = make_empty_struct_from_cell(directions, nan(size(OMI.Longitude,1), size(OMI.Longitude,2), numel(fdates)*n_swath));
                aw = make_empty_struct_from_cell(directions, nan(size(OMI.Longitude,1), size(OMI.Longitude,2), numel(fdates)*n_swath));
                lon = OMI.Longitude;
                lat = OMI.Latitude;
                debug_cell = cell(numel(fdates)*n_swath,1);
                create_array = false;
            end
            
            debug_cell{i} = sprintf('%s: swath %d', datestr(fdates(d)), e);
            
            
            % Find which bin this belongs in
            if theta(d) < theta_bin_edges(2) || theta(d) > theta_bin_edges(end)
                bin = 'W';
            else
                theta_xx = theta_bin_edges(1:end-1) <= theta(d) & theta_bin_edges(2:end) > theta(d);
                bin = directions{theta_xx};
            end
            
            OMI.BEHRColumnAmountNO2Trop(~xx) = nan;
            OMI.Areaweight(~xx) = nan;
            nox.(bin)(:,:,i) = OMI.BEHRColumnAmountNO2Trop*nox_no2_scale;
            aw.(bin)(:,:,i) = OMI.Areaweight;
        else
            [nox_slow, nox_fast] = grid_by_interp(D.Data(e), estimated_ntimes, theta(d), windvel(d), windsepcrit, Grid, nox_slow, nox_fast);
        end
    end
end

if strcmpi(grid_method,'rotate')
    
    %no2_mean = nanmean(nox,3);
    for f=1:numel(directions)
        no2_mean.(directions{f}) = nansum2(nox.(directions{f}) .* aw.(directions{f}), 3) ./ nansum2(aw.(directions{f}), 3);
        num_valid_obs = sum( ~isnan(nox.(directions{f})) & aw.(directions{f}) > 0, 3);
        % Calculate the weighted standard deviation (c.f.
        % https://en.wikipedia.org/wiki/Mean_square_weighted_deviation)
        no2_var.(directions{f}) = (nansum2(aw.(directions{f}) .* nox.(directions{f}).^2, 3) .* nansum2(aw.(directions{f}),3) - (nansum2(aw.(directions{f}) .* nox.(directions{f}), 3)).^2)./(nansum2(aw.(directions{f}),3).^2 - nansum(aw.(directions{f}).^2,3));
        no2_std.(directions{f}) = sqrt(no2_var.(directions{f}));
        
        % Calculate the line density. See de Foy et al., Atmos. Environ. (2014) p.
        % 66. Basically an integration along the line perpendicular to the plume.
        
        if DEBUG_LEVEL > 0; disp('Calculating line density'); end
        no2_linedens.(directions{f}) = zeros(1,size(lon,2));
        no2_linedens_std.(directions{f}) = zeros(1,size(lon,2));
        no2_x.(directions{f}) = nan(1,size(lon,2));
        d_cm = nan(size(lon));
        for a=1:numel(no2_linedens.(directions{f}))
            for b=1:size(lon,1)-1 % because we need a latitudinal difference to average over, we can't do the last row.
                d_cm(b,a) = m_lldist(lon(b:b+1,a),lat(b:b+1,a)) * 1e5;
                if ~isnan(no2_mean.(directions{f})(b,a))
                    no2_linedens.(directions{f})(a) = no2_linedens.(directions{f})(a) + no2_mean.(directions{f})(b,a) * d_cm(b,a);
                    % add the uncertainties in quadrature.
                    no2_linedens_std.(directions{f})(a) = no2_linedens_std.(directions{f})(a) + no2_std.(directions{f})(b,a).^2 * d_cm(b,a);
                end
            end
            % Calculate x in km distant from the center lon/lat. OMI is a gridded
            % representation so OMI.Longitude(:,a) are all the same, as are
            % OMI.Latitude(b,:).
            no2_x.(directions{f})(a) = m_lldist([lon(1,a), center_lon], [center_lat, center_lat]) * sign(lon(1,a) - center_lon);
        end
        
        % Remove any values that never got anything added to them b/c there were no
        % non-nan values for that transect.
        rr = ~all(isnan(no2_mean.(directions{f})),1);
        no2_x.(directions{f}) = no2_x.(directions{f})(rr);
        no2_linedens.(directions{f}) = no2_linedens.(directions{f})(rr);
        no2_linedens_std.(directions{f}) = no2_linedens_std.(directions{f})(rr);
        
        % Finalize the uncertainties
        no2_linedens_std.(directions{f}) = sqrt(no2_linedens_std.(directions{f}));
        
        % Convert line density from molec/cm to moles/km
        no2_linedens.(directions{f}) = no2_linedens.(directions{f}) * 1e5 / 6.022e23;
        no2_linedens_std.(directions{f}) = no2_linedens_std.(directions{f}) * 1e5 / 6.022e23;
    end
else
    % Do the integration by interpolating to a rotated grid, then adding up
    % across that grid as before
    [ld_fast, ld_slow] = integrate_by_interp(nox_fast, nox_slow); 
end
end

function xx = find_useful_swaths(Data, clon, clat, n_edge_rows)
end_rows = [1+n_edge_rows, 60-n_edge_rows];
xx = false(size(Data));
for a=1:numel(Data)
    westlon = interp1(Data(a).Latitude(end_rows(1),:), Data(a).Longitude(end_rows(1),:), clat);
    eastlon = interp1(Data(a).Latitude(end_rows(2),:), Data(a).Longitude(end_rows(2),:), clat);
    % Note: this will fail for center lons around the international date
    % line
    if westlon < clon && eastlon > clon
        xx(a) = true;
    end
end
end

function D = load_files_for_day(dnum, datatype)
switch datatype
    case 'domino'
        D.Data = load_domino_for_day(dnum);
    otherwise
        E.badinput('%s is not a recognized value for the parameter ''datatype''', datatype);
end
end

function Data = load_domino_for_day(dnum)
domino_root = '/Volumes/share-sat/SAT/OMI/DOMINOv2.0';
year_str = datestr(dnum, 'yyyy');
month_str = datestr(dnum, 'mm');
f_str = sprintf('OMI-Aura_L2-OMDOMINO_%04dm%02d%02d*.he5', year(dnum), month(dnum), day(dnum));
F = dir(fullfile(domino_root, year_str, month_str, f_str));
fnames = cell(size(F));
for a=1:numel(F)
    fnames{a} = fullfile(domino_root, year_str, month_str, F(a).name);
end

dom_vars = {'Longitude', 'Latitude', 'TroposphericVerticalColumn', 'TroposphericColumnFlag', 'CloudFraction', 'SurfaceAlbedo'};

addpath('/Users/Josh/Documents/MATLAB/Non BEHR Satellite/OMI Utils');
Data = read_domino_simple(fnames{:}, 'variables', dom_vars);
for a=1:numel(Data)
    badvals = domino_pixel_reject(Data(a), 'cloudfrac', 0.3, 'nedgerows', 10);
    Data(a).TroposphericVerticalColumn(badvals) = nan;
end
end

function G = design_dom_grid(res, box_length, clon, clat)
% Ask for a domain that puts the center lon and lat in the center box and
% is long enough for the line densities requested. Add a number of grid
% cells around the edge to give us some wiggle room, just in case.
buffer = 5;
box_size = repmat(box_length + res*buffer, 1, 4);

G = GlobeGrid(res);
G.SetDomainByCenterWidth([clon, clat], box_size)
end


function [nox_slow, nox_fast] = grid_by_interp(Data, est_ntimes, theta, windvel, fast_slow_sep, Grid, nox_slow, nox_fast)
% The name of this function is perhaps misleading; this is a gridding
% approach that works with integrating the column densities to line
% densities later on in the code. The key point is that for each wind
% direction, NO2 column densities all around the city are stored; the
% integration domain will be extracted later.

E = JLLErrors;
% First, grid the NO2 data
no2_grid = nan(size(Grid.GridLon));
count = zeros(size(Grid.GridLon));
for a=1:numel(Data.Longitude)
    % Debugging
    if Data.Longitude(a) > min(Grid.GridLon(:)) && Data.Longitude(a) < max(Grid.GridLon(:)) && ...
            Data.Latitude(a) > min(Grid.GridLat(:)) && Data.Latitude(a) < max(Grid.GridLat(:)) && ...
            ~isnan(Data.TroposphericVerticalColumn(a))
        dum=1;
    end
    
    [indx, indy] = Grid.GridcellContains(Data.Longitude(a), Data.Latitude(a));
    if isempty(indx)
        continue
    end
    
    no2_grid(indx, indy) = nansum2([no2_grid(indx, indy) * count(indx, indy), Data.TroposphericVerticalColumn(a)]) / (count(indx, indy)+1);
    count(indx, indy) = count(indx, indy) + 1;
end

amclass = [isa(nox_slow, 'NO2WindSectors'), isa(nox_fast, 'NO2WindSectors')];
if any(amclass) && ~all(amclass)
    E.badinput('Do not input some but not all of the nox NO2WindSectors instances initialized')
elseif ~all(amclass)
    nox_slow = NO2WindSectors(Grid, est_ntimes);
    nox_fast = NO2WindSectors(Grid, est_ntimes);
end

if windvel < fast_slow_sep
    nox_slow.AddDataToDirection(no2_grid, theta, windvel);
else
    nox_fast.AddDataToDirection(no2_grid, theta, windvel);
end

end


function ld = integrate_by_interp(nox, box_edges)
% nox must be an instance of NO2WindSectors. This subfunction will take
% each direction in turn, interpolate it to a rotated grid, then integrate
% perpendicular to the wind direction

clon = nox.Grid.GridCenterLon;
clat = nox.Grid.GridCenterLat;
GridRot = GlobeGrid(nox.Grid.LonRes, 'projection', 'equirect-rotated');
GridRot.SetDomain([clon, clat], box_edges)

ld = make_empty_struct_from_cell(nox.directions, struct('no2_x',[],'no2_ld',[], 'no2_cd', 'no2_lon', 'no2_lat'));

for a=1:numel(nox.directions)
    GridRot.Rotation = nox.theta_bin_centers(a);
    no2_cd = interp2(nox.Grid.GridLon, nox.Grid.GridLat, nox.(nox.directions{a}), GridRot.GridLon, GridRot.GridLat);
    
    ld.(nox.directions{a}).no2_cd = no2_cd;
    ld.(nox.directions{a}).no2_lon = GridRot.GridLon;
    ld.(nox.directions{a}).no2_lat = GridRot.GridLat;
    
    ld.(nox.directions{a}).no2_x = do_along_wind_distance(GridRot.GridLon, GridRot.GridLat, clon, clat);
    ld.(nox.directions{a}).no2_ld = do_cross_wind_integration(GridRot.GridLoncorn, GridRot.GridLatcorn, no2_cd);
end

end


function ld = do_cross_wind_integration(loncorn, latcorn, no2_mean)
% Get the coordinates at the midpoints of the grid cells in the along wind
% direction. Keep the coordinates at the edges of the cell in the across
% track direction since we'll use that to get the length of the grid cell
% in that direction.
loncorn = (loncorn(:, 1:end-1) + loncorn(:, 2:end))/2;
latcorn = (latcorn(:, 1:end-1) + latcorn(:, 2:end))/2;

ld = nan(1, size(no2_mean,2));
d_cm = nan(size(no2_mean));
for a=1:size(no2_mean,2)
    for b=1:size(no2_mean,1)
        d_cm(b,a) = m_lldist(loncorn(b:b+1,a),latcorn(b:b+1,a)) * 1e5; 
        if ~isnan(no2_mean(b,a))
            ld(a) = ld(a) + no2_mean(b,a) * d_cm(b,a);
        end
    end
end
end

function x = do_along_wind_distance(lon, lat, center_lon, center_lat)

x = nan(1, size(lon,2));
mid = round(size(lon,1));

for a=1:size(lon,2)
    x(a) = m_lldist([lon(mid,a), center_lon], [lat(mid, a), center_lat]) * sign(lon(mid,a) - center_lon);
end
end
