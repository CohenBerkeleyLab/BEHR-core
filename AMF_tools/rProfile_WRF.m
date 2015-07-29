function [ no2_bins ] = rProfile_WRF( date_in, avg_mode, loncorns, latcorns, surfPres, pressures )
%RPROFILE_WRF Reads WRF NO2 profiles and averages them to pixels.
%   This function is the successor to rProfile_US and serves essentially
%   the same purpose - read in WRF-Chem NO2 profiles to use as the a priori
%   for the BEHR OMI NO2 retrieval algorithm. The key difference is that
%   this no longer assumes every NO2 profile from WRF is on the same
%   pressure grid. Therefore, it recognizes that the profiles will need to
%   be interpolated to a uniform pressure grid. It will therefore
%   interpolate all NO2 profiles to be averaged together for a pixel to the
%   same pressures before averaging, and will ensure that one bin below the
%   surface is calculated; that way omiAmfAK2 and integPr2 will be able to
%   calculate the best approximation of the surface concentration for
%   integration.
%
%   Further this will read the profiles directly from the netCDF files
%   containing the WRF output, processed by the slurmrun_wrf_output.sh
%   utility. It will look for WRF_BEHR files. These should have variables
%   no2, XLAT, XLONG, and pres. (pres is calculated as the sum of the
%   original WRF variables P and PB, perturbation and base pressure. See
%   calculated_quantities.nco and read_wrf_output.sh in the WRF_Utils
%   folder.)
%
%   This can also use different profiles: monthly, daily, or hourly. These
%   should be saved under folders so named in the main BEHR WRF output
%   folder.  The averaging will have been done in the processing by
%   slurmrun_wrf_output.sh, so the only difference to the execution of this
%   function will be that if hourly is chosen, it will need to select the
%   hour closest to 1400 local standard time for each pixel - this tries to
%   get as close to overpass as possible, although until the fraction of
%   the orbit passed is imported as well, an exact calculation of overpass
%   time will not be possible.
%
%   The inputs to this function are:
%       date_in: The date being processed, used to find the right file.
%
%       avg_mode: Should be 'hourly', 'daily', or 'monthly'. Chooses which
%       level of averaging to be used. See read_wrf_output.sh and
%       lonweight.nco in the WRF_Utils folder for the averaging process.
%
%       loncorns & latcorns: arrays containing the longitude and latitude
%       corners of the pixels. Can be 2- or 3-D, so long as the first
%       dimension has size 4 (i.e. loncorn(:,a) would give all 4 corners
%       for pixel a).
%
%       surfPres: a 1- or 2-D array containing the GLOBE surface pressures
%       for the pixels. This will be used to ensure enough bins are
%       calculated that omiAmfAK2 and integPr2 can interpolate to the
%       surface pressure.
%
%       pressure: the vector of pressures that the NO2 profiles should be
%       interpolated to.
%
%   Josh Laughner <joshlaugh5@gmail.com> 22 Jul 2015

DEBUG_LEVEL = 1;
E = JLLErrors;

% Defining custom errors
% Error for failing to find netCDF variable that is more descriptive about
% what likely went wrong. The error will take two additional arguments when
% called: the variable name and the file name.
E.addCustomError('ncvar_not_found','The variable %s is not defined in the file %s. Likely this file was not processed with (slurm)run_wrf_output.sh, or the processing failed before writing the calculated quantites.');

%%%%%%%%%%%%%%%%%%%%%
%%%%% CONSTANTS %%%%%
%%%%%%%%%%%%%%%%%%%%%

% The main folder for the WRF output, should contain subfolders 'monthly',
% 'daily', and 'hourly'. This will need modified esp. if trying to run on a
% PC.
wrf_output_path = fullfile('/Volumes','share2','USERS','LaughnerJ','WRF','SE_US_BEHR','NEI11Emis');


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

allowed_avg_modes = {'hourly','daily','monthly'};
avg_mode = lower(avg_mode);
if ~ischar(avg_mode) || ~ismember(avg_mode, allowed_avg_modes);
    E.badinput('avg_mode must be one of %s',strjoin(allowed_avg_modes,', '));
end



if size(loncorns,1) ~= 4 || size(latcorns,1) ~= 4
    E.badinput('loncorns and latcorns must have the corners along the first dimension (i.e. size(loncorns,1) == 4')
elseif ndims(loncorns) ~= ndims(latcorns) || ~all(size(loncorns) == size(latcorns)) 
    E.badinput('loncorns and latcorns must have the same dimensions')
end

sz_corners = size(loncorns);
sz_surfPres = size(surfPres);

% Check that the corner arrays and the surfPres array represent the same
% number of pixels
if ndims(loncorns)-1 ~= ndims(surfPres) || ~all(sz_corners(2:end) == sz_surfPres(1:end))
    E.badinput('The size of the surfPres array must be the same as the corner arrays without their first dimension (size(surfPres,1) == size(loncorns,2, etc)')
end

% pressures should just be a vector, monotonically decreasing
if ~isvector(pressures) || any(diff(pressures)>0)
    E.badinput('pressures must be a monotonically decreasing vector')
end

% Make sure that the input for the date can be understood by MATLAB as a
% date
try
    date_num_in = datenum(date_in);
catch err
    if strcmp(err.identifier,'MATLAB:datenum:ConvertDateString')
        E.badinput('%s could not be understood by MATLAB as a date')
    else
        rethrow(err);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD netCDF and READ VARIABLES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Redefine the path to the WRF data to include the averaging mode
wrf_output_path = fullfile(wrf_output_path,avg_mode);

% Find the file for this day or this month
year_in = year(date_num_in);
month_in = month(date_num_in);
day_in = day(date_num_in);
if strcmp(avg_mode,'monthly');
    file_pat = sprintf('WRF_BEHR_%s_%04d-%02d-*.nc', avg_mode, year_in, month_in);
else
    file_pat = sprintf('WRF_BEHR_%s_%04d-%02d-%02d.nc', avg_mode, year_in, month_in, day_in);
end

F = dir(fullfile(wrf_output_path,file_pat));
if numel(F) < 1
    E.filenotfound(file_pat);
elseif numel(F) > 1
    E.toomanyfiles(file_pat);
else
    wrf_info = ncinfo(fullfile(wrf_output_path,F(1).name));
end

wrf_vars = {wrf_info.Variables.Name};

% Load NO2 and check what unit it is - we'll use that to convert to
% parts-per-part later. Make the location of units as general as possible
try
    wrf_no2 = ncread(wrf_info.Filename, 'no2');
catch err
    if strcmp(err.identifier,'MATLAB:imagesci:netcdf:unknownLocation')
        E.callCustomError('ncvar_not_found','no2',F(1).name);
    else
        rethrow(err);
    end
end
ww = strcmp('no2',wrf_vars);
no2_attr = {wrf_info.Variables(ww).Attributes.Name};
aa = strcmp('units',no2_attr);
wrf_no2_units = wrf_info.Variables(ww).Attributes(aa).Value;

if isempty(wrf_no2_units)
    if DEBUG_LEVEL > 1; fprintf('\tWRF NO2 unit not identified. Assuming ppm\n'); end
    wrf_no2_units = 'ppm';
end

% Convert to be an unscaled mixing ratio (parts-per-part)
wrf_no2 = convert_units(wrf_no2, wrf_no2_units, 'ppp');

% Load the remaining variables
try
    varname = 'pres';
    wrf_pres = ncread(wrf_info.Filename, varname);
    varname = 'XLONG';
    wrf_lon = ncread(wrf_info.Filename, varname);
    varname = 'XLAT';
    wrf_lat = ncread(wrf_info.Filename, varname);
catch err
    if strcmp(err.identifier,'MATLAB:imagesci:netcdf:unknownLocation')
        E.callCustomError('ncvar_not_found',varname,F(1).name);
    else
        rethrow(err);
    end
end

% Finally, if the mode is "hourly" we need to take the correct WRF output
% time. We can use the utchr variable to do that.

if strcmp(avg_mode,'hourly')
    try 
        utchr = ncread(wrf_info.Filename, 'utchr');
    catch err
        if strcmp(err.identifier, 'MATLAB:imagesci:netcdf:unknownLocation')
            E.callCustomError('ncvar_not_found','utchr',F(1).name);
        else 
            rethrow(err)
        end
    end
    
    utc_offset = round(nanmean(loncorns(:))/15);
    % 14 - utc_offset will give 1400 local std. time in UTC, finding the
    % minimum between that and utchr indicates which WRF profile is closest
    % to overpass
    [~,uu] = min(abs(14 - utc_offset - utchr)); 
    % These two variables should have dimensions west_east, south_north,
    % bottom_top, Time
    wrf_no2 = wrf_no2(:,:,:,uu);
    wrf_pres = wrf_pres(:,:,:,uu);
    % These should have west_east, south_north, Time
    wrf_lon = wrf_lon(:,:,uu);
    wrf_lat = wrf_lat(:,:,uu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BIN PROFILES TO PIXELS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape the NO2 profiles and pressures such that the profiles are along
% the first dimension
num_profs = numel(wrf_lon);
prof_length = size(wrf_no2,3);

% Reorder dimensions. This will make the perm_vec be [3, 1, 2, (4:end)]
perm_vec = 1:ndims(wrf_no2);
perm_vec(3) = [];
perm_vec = [3, perm_vec];

wrf_no2 = permute(wrf_no2, perm_vec);
wrf_pres = permute(wrf_pres, perm_vec);

wrf_no2 = reshape(wrf_no2, prof_length, num_profs);
wrf_pres = reshape(wrf_pres, prof_length, num_profs);
wrf_lon = reshape(wrf_lon, 1, num_profs);
wrf_lat = reshape(wrf_lat, 1, num_profs);


num_pix = numel(surfPres);
no2_bins = nan(length(pressures), size(surfPres,1), size(surfPres,2));
for p=1:num_pix
    xall = loncorns(:,p);
    xall(5) = xall(1);
    
    yall = latcorns(:,p);
    yall(5) = yall(1);
    
    % Try to speed this up by removing profiles outside a rectangle around
    % the pixel first, then deal with the fact that the pixel is angled
    % relative to lat/lon.
    
    xx = wrf_lon < max(xall) & wrf_lon > min(xall) & wrf_lat < max(yall) & wrf_lat > min(yall);
    tmp_no2 = wrf_no2(:,xx);
    tmp_pres = wrf_pres(:,xx);
    tmp_lon = wrf_lon(xx);
    tmp_lat = wrf_lat(xx);
    
    yy = inpolygon(tmp_lon, tmp_lat, xall, yall);
    
    if sum(yy) < 1
        E.callError('no_prof','WRF Profile not found for pixel near %.1, %.1f',mean(xall),mean(yall));
    end
    
    tmp_no2(:,~yy) = [];
    tmp_pres(:,~yy) = [];
    
    % Interpolate all the NO2 profiles to the input pressures, then average
    % them together. Extrapolate so that later we can be sure to have one
    % bin below the surface pressure for omiAmfAK2 and integPr2.
    % Interpolate in log-log space to account for the exponential
    % dependence of pressure on altitude and the often exponential decrease
    % of concentration with altitude.
    
    interp_no2 = nan(length(pressures), size(tmp_no2,2));
    
    if ~iscolumn(pressures); pressures = pressures'; end
    
    for a=1:size(tmp_no2,2)
        interp_no2(:,a) = interp1(log(tmp_pres(:,a)), log(tmp_no2(:,a)), log(pressures), 'linear', 'extrap');
    end
    
    interp_no2 = exp(interp_no2);
    
    last_below_surf = find(pressures > surfPres(p),1,'last')-1;
    interp_no2(1:last_below_surf,:) = nan;
    
    no2_bins(:,p) = nanmean(interp_no2,2);
end




