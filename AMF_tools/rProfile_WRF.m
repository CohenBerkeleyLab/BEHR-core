function [ no2_bins, wrf_file ] = rProfile_WRF( date_in, profile_mode, loncorns, latcorns, omi_time, surfPres, pressures, wrf_output_path )
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
%       profile_mode: Should be 'daily' or 'monthly'. Chooses which WRF
%       profiles are used, the daily or monthly outputs.
%
%       loncorns & latcorns: arrays containing the longitude and latitude
%       corners of the pixels. Can be 2- or 3-D, so long as the first
%       dimension has size 4 (i.e. loncorn(:,a) would give all 4 corners
%       for pixel a).
%
%       omi_time: the starting time of the OMI swath in UTC. Used to match
%       up daily profiles to the OMI swath.
%
%       surfPres: a 1- or 2-D array containing the GLOBE surface pressures
%       for the pixels. This will be used to ensure enough bins are
%       calculated that omiAmfAK2 and integPr2 can interpolate to the
%       surface pressure.
%
%       pressure: the vector of pressures that the NO2 profiles should be
%       interpolated to.
%
%       wrf_output_path: optional, if provided, overrides which directory
%       this function will look for WRF output in. If passed an empty
%       string, the proper WRF directory will be found, just as if this
%       input was omitted.
%
%   Josh Laughner <joshlaugh5@gmail.com> 22 Jul 2015

DEBUG_LEVEL = 1;
E = JLLErrors;

% Defining custom errors
% Error for failing to find netCDF variable that is more descriptive about
% what likely went wrong. The error will take two additional arguments when
% called: the variable name and the file name.
E.addCustomError('ncvar_not_found','The variable %s is not defined in the file %s. Likely this file was not processed with (slurm)run_wrf_output.sh, or the processing failed before writing the calculated quantites.');

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% Get the WRF output path - this function will itself throw an error if the
% profile mode is wrong or the path does not exist.
if ~exist('wrf_output_path', 'var') || isempty(wrf_output_path)
    wrf_output_path = find_wrf_path(profile_mode, date_in);
else
    if ~ischar(wrf_output_path)
        E.badinput('WRF_OUTPUT_PATH must be a string');
    elseif ~exist(wrf_output_path, 'dir')
        E.dir_dne(wrf_output_path);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% LOAD netCDF and READ VARIABLES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[wrf_no2, wrf_pres, wrf_lon, wrf_lat, wrf_file] = load_wrf_vars();


num_profs = numel(wrf_lon);
prof_length = size(wrf_no2,3);

num_pix = numel(surfPres);
no2_bins = nan(length(pressures), size(surfPres,1), size(surfPres,2));

if any(size(wrf_lon) < 2) || any(size(wrf_lat) < 2)
    error('rProfile_WRF:wrf_dim','wrf_lon and wrf_lat should be 2D');
end
wrf_lon_bnds = [wrf_lon(1,1), wrf_lon(1,end), wrf_lon(end,end), wrf_lon(end,1)];
wrf_lat_bnds = [wrf_lat(1,1), wrf_lat(1,end), wrf_lat(end,end), wrf_lat(end,1)];
    

% If the WRF profiles are spaced at intervals larger than the smallest dimension
% of OMI pixels, interpolate instead of averaging b/c we will likely have at least
% some pixels with no profiles within them.

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BIN PROFILES TO PIXELS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape the NO2 profiles and pressures such that the profiles are along
% the first dimension

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


lons = squeeze(nanmean(loncorns,1));
lats = squeeze(nanmean(latcorns,1));
for p=1:num_pix
    if ~inpolygon(lons(p), lats(p), wrf_lon_bnds, wrf_lat_bnds)
        continue
    end

    no2_bins(:,p) = avg_apriori();
end

    function no2_vec = avg_apriori()
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
            %E.callError('no_prof','WRF Profile not found for pixel near %.1, %.1f',mean(xall),mean(yall));
            no2_vec = nan(length(pressures),1);
            return
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
        
        no2_vec = nanmean(interp_no2,2);
    end

    function [wrf_no2, wrf_pres, wrf_lon, wrf_lat, file_name] = load_wrf_vars()
        % Find the file for this day and the nearest hour May be "wrfout" or
        % "wrfout_subset"
        year_in = year(date_num_in);
        month_in = month(date_num_in);
        day_in = day(date_num_in);
        
        omi_utc_mean = omi_time_conv(nanmean(omi_time(:)));
        utc_hr = round(hour(omi_utc_mean));
        if strcmpi(profile_mode, 'daily')
            file_name = sprintf('wrfout_*_%04d-%02d-%02d_%02d-00-00', year_in, month_in, day_in, utc_hr);
        elseif strcmpi(profile_mode, 'monthly')
            file_name = sprintf('WRF_BEHR_monthly_%04d-%02d.nc', year_in, month_in);
        end
        
        F = dir(fullfile(wrf_output_path,file_name));
        if numel(F) < 1
            E.filenotfound(file_name);
        elseif numel(F) > 1
            E.toomanyfiles(file_name);
        else
            wrf_info = ncinfo(fullfile(wrf_output_path,F(1).name));
        end

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
        
        try
            wrf_no2_units = ncreadatt(wrf_info.Filename, 'no2', 'units');
        catch err
            % If we cannot find the attribute "units" in the file for some
            % reason
            if strcmp(err.identifier, 'MATLAB:imagesci:netcdf:libraryFailure')
                if DEBUG_LEVEL > 1; fprintf('\tWRF NO2 unit not identified. Assuming ppm\n'); end
                wrf_no2_units = 'ppm';
            else
                rethrow(err);
            end
        end
        
        % Convert to be an unscaled mixing ratio (parts-per-part)
        wrf_no2 = convert_units(wrf_no2, wrf_no2_units, 'ppp');
        
        wrf_vars = {wrf_info.Variables.Name};
        pres_precomputed = ismember('pres', wrf_vars);
        
        % Load the remaining variables
        try
            if pres_precomputed % P and PB are already combined into the 'pres' variable in the monthly files
                varname = 'pres';
                p_tmp = ncread(wrf_info.Filename, varname);
                pb_tmp = 0; % Allows us to skip a second logical test later
                p_units = ncreadatt(wrf_info.Filename, 'pres', 'units');
                pb_units = p_units;
            else
                varname = 'P';
                p_tmp = ncread(wrf_info.Filename, varname);
                p_units = ncreadatt(wrf_info.Filename, 'P', 'units');
                varname = 'PB';
                pb_tmp = ncread(wrf_info.Filename, varname);
                pb_units = ncreadatt(wrf_info.Filename, 'PB', 'units');
            end
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
        
        if ~strcmp(pb_units, p_units)
            E.callError('unit_mismatch', 'Units for P and PB in %s do not match', wrf_info.Filename);
        end
        wrf_pres = convert_units(p_tmp + pb_tmp, p_units, 'hPa'); % convert from units in the wrfout files to hPa
    end
end


