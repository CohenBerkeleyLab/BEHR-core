function [ no2_bins ] = rProfile_WRF( date_in, hour_in, avg_mode, loncorns, latcorns, surfPres, pressures )
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
wrf_output_path = fullfile('/Volumes','share2','USERS','LaughnerJ','WRF','SE_US_TEMPO','NEI11Emis');


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

allowed_avg_modes = {'hourly','monthly','hybrid'};
avg_mode = lower(avg_mode);
if ~ischar(avg_mode) || ~ismember(avg_mode, allowed_avg_modes);
    E.badinput('avg_mode must be one of %s',strjoin(allowed_avg_modes,', '));
end

if ~isscalar(hour_in) || ~isnumeric(hour_in)
    E.badinput('hour_in must be a numeric scalar')
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

if strcmpi(avg_mode,'hybrid')
    [wrf_no2_m, wrf_pres_m, wrf_lon, wrf_lat] = load_wrf_vars('monthly');
    [wrf_no2_h, wrf_pres] = load_wrf_vars('hourly');
    wrf_no2 = combine_wrf_profiles(wrf_no2_h, wrf_pres, wrf_no2_m, wrf_pres_m);
else
    [wrf_no2, wrf_pres, wrf_lon, wrf_lat] = load_wrf_vars(avg_mode);
end
    
% Grid to coarser resolution if requested
if coarse > 1
    if ndims(wrf_no2) > 3
        E.notimplemented('coarse regridding if NO2 matrix has >1 time coordinate')
    end
    coarse_size = ceil(size(wrf_lon)/coarse);
    wrf_no2_c = nan([coarse_size, size(wrf_no2,3)]);
    wrf_pres_c = nan([coarse_size, size(wrf_pres,3)]);
    wrf_lon_c = nan(coarse_size);
    wrf_lat_c = nan(coarse_size);
    
    for a=1:size(wrf_lon_c,1)
        for b=1:size(wrf_lon_c,2)
            a1 = (a-1)*coarse + 1;
            b1 = (b-1)*coarse + 1;
            a2 = min(a1+coarse-1, size(wrf_lon,1));
            b2 = min(b1+coarse-1, size(wrf_lon,2));
            
            num_profs = (a2-a1+1)*(b2-b1+1);
            
            no2_tmp = reshape(wrf_no2(a1:a2,b1:b2,:),num_profs,[]);
            wrf_no2_c(a,b,:) = nanmean(no2_tmp,1);
            pres_tmp = reshape(wrf_pres(a1:a2,b1:b2,:),num_profs,[]);
            wrf_pres_c(a,b,:) = nanmean(pres_tmp,1);
            wrf_lon_c(a,b) = nanmean(reshape(wrf_lon(a1:a2,b1:b2),num_profs,[]));
            wrf_lat_c(a,b) = nanmean(reshape(wrf_lat(a1:a2,b1:b2),num_profs,[]));
        end
    end
end

wrf_no2 = wrf_no2_c;
wrf_pres = wrf_pres_c;
wrf_lon = wrf_lon_c;
wrf_lat = wrf_lat_c;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BIN PROFILES TO PIXELS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate corner points for the WRF grid cells.
%[wrf_loncorn, wrf_latcorn] = wrf_corners(wrf_lon, wrf_lat);

% Reshape the NO2 profiles and pressures such that the profiles are along
% the first dimension
num_profs = numel(wrf_lon);
prof_length = size(wrf_no2,3);


if interp_bool
    wrf_lon = repmat(wrf_lon,1,1,prof_length);
    wrf_lat = repmat(wrf_lat,1,1,prof_length);
    wrf_z = repmat(permute(1:prof_length,[1 3 2]),size(wrf_lon,1),size(wrf_lon,2),1);
    wrf_no2_F = scatteredInterpolant(wrf_lon(:), wrf_lat(:), wrf_z(:), wrf_no2(:));
    wrf_pres_F = scatteredInterpolant(wrf_lon(:), wrf_lat(:), wrf_z(:), wrf_pres(:));
else
    % Reorder dimensions. This will make the perm_vec be [3, 1, 2, (4:end)]
    perm_vec = 1:ndims(wrf_no2);
    perm_vec(3) = [];
    perm_vec = [3, perm_vec];
    
    wrf_no2 = permute(wrf_no2, perm_vec);
    wrf_pres = permute(wrf_pres, perm_vec);

    wrf_no2 = reshape(wrf_no2, prof_length, num_profs);
    wrf_pres = reshape(wrf_pres, prof_length, num_profs);
    wrf_lon = reshape(wrf_lon, num_profs, 1);
    wrf_lat = reshape(wrf_lat, num_profs, 1);
end

    no2_bins = nan(length(pressures), size(surfPres,1), size(surfPres,2));


% For each pixel, just find the closest WRF profile. Error if the distance
% to the profile is greater than about the pixel resolution. Will also
% implement interpolation to compare results.
for x = 1:size(lons,1)
    for y = 1:size(lons,2)
        if ~interp_bool
            [~,I] = min(sqrt((lons(x,y)-wrf_lon).^2 + (lats(x,y)-wrf_lat).^2));
            interp_no2 = exp(interp1(log(wrf_pres(:,I)), log(wrf_no2(:,I)), log(pressures),'linear','extrap'));
        else
            interp_no2 = nan(prof_length,1);
            interp_pres = nan(prof_length,1);
            for i=1:prof_length
                interp_no2(i) = wrf_no2_F(lons(x,y), lats(x,y),i);
                interp_pres(i) = wrf_pres_F(lons(x,y), lats(x,y),i);   
            end
            interp_no2 = exp(interp1(log(interp_pres), log(interp_no2), log(pressures),'linear','extrap'));
        end

        last_below_surf = find(pressures > surfPres(x,y),1,'last');
        no2_bins(last_below_surf:end,x,y) = interp_no2(last_below_surf:end);
    end
end

    function [wrf_no2, wrf_pres, wrf_lon, wrf_lat] = load_wrf_vars(avg_mode)
        % Redefine the path to the WRF data to include the averaging mode
        wrf_output_mode_path = fullfile(wrf_output_path,avg_mode);
        
        % Find the file for this day or this month
        year_in = year(date_num_in);
        month_in = month(date_num_in);
        day_in = day(date_num_in);
        if strcmp(avg_mode,'monthly');
            file_pat = sprintf('WRF_TEMPO_%s_%04d-%02d.nc', avg_mode, year_in, month_in);
        else
            file_pat = sprintf('WRF_TEMPO_%s_%04d-%02d-%02d.nc', avg_mode, year_in, month_in, day_in);
        end
        
        F = dir(fullfile(wrf_output_mode_path,file_pat));
        if numel(F) < 1
            E.filenotfound(file_pat);
        elseif numel(F) > 1
            E.toomanyfiles(file_pat);
        else
            wrf_info = ncinfo(fullfile(wrf_output_mode_path,F(1).name));
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
        
        % Finally, with either the monthly or hourly files (daily makes no
        % sense for TEMPO, why would we average over the course of a day
        % when the satellite takes repeated measurements?) we need to cut
        % down to the proper hour for this scan, but the way the averaging
        % occurs puts the hour along a different dimension in the two
        % files.
        try
            utchr = ncread(wrf_info.Filename, 'utchr');
        catch err
            if strcmp(err.identifier, 'MATLAB:imagesci:netcdf:unknownLocation')
                E.callCustomError('ncvar_not_found','utchr',F(1).name);
            else
                rethrow(err)
            end
        end
        
        uu = utchr == hour_in;
        if sum(uu) < 1
            E.callError('hour_not_avail','The hour %d is not available in the WRF chem output file %s.', hour_in, F(1).name);
        end
            
        if strcmp(avg_mode,'hourly')
            % These two variables should have dimensions west_east, south_north,
            % bottom_top, Time (i.e. day), and hour_index
            wrf_no2 = squeeze(wrf_no2(:,:,:,:,uu));
            wrf_pres = squeeze(wrf_pres(:,:,:,:,uu));
            % These should have west_east, south_north, Time
            wrf_lon = wrf_lon(:,:,1);
            wrf_lat = wrf_lat(:,:,1);
        elseif strcmp(avg_mode,'monthly')
            % In the monthly files these will not have the Time dimension
            wrf_no2 = squeeze(wrf_no2(:,:,:,uu));
            wrf_pres = squeeze(wrf_pres(:,:,:,uu));
            % the lon/lat coordinates are already cut down to the right
            % number of dimensions.
        end
    end
end

function wrf_no2 = combine_wrf_profiles(wrf_no2_h, wrf_pres_h, wrf_no2_m, ~)
% For now this will be very simple, it will just replace any part of the
% hourly profile above 750 hPa with the monthly profile.  I chose 750 hPa
% to start with because that seems to be the pressure where the NO2
% profiles around Atlanta get into free troposphere NO2, with less
% influence from surface winds (i.e. no longer an exponential decay with
% altitude). This is purely qualitative and intended (currently) to test if
% removing FT profile changes from the AMF calc will simplify the
% relationship between AMF and wind around Atlanta.
%
% In the future, the better way might be to use the monthly average PBL
% height from WRF-Chem..
E = JLLErrors;

if any(size(wrf_no2_h)~=size(wrf_no2_m))
    E.sizeMismatch('wrf_no2_h','wrf_no2_m');
elseif any(size(wrf_no2_h)~=size(wrf_pres_h))
    E.sizeMismatch('wrf_no2_h','wrf_pres_h');
end

pp = wrf_pres_h < 750;
wrf_no2 = wrf_no2_h;
wrf_no2(pp) = wrf_no2_m(pp);

end


