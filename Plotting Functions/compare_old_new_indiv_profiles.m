function [  ] = compare_old_new_indiv_profiles( date_in, lon, lat, wrf_dir )
%COMPARE_OLD_NEW_PROFILES Compare v2 BEHR profiles with new WRF profiles
%   COMPARE_OLD_NEW_PROFILES( DATE_IN, LON, LAT, WRFDIR ) will plot NO2
%   profiles for the DATE_IN, interpolated to the given LON and LAT using
%   wrfout files located in WRFDIR and the standard monthly average
%   profiles for version 2 of BEHR. For the new WRF profiles, it chooses
%   the time based on the longitude (to calculate the UTC offset) and an
%   assumed OMI overpass time of 1330 local. LON and LAT can be vectors, in
%   which case one plot per element will be generated

E = JLLErrors;
if ~exist(wrf_dir,'dir')
    E.dir_dne(wrf_dir)
end

% Load old profiles
old_prof_file = sprintf('m%02d_NO2_profile.mat',month(date_in));
O = load(fullfile(BEHR_paths('no2_profile_path'), old_prof_file));
old_prof = O.PROFILE;
old_prof.NO2_profile = permute(old_prof.NO2_profile,[2 3 1]);

% Load new profiles
utc_offset = round(lon/15);
utc_hour = 13 - utc_offset; % we just need the hour, we'll always choose the file on the half hour if it exists
wrf_pattern = sprintf('wrfout*%s_%02d-*',datestr(date_in, 'yyyy-mm-dd'),utc_hour);
F = dir(fullfile(wrf_dir, wrf_pattern));
if isempty(F)
    E.filenotfound('File matching pattern %s in %s', wrf_pattern, wrf_dir);
end
wrf_file = glob({F.name},sprintf('.*%02d[:-]30[:-]00', utc_hour));
if isempty(wrf_file)
    warning('File for %02d:30 UTC not found, using %s instead', utc_hour, F(1).name);
    wrf_file = F(1).name;
end
wi = ncinfo(fullfile(wrf_dir, wrf_file{1}));

new_prof.Longitude = ncread(wi.Filename, 'XLONG');
new_prof.Latitude = ncread(wi.Filename, 'XLAT');
new_prof.Pressure = ncread(wi.Filename, 'pres');
new_prof.NO2_profile = ncread(wi.Filename, 'no2');



% Interpolate and plot. Old profiles all use the same pressures, so those
% don't need interpolated
for a=1:numel(lon)
    this_new_prof = interp_to_lat_lon(new_prof.Longitude, new_prof.Latitude, new_prof.NO2_profile, lon(a), lat(a));
    this_new_pres = interp_to_lat_lon(new_prof.Longitude, new_prof.Latitude, new_prof.Pressure, lon(a), lat(a));
    this_old_prof = interp_to_lat_lon(old_prof.Longitude, old_prof.Latitude, old_prof.NO2_profile, lon(a), lat(a));
    
    if all(isnan(this_new_prof)) 
        fprintf('No new profile available at long = %f, lat = %f\n', lon(a), lat(a));
        continue
    elseif all(isnan(this_old_prof))
        fprintf('No old profile available at long = %f, lat = %f\n', lon(a), lat(a));
        continue
    end
    
    figure; 
    plot(this_old_prof, old_prof.Pressure, this_new_prof, this_new_pres);
    set(gca,'ydir','reverse','fontsize',16)
    xlabel('[NO_2] (ppmv)')
    ylabel('Pressure (hPa)')
    title(sprintf('Long = %f  Lat = %f',lon(a), lat(a)));
end

end

function prof_out = interp_to_lat_lon(prof_lon, prof_lat, profs, lon, lat)
prof_out = nan(size(profs,3),1);
for a=1:size(profs,3)
    I = scatteredInterpolant(double(prof_lon(:)), double(prof_lat(:)), double(reshape(profs(:,:,a),[],1)));
    prof_out(a) = I(lon, lat);
end
end