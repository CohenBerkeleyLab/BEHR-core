function [  ] = compare_old_new_profiles( date_in, wrf_dir )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

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
user_ans = ask_multichoice('Daily or monthly profiles?', {'daily','monthly'}, 'list', true);

if strcmpi(user_ans, 'daily')
    utc_hour = ask_number('Enter the UTC hour to pick the new profiles from', 'testfxn', @(x) x >= 0 && x <= 23 && mod(x,1) == 0, 'testmsg', 'Value must be a whole number between 0 and 23');
    wrf_pattern = sprintf('wrfout*%s_%02d-*',datestr(date_in, 'yyyy-mm-dd'),utc_hour);
    
    F = dir(fullfile(wrf_dir, wrf_pattern));
    if isempty(F)
        E.filenotfound('File matching pattern %s in %s', wrf_pattern, wrf_dir);
    end
    wrf_file = glob({F.name},sprintf('.*%02d[:-]30[:-]00', utc_hour));
    if isempty(wrf_file)
        warning('File for %02d:30 UTC not found, using %s instead', utc_hour, F(1).name);
        wrf_file = {F(1).name};
    end
else
    wrf_pattern = sprintf('WRF_BEHR_monthly_%s-*.nc', datestr(date_in, 'yyyy-mm'));
    F = dir(fullfile(wrf_dir, wrf_pattern));
    if isempty(F)
        E.filenotfound('File matching pattern %s in %s', wrf_pattern, wrf_dir);
    elseif numel(F) > 1
        E.toomanyfiles('Too many files found matching pattern %s in %s', wrf_pattern, wrf_dir);
    end
    wrf_file = {F(1).name};
end
wi = ncinfo(fullfile(wrf_dir, wrf_file{1}));


new_prof.Longitude = double(ncread(wi.Filename, 'XLONG'));
new_prof.Latitude = double(ncread(wi.Filename, 'XLAT'));
new_prof.Pressure = double(ncread(wi.Filename, 'pres'));
new_prof.NO2_profile = double(ncread(wi.Filename, 'no2'));


% Plot the old profiles without interpolating
[slon, slat] = state_outlines('not', 'ak', 'hi');
plot_slice_gui(old_prof.NO2_profile, old_prof.Longitude, old_prof.Latitude, slon, slat, 'Old profiles, no interp');

% Interpolate the old profiles to the new coordinates
old_prof_interp = interp_profiles(old_prof, new_prof);
plot_slice_gui(old_prof_interp.NO2_profile, old_prof_interp.Longitude, old_prof_interp.Latitude, slon, slat, 'Old profiles, after interp');

% Plot the new profiles and the differences
plot_slice_gui(new_prof.NO2_profile, new_prof.Longitude, new_prof.Latitude, slon, slat, 'New profiles');
plot_slice_gui(new_prof.NO2_profile - old_prof_interp.NO2_profile, old_prof_interp.Longitude, old_prof_interp.Latitude, slon, slat, 'New profiles - old interp profiles');

if strcmpi(ask_multichoice('Add new_prof and old_prof_interp to workspace?',{'y','n'}),'y');
    putvar(new_prof, old_prof_interp);
end
end

function op_interp = interp_profiles(old_prof, new_prof)
op_interp.Longitude = new_prof.Longitude;
op_interp.Latitude = new_prof.Latitude;
op_interp.Pressure = new_prof.Pressure;

lonv = old_prof.Longitude;
latv = old_prof.Latitude;
presv = old_prof.Pressure;

lonq = new_prof.Longitude;
latq = new_prof.Latitude;
presq = new_prof.Pressure;

sz = size(new_prof.NO2_profile);
op_no2_tmp = nan(sz(1), sz(2), size(old_prof.NO2_profile,3));

% First interpolate each level to lon/lat coordinates
for a=1:size(old_prof.NO2_profile,3)
    fprintf('Interpolating to lat/lon - level %d\n', a)
    F = scatteredInterpolant(lonv(:), latv(:), reshape(old_prof.NO2_profile(:,:,a),[],1));
    op_no2_tmp(:,:,a) = F(lonq, latq);
end

% Now interpolate in log-log space for the pressure coordinate
op_interp.NO2_profile = nan(size(new_prof.NO2_profile));
% Ensure no negative values exist for the log
op_no2_tmp(op_no2_tmp<0) = nan;
p_vec = log(presv(:));
for a=1:prod(sz(1:2))
    fprintf('Vertical interp: %d of %d\n', a, prod(sz));
    [x,y] = ind2sub(sz(1:2),a);
    no2_vec = log(squeeze(op_no2_tmp(x,y,:)));
    pq_vec = log(squeeze(presq(x,y,:)));
    
    no2_vec_tmp = interp1(p_vec, no2_vec, pq_vec);
    op_interp.NO2_profile(x,y,:) = exp(no2_vec_tmp);
end


end