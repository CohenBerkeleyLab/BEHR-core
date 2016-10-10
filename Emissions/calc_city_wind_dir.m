function [ windvel, theta, dnums, city_lon, city_lat ] = calc_city_wind_dir( city )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

timemode='hour'; %instant, hour, or avg

switch city
    case 'Atlanta'
        city_lon = -84.39; city_lat = 33.775;
        coast = 'E_US_BEHR';
    case 'Birmingham'
        city_lon = -86.80; city_lat = 33.52;
        coast = 'E_US_BEHR';
    case 'Montgomery'
        city_lon = -86.3; city_lat = 32.37;
        coast = 'E_US_BEHR';
    case 'SF'
        city_lon = -122.42; city_lat = 37.77;
        coast = 'W_US_BEHR';
    otherwise
        error('emg:city','City %s not recognized',city);
end

fpath = fullfile('/Volumes','share2','USERS','LaughnerJ','WRF',coast,'hourly');
F = dir(fullfile(fpath,'*.nc'));
dnums = [];
for a=1:numel(F)
    [s,e] = regexp(F(a).name, '\d\d\d\d-\d\d-\d\d');
    dnums = cat(1,dnums, datenum(F(a).name(s:e),'yyyy-mm-dd'));
end
% Remove any days before June 1st - allow WRF spinup
F = F(dnums >= datenum('2013-06-01'));
dnums = dnums(dnums >= datenum('2013-06-01'));
[XLON, XLAT, U, V, COSALPHA, SINALPHA, utchr] = read_wrf_vars(fpath, F, {'XLONG', 'XLAT', 'U', 'V', 'COSALPHA', 'SINALPHA', 'utchr'},false,0);
% We can take just one 2D slice of lon, lat, cos, and sin because these do
% not change in time. U and V we will take surface averaged over utchrs 19
% and 20 as these are closest to OMI overpass at Atlanta (13-14 local time)
utchr = double(utchr); % imported as int
for a=1:size(utchr,1)
    if any(utchr(a,:) ~= utchr(a,1))
        error('do_all_emg:utchr', 'Not all files have the same set of utc hours')
    end
end
if strcmpi(timemode,'instant')
    tt = find(utchr(:,1) == 19,1,'first');
elseif strcmpi(timemode,'hour')
    tt = utchr(:,1) == 19;
else
    tt = utchr(:,1) >= 19 & utchr(:,1) < 21;
end
XLON = XLON(:,:,1,1);
XLAT = XLAT(:,:,1,1);
COSALPHA = COSALPHA(:,:,1,1);
SINALPHA = SINALPHA(:,:,1,1);
% Lu 2015 used winds across the bottom 500 m or so
Ubar = squeeze(nanmean(nanmean(U(:,:,1:5,tt,:),3),4));
Vbar = squeeze(nanmean(nanmean(V(:,:,1:5,tt,:),3),4));
[Ue, Ve] = wrf_winds_transform(Ubar, Vbar, COSALPHA, SINALPHA);

[windvel, theta] = misc_behr_wind_plots('calcavgwind', XLON, XLAT, Ue, Ve, city_lon, city_lat);

end

