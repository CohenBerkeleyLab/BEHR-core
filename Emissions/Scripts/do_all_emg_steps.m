% This will generate all the line densities and the statistics on them for
% the analysis of EMG fits to line densities derived from various a priori
% around atlanta.
%
% Generating the line densities and analyzing them are in separate sections
% if you already have the line densities.

%% calc wind speed/direction
timemode='hour'; %instant, hour, or avg
city = 'Montgomery'; %Atlanta, Birmingham, or Montgomery

fpath = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US_BEHR/hourly';
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

switch city
    case 'Atlanta'
        city_lon = -84.39; city_lat = 33.775;
    case 'Birmingham'
        city_lon = -86.80; city_lat = 33.52;
    case 'Montgomery'
        city_lon = -86.3; city_lat = 32.37;
    otherwise
        error('emg:city','City %s not recognized',city);
end

[windvel, theta] = misc_behr_wind_plots('calcavgwind', XLON, XLAT, Ubar, Vbar, city_lon, city_lat);
%% gen line densities
fprintf('Calculating line densities\n');
if ~exist('city','var')
    city = 'Atlanta';
end
switch city
    case 'Atlanta'
        wind_file = 'Atlanta-Wind-Conditions-1900UTC.mat';
    case 'Birmingham'
        wind_file = 'Birmingham-Wind-Conditions-1900UTC.mat';
    case 'Montgomery'
        wind_file = 'Montgomery-Wind-Conditions-1900UTC.mat';
    otherwise
        error('emg:city','City %s not recognized',city);
end

interp_bool = false;
box = [1 2 0.5 0.5]; %[0.5 2 0.5 0.5];

homedir = getenv('HOME');
behr_work_dir = fullfile(homedir,'Documents','MATLAB','BEHR','Workspaces','Wind speed');
hybrid_dir = 'SE US BEHR Hybrid - No ghost';
monthly_dir = 'SE US BEHR Monthly - No ghost';
coarse_mn_dir = 'SE US BEHR Monthly - No ghost - Coarse WRF';


% loads theta and windvel plus city lat and lon
load(fullfile(behr_work_dir,wind_file));

gtcrit = windvel >= 3.7;% & (theta < 0 | theta > 60);
ltcrit = windvel < 3.7;% & (theta < 0 | theta > 60);

F = dir(fullfile(behr_work_dir, hybrid_dir, 'OMI_BEHR_*.mat'));
fdnums = nan(size(F));
for a=1:numel(fdnums)
    [s,e] = regexp(F(a).name,'\d\d\d\d\d\d\d\d');
    fdnums(a) = datenum(F(a).name(s:e), 'yyyymmdd');
end
F = F(fdnums >= datenum('2013-06-01'));

fprintf('Fast hybrid\n')
[no2x_hygt5, no2ld_hygt5, no2ldstd_hygt5, lon_hygt5, lat_hygt5, no2cd_hygt5] = calc_line_density(fullfile(behr_work_dir, hybrid_dir),F,city_lon,city_lat,theta,'crit_logical',gtcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);
fprintf('Slow hybrid\n')
[no2x_hylt5, no2ld_hylt5, no2ldstd_hylt5, lon_hylt5, lat_hylt5, no2cd_hylt5] = calc_line_density(fullfile(behr_work_dir, hybrid_dir),F,city_lon,city_lat,theta,'crit_logical',ltcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);
fprintf('Fast monthly\n')
[no2x_mngt5, no2ld_mngt5, no2ldstd_mngt5, lon_mngt5, lat_mngt5, no2cd_mngt5] = calc_line_density(fullfile(behr_work_dir, monthly_dir),F,city_lon,city_lat,theta,'crit_logical',gtcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);
fprintf('Slow monthly\n')
[no2x_mnlt5, no2ld_mnlt5, no2ldstd_mnlt5, lon_mnlt5, lat_mnlt5, no2cd_mnlt5] = calc_line_density(fullfile(behr_work_dir, monthly_dir),F,city_lon,city_lat,theta,'crit_logical',ltcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);
fprintf('Fast coarse monthly\n')
[no2x_mn108gt5, no2ld_mn108gt5, no2ldstd_mn108gt5, lon_mn108gt5, lat_mn108gt5, no2cd_mn108gt5] = calc_line_density(fullfile(behr_work_dir, coarse_mn_dir),F,city_lon,city_lat,theta,'crit_logical',gtcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);
fprintf('Slow coarse monthly\n')
[no2x_mn108lt5, no2ld_mn108lt5, no2ldstd_mn108lt5, lon_mn108lt5, lat_mn180lt5, no2cd_mn108lt5] = calc_line_density(fullfile(behr_work_dir, coarse_mn_dir),F,city_lon,city_lat,theta,'crit_logical',ltcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);

%% Basic fitting
fprintf('Doing simple fitting of the line densities\n');
fmincon_out_level = 'none'; % could be 'none', 'final', 'iter'
fitfxn = 'ssresid'; % can be 'ssresid' or 'unexvar'
[f_hygt5.ffit, f_hygt5.emgfit, f_hygt5.stats, f_hygt5.f0, ~, f_hygt5.fitresults] = fit_line_density(no2x_hygt5, no2ld_hygt5, fmincon_out_level,'fittype',fitfxn);
[f_hylt5.ffit, f_hylt5.emgfit, f_hylt5.stats, f_hylt5.f0, ~, f_hylt5.fitresults] = fit_line_density(no2x_hylt5, no2ld_hylt5, fmincon_out_level,'fittype',fitfxn);
[f_mngt5.ffit, f_mngt5.emgfit, f_mngt5.stats, f_mngt5.f0, ~, f_mngt5.fitresults] = fit_line_density(no2x_mngt5, no2ld_mngt5, fmincon_out_level,'fittype',fitfxn);
[f_mnlt5.ffit, f_mnlt5.emgfit, f_mnlt5.stats, f_mnlt5.f0, ~, f_mnlt5.fitresults] = fit_line_density(no2x_mnlt5, no2ld_mnlt5, fmincon_out_level,'fittype',fitfxn);
[f_mn108gt5.ffit, f_mn108gt5.emgfit, f_mn108gt5.stats, f_mn108gt5.f0, ~, f_mn108gt5.fitresults] = fit_line_density(no2x_mn108gt5, no2ld_mn108gt5, fmincon_out_level,'fittype',fitfxn);
[f_mn108lt5.ffit, f_mn108lt5.emgfit, f_mn108lt5.stats, f_mn108lt5.f0, ~, f_mn108lt5.fitresults] = fit_line_density(no2x_mn108lt5, no2ld_mn108lt5, fmincon_out_level,'fittype',fitfxn);

%% Variational fitting
fprintf('Sampling the residuals using the standard deviations derived from the Hessian\n');
[Ffix_hygt5, F_hygt5] = fit_line_density_variation(no2x_hygt5, no2ld_hygt5, 30,'fittype',fitfxn);
[Ffix_hylt5, F_hylt5] = fit_line_density_variation(no2x_hylt5, no2ld_hylt5, 30,'fittype',fitfxn);
[Ffix_mngt5, F_mngt5] = fit_line_density_variation(no2x_mngt5, no2ld_mngt5, 30,'fittype',fitfxn);
[Ffix_mnlt5, F_mnlt5] = fit_line_density_variation(no2x_mnlt5, no2ld_mnlt5, 30,'fittype',fitfxn);
[Ffix_mn108gt5, F_mn108gt5] = fit_line_density_variation(no2x_mn108gt5, no2ld_mn108gt5, 30,'fittype',fitfxn);
[Ffix_mn108lt5, F_mn108lt5] = fit_line_density_variation(no2x_mn108lt5, no2ld_mn108lt5, 30,'fittype',fitfxn);

%% Montecarlo sampling
% Be warned - this'll take awhile.
fprintf('Now running Monte Carlo sampling of the fmincon fit function - be prepared to be bored for a bit.')
if ~exist('mc_hygt5','var')
    fprintf('\t Hybrid >= 5\n');
    [mc_hygt5.x, mc_hygt5.f0, mc_hygt5.mcacc] = emg_fit_montecarlo(no2x_hygt5, no2ld_hygt5);
    sd_mc_hygt5 = nan(1,5);
    for a=1:5; sd_mc_hygt5(a) = std(mc_hygt5.x(:,a)); end
end

if ~exist('mc_hylt5','var')
    fprintf('\t Hybrid < 5\n');
    [mc_hylt5.x, mc_hylt5.f0, mc_hylt5.mcacc] = emg_fit_montecarlo(no2x_hylt5, no2ld_hylt5);
    sd_mc_hylt5 = nan(1,5);
    for a=1:5; sd_mc_hylt5(a) = std(mc_hylt5.x(:,a)); end
end

if ~exist('mc_mngt5','var')
    fprintf('\t Monthly >= 5\n');
    [mc_mngt5.x, mc_mngt5.f0, mc_mngt5.mcacc] = emg_fit_montecarlo(no2x_mngt5, no2ld_mngt5);
    sd_mc_mngt5 = nan(1,5);
    for a=1:5; sd_mc_mngt5(a) = std(mc_mngt5.x(:,a)); end
end

if ~exist('mc_mnlt5','var')
    fprintf('\t Monthly < 5\n');
    [mc_mnlt5.x, mc_mnlt5.f0, mc_mnlt5.mcacc] = emg_fit_montecarlo(no2x_mnlt5, no2ld_mnlt5);
    sd_mc_mnlt5 = nan(1,5);
    for a=1:5; sd_mc_mnlt5(a) = std(mc_mnlt5.x(:,a)); end
end

if ~exist('mc_mn108gt5','var')
    fprintf('\t Coarse >= 5\n');
    [mc_mn108gt5.x, mc_mn108gt5.f0, mc_mn108gt5.mcacc] = emg_fit_montecarlo(no2x_mn108gt5, no2ld_mn108gt5);
    sd_mc_mn108gt5 = nan(1,5);
    for a=1:5; sd_mc_mn108gt5(a) = std(mc_mn108gt5.x(:,a)); end
end

if ~exist('mc_mn108lt5','var')
    fprintf('\t Coarse < 5\n');
    [mc_mn108lt5.x, mc_mn108lt5.f0, mc_mn108lt5.mcacc] = emg_fit_montecarlo(no2x_mn108lt5, no2ld_mn108lt5);
    sd_mc_mn108lt5 = nan(1,5);
    for a=1:5; sd_mc_mn108lt5(a) = std(mc_mn108lt5.x(:,a)); end
end

%% Variational sampling over the range of values from Monte Carlo
[Ffix_hygt5_mc, F_hygt5_mc] = fit_line_density_variation(no2x_hygt5, no2ld_hygt5, 30, 'sd', sd_mc_hygt5);
[Ffix_hylt5_mc, F_hylt5_mc] = fit_line_density_variation(no2x_hylt5, no2ld_hylt5, 30, 'sd', sd_mc_hylt5);
[Ffix_mngt5_mc, F_mngt5_mc] = fit_line_density_variation(no2x_mngt5, no2ld_mngt5, 30, 'sd', sd_mc_mngt5);
[Ffix_mnlt5_mc, F_mnlt5_mc] = fit_line_density_variation(no2x_mnlt5, no2ld_mnlt5, 30, 'sd', sd_mc_mnlt5);
[Ffix_mn108gt5_mc, F_mn108gt5_mc] = fit_line_density_variation(no2x_mn108gt5, no2ld_mn108gt5, 30, 'sd', sd_mc_mn108gt5);
[Ffix_mn108lt5_mc, F_mn108lt5_mc] = fit_line_density_variation(no2x_mn108lt5, no2ld_mn108lt5, 30, 'sd', sd_mc_mn108lt5);

%% Average correlation among fit params
corr_hygt5 = corrcov(inv(f_hygt5.fitresults.fminunc_hessian));
corr_hylt5 = corrcov(inv(f_hylt5.fitresults.fminunc_hessian));
corr_mngt5 = corrcov(inv(f_mngt5.fitresults.fminunc_hessian));
corr_mnlt5 = corrcov(inv(f_mnlt5.fitresults.fminunc_hessian));
corr_mn108gt5 = corrcov(inv(f_mn108gt5.fitresults.fminunc_hessian));
corr_mn108lt5 = corrcov(inv(f_mn108lt5.fitresults.fminunc_hessian));

avg_corr_mat = nanmean(cat(3,corr_hygt5,corr_hylt5,corr_mngt5,corr_mnlt5,corr_mn108gt5,corr_mn108lt5),3);

%% Montecarlo sampling starting at opt solution
% Be warned - this'll take awhile.
fprintf('Now running Monte Carlo sampling of the fmincon fit function - be prepared to be bored for a bit.')
if ~exist('mc_hygt5_fixf0','var')
    fprintf('\t Hybrid >= 5\n');
    fixf0 = [f_hygt5.ffit.a, f_hygt5.ffit.x_0, f_hygt5.ffit.mu_x, f_hygt5.ffit.sigma_x, f_hygt5.ffit.B]; 
    [mc_hygt5_fixf0.x, mc_hygt5_fixf0.f0, mc_hygt5_fixf0.mcacc] = emg_fit_montecarlo(no2x_hygt5, no2ld_hygt5, 'f0', fixf0);
    sd_mc_hygt5_fixf0 = nan(1,5);
    for a=1:5; sd_mc_hygt5_fixf0(a) = std(mc_hygt5_fixf0.x(:,a)); end
end

if ~exist('mc_hylt5_fixf0','var')
    fprintf('\t Hybrid < 5\n');
    fixf0 = [f_hylt5.ffit.a, f_hylt5.ffit.x_0, f_hylt5.ffit.mu_x, f_hylt5.ffit.sigma_x, f_hylt5.ffit.B];
    [mc_hylt5_fixf0.x, mc_hylt5_fixf0.f0, mc_hylt5_fixf0.mcacc] = emg_fit_montecarlo(no2x_hylt5, no2ld_hylt5, 'f0', fixf0);
    sd_mc_hylt5_fixf0 = nan(1,5);
    for a=1:5; sd_mc_hylt5_fixf0(a) = std(mc_hylt5_fixf0.x(:,a)); end
end

if ~exist('mc_mngt5_fixf0','var')
    fprintf('\t Monthly >= 5\n');
    fixf0 = [f_mngt5.ffit.a, f_mngt5.ffit.x_0, f_mngt5.ffit.mu_x, f_mngt5.ffit.sigma_x, f_mngt5.ffit.B];
    [mc_mngt5_fixf0.x, mc_mngt5_fixf0.f0, mc_mngt5_fixf0.mcacc] = emg_fit_montecarlo(no2x_mngt5, no2ld_mngt5, 'f0', fixf0);
    sd_mc_mngt5_fixf0 = nan(1,5);
    for a=1:5; sd_mc_mngt5_fixf0(a) = std(mc_mngt5_fixf0.x(:,a)); end
end

if ~exist('mc_mnlt5_fixf0','var')
    fprintf('\t Monthly < 5\n');
    fixf0 = [f_mnlt5.ffit.a, f_mnlt5.ffit.x_0, f_mnlt5.ffit.mu_x, f_mnlt5.ffit.sigma_x, f_mnlt5.ffit.B];
    [mc_mnlt5_fixf0.x, mc_mnlt5_fixf0.f0, mc_mnlt5_fixf0.mcacc] = emg_fit_montecarlo(no2x_mnlt5, no2ld_mnlt5, 'f0', fixf0);
    sd_mc_mnlt5_fixf0 = nan(1,5);
    for a=1:5; sd_mc_mnlt5_fixf0(a) = std(mc_mnlt5_fixf0.x(:,a)); end
end

if ~exist('mc_mn108gt5_fixf0','var')
    fprintf('\t Coarse >= 5\n');
    fixf0 = [f_mn108gt5.ffit.a, f_mn108gt5.ffit.x_0, f_mn108gt5.ffit.mu_x, f_mn108gt5.ffit.sigma_x, f_mn108gt5.ffit.B];
    [mc_mn108gt5_fixf0.x, mc_mn108gt5_fixf0.f0, mc_mn108gt5_fixf0.mcacc] = emg_fit_montecarlo(no2x_mn108gt5, no2ld_mn108gt5, 'f0', fixf0);
    sd_mc_mn108gt5_fixf0 = nan(1,5);
    for a=1:5; sd_mc_mn108gt5_fixf0(a) = std(mc_mn108gt5_fixf0.x(:,a)); end
end

if ~exist('mc_mn108lt5_fixf0','var')
    fprintf('\t Coarse < 5\n');
    fixf0 = [f_mn108lt5.ffit.a, f_mn108lt5.ffit.x_0, f_mn108lt5.ffit.mu_x, f_mn108lt5.ffit.sigma_x, f_mn108lt5.ffit.B];
    [mc_mn108lt5_fixf0.x, mc_mn108lt5_fixf0.f0, mc_mn108lt5_fixf0.mcacc] = emg_fit_montecarlo(no2x_mn108lt5, no2ld_mn108lt5, 'f0', fixf0);
    sd_mc_mn108lt5_fixf0 = nan(1,5);
    for a=1:5; sd_mc_mn108lt5_fixf0(a) = std(mc_mn108lt5_fixf0.x(:,a)); end
end
