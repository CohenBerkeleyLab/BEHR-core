% This will generate all the line densities and the statistics on them for
% the analysis of EMG fits to line densities derived from various a priori
% around atlanta.
%
% Generating the line densities and analyzing them are in separate sections
% if you already have the line densities.

%% gen line densities
fprintf('Calculating line densities\n');
city_lon = -84.39;
city_lat = 33.755;

homedir = getenv('HOME');
behr_work_dir = fullfile(homedir,'Documents','MATLAB','BEHR','Workspaces','Wind speed');
hybrid_dir = 'SE US BEHR Hybrid - No ghost';
monthly_dir = 'SE US BEHR Monthly - No ghost';
coarse_mn_dir = 'SE US BEHR Monthly - No ghost - 288 km resolution';

% loads theta and windvel
load(fullfile(behr_work_dir,'Atlanta-Wind-Conditions.mat'));

[no2x_hygt5, no2ld_hygt5] = calc_line_density(fullfile(behr_work_dir, hybrid_dir),'OMI_BEHR_*.mat',2,city_lon,city_lat,theta,'windvel',windvel,'windcrit',5,'windop','>=','rel_box_corners',[0.5 2 0.5 0.5]);
[no2x_hylt5, no2ld_hylt5] = calc_line_density(fullfile(behr_work_dir, hybrid_dir),'OMI_BEHR_*.mat',2,city_lon,city_lat,theta,'windvel',windvel,'windcrit',5,'windop','<','rel_box_corners',[0.5 2 0.5 0.5]);
[no2x_mngt5, no2ld_mngt5] = calc_line_density(fullfile(behr_work_dir, monthly_dir),'OMI_BEHR_*.mat',2,city_lon,city_lat,theta,'windvel',windvel,'windcrit',5,'windop','>=','rel_box_corners',[0.5 2 0.5 0.5]);
[no2x_mnlt5, no2ld_mnlt5] = calc_line_density(fullfile(behr_work_dir, monthly_dir),'OMI_BEHR_*.mat',2,city_lon,city_lat,theta,'windvel',windvel,'windcrit',5,'windop','<','rel_box_corners',[0.5 2 0.5 0.5]);
[no2x_mn288gt5, no2ld_mn288gt5] = calc_line_density(fullfile(behr_work_dir, coarse_mn_dir),'OMI_BEHR_*.mat',2,city_lon,city_lat,theta,'windvel',windvel,'windcrit',5,'windop','>=','rel_box_corners',[0.5 2 0.5 0.5]);
[no2x_mn288lt5, no2ld_mn288lt5] = calc_line_density(fullfile(behr_work_dir, coarse_mn_dir),'OMI_BEHR_*.mat',2,city_lon,city_lat,theta,'windvel',windvel,'windcrit',5,'windop','<','rel_box_corners',[0.5 2 0.5 0.5]);

%% Basic fitting
fprintf('Doing simple fitting of the line densities\n');
fmincon_out_level = 'none'; % could be 'none', 'final', 'iter'
fitfxn = 'ssresid'; % can be 'ssresid' or 'unexvar'
[f_hygt5.ffit, f_hygt5.emgfit, f_hygt5.stats, f_hygt5.f0, ~, f_hygt5.fitresults] = fit_line_density(no2x_hygt5, no2ld_hygt5, fmincon_out_level,'fittype',fitfxn);
[f_hylt5.ffit, f_hylt5.emgfit, f_hylt5.stats, f_hylt5.f0, ~, f_hylt5.fitresults] = fit_line_density(no2x_hylt5, no2ld_hylt5, fmincon_out_level,'fittype',fitfxn);
[f_mngt5.ffit, f_mngt5.emgfit, f_mngt5.stats, f_mngt5.f0, ~, f_mngt5.fitresults] = fit_line_density(no2x_mngt5, no2ld_mngt5, fmincon_out_level,'fittype',fitfxn);
[f_mnlt5.ffit, f_mnlt5.emgfit, f_mnlt5.stats, f_mnlt5.f0, ~, f_mnlt5.fitresults] = fit_line_density(no2x_mnlt5, no2ld_mnlt5, fmincon_out_level,'fittype',fitfxn);
[f_mn288gt5.ffit, f_mn288gt5.emgfit, f_mn288gt5.stats, f_mn288gt5.f0, ~, f_mn288gt5.fitresults] = fit_line_density(no2x_mn288gt5, no2ld_mn288gt5, fmincon_out_level,'fittype',fitfxn);
[f_mn288lt5.ffit, f_mn288lt5.emgfit, f_mn288lt5.stats, f_mn288lt5.f0, ~, f_mn288lt5.fitresults] = fit_line_density(no2x_mn288lt5, no2ld_mn288lt5, fmincon_out_level,'fittype',fitfxn);

%% Variational fitting
fprintf('Sampling the residuals using the standard deviations derived from the Hessian\n');
[Ffix_hygt5, F_hygt5] = fit_line_density_variation(no2x_hygt5, no2ld_hygt5, 30,'fittype',fitfxn);
[Ffix_hylt5, F_hylt5] = fit_line_density_variation(no2x_hylt5, no2ld_hylt5, 30,'fittype',fitfxn);
[Ffix_mngt5, F_mngt5] = fit_line_density_variation(no2x_mngt5, no2ld_mngt5, 30,'fittype',fitfxn);
[Ffix_mnlt5, F_mnlt5] = fit_line_density_variation(no2x_mnlt5, no2ld_mnlt5, 30,'fittype',fitfxn);
[Ffix_mn288gt5, F_mn288gt5] = fit_line_density_variation(no2x_mn288gt5, no2ld_mn288gt5, 30,'fittype',fitfxn);
[Ffix_mn288lt5, F_mn288lt5] = fit_line_density_variation(no2x_mn288lt5, no2ld_mn288lt5, 30,'fittype',fitfxn);

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

if ~exist('mc_mn288gt5','var')
    fprintf('\t Coarse >= 5\n');
    [mc_mn288gt5.x, mc_mn288gt5.f0, mc_mn288gt5.mcacc] = emg_fit_montecarlo(no2x_mn288gt5, no2ld_mn288gt5);
    sd_mc_mn288gt5 = nan(1,5);
    for a=1:5; sd_mc_mn288gt5(a) = std(mc_mn288gt5.x(:,a)); end
end

if ~exist('mc_mn288lt5','var')
    fprintf('\t Coarse < 5\n');
    [mc_mn288lt5.x, mc_mn288lt5.f0, mc_mn288lt5.mcacc] = emg_fit_montecarlo(no2x_mn288lt5, no2ld_mn288lt5);
    sd_mc_mn288lt5 = nan(1,5);
    for a=1:5; sd_mc_mn288lt5(a) = std(mc_mn288lt5.x(:,a)); end
end

%% Variational sampling over the range of values from Monte Carlo
[Ffix_hygt5_mc, F_hygt5_mc] = fit_line_density_variation(no2x_hygt5, no2ld_hygt5, 30, 'sd', sd_mc_hygt5);
[Ffix_hylt5_mc, F_hylt5_mc] = fit_line_density_variation(no2x_hylt5, no2ld_hylt5, 30, 'sd', sd_mc_hylt5);
[Ffix_mngt5_mc, F_mngt5_mc] = fit_line_density_variation(no2x_mngt5, no2ld_mngt5, 30, 'sd', sd_mc_mngt5);
[Ffix_mnlt5_mc, F_mnlt5_mc] = fit_line_density_variation(no2x_mnlt5, no2ld_mnlt5, 30, 'sd', sd_mc_mnlt5);
[Ffix_mn288gt5_mc, F_mn288gt5_mc] = fit_line_density_variation(no2x_mn288gt5, no2ld_mn288gt5, 30, 'sd', sd_mc_mn288gt5);
[Ffix_mn288lt5_mc, F_mn288lt5_mc] = fit_line_density_variation(no2x_mn288lt5, no2ld_mn288lt5, 30, 'sd', sd_mc_mn288lt5);

%% Average correlation among fit params
corr_hygt5 = corrcov(inv(f_hygt5.fitresults.fminunc_hessian));
corr_hylt5 = corrcov(inv(f_hylt5.fitresults.fminunc_hessian));
corr_mngt5 = corrcov(inv(f_mngt5.fitresults.fminunc_hessian));
corr_mnlt5 = corrcov(inv(f_mnlt5.fitresults.fminunc_hessian));
corr_mn288gt5 = corrcov(inv(f_mn288gt5.fitresults.fminunc_hessian));
corr_mn288lt5 = corrcov(inv(f_mn288lt5.fitresults.fminunc_hessian));

avg_corr_mat = nanmean(cat(3,corr_hygt5,corr_hylt5,corr_mngt5,corr_mnlt5,corr_mn288gt5,corr_mn288lt5),3);

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

if ~exist('mc_mn288gt5_fixf0','var')
    fprintf('\t Coarse >= 5\n');
    fixf0 = [f_mn288gt5.ffit.a, f_mn288gt5.ffit.x_0, f_mn288gt5.ffit.mu_x, f_mn288gt5.ffit.sigma_x, f_mn288gt5.ffit.B];
    [mc_mn288gt5_fixf0.x, mc_mn288gt5_fixf0.f0, mc_mn288gt5_fixf0.mcacc] = emg_fit_montecarlo(no2x_mn288gt5, no2ld_mn288gt5, 'f0', fixf0);
    sd_mc_mn288gt5_fixf0 = nan(1,5);
    for a=1:5; sd_mc_mn288gt5_fixf0(a) = std(mc_mn288gt5_fixf0.x(:,a)); end
end

if ~exist('mc_mn288lt5_fixf0','var')
    fprintf('\t Coarse < 5\n');
    fixf0 = [f_mn288lt5.ffit.a, f_mn288lt5.ffit.x_0, f_mn288lt5.ffit.mu_x, f_mn288lt5.ffit.sigma_x, f_mn288lt5.ffit.B];
    [mc_mn288lt5_fixf0.x, mc_mn288lt5_fixf0.f0, mc_mn288lt5_fixf0.mcacc] = emg_fit_montecarlo(no2x_mn288lt5, no2ld_mn288lt5, 'f0', fixf0);
    sd_mc_mn288lt5_fixf0 = nan(1,5);
    for a=1:5; sd_mc_mn288lt5_fixf0(a) = std(mc_mn288lt5_fixf0.x(:,a)); end
end
