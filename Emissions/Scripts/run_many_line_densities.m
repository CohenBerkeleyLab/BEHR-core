function run_many_line_densities
% I'm not sitting around while these things run!

%cities = {'Atlanta','Birmingham','Montgomery'};
cities = {'Atlanta', 'Birmingham'};
wind_crits = {3,4,5};%{'mean',3,5};
box = [1.0 2.0 0.5 0.5];%[1.0 2.0 0.5 0.5];
save_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/EMG fits/Autorun/FullDaily-NumObs/50km-side-earthrel-no0to-112-lonwt14-1822UTC';
for a=1:numel(cities)
    for b=1:numel(wind_crits)
        if isnumeric(wind_crits{b})
            wind_crit_str = sprintf('%.1f',wind_crits{b});
        else
            wind_crit_str = wind_crits{b};
        end
        fprintf('%s: wind crit = %s\n',cities{a}, wind_crit_str);
        [S, wind_crit] = gen_line_densities( cities{a}, wind_crits{b}, box );
        save_line_densities(S, wind_crit, save_dir, cities{a}, box);
    end
end
end


function [ S, wind_crit ] = gen_line_densities( city, wind_crit, box )

E = JLLErrors;

switch city
    case 'Atlanta'
        wind_file = 'Atlanta-Wind-Conditions-1900UTC-5layers-earthrel.mat';
    case 'Birmingham'
        wind_file = 'Birmingham-Wind-Conditions-1900UTC-5layers-earthrel.mat';
    case 'Montgomery'
        wind_file = 'Montgomery-Wind-Conditions-1900UTC-5layers-earthrel.mat';
    otherwise
        E.badinput('City %s not recognized',city);
end

interp_bool = false;
%box = [1 2 0.5 0.5]; %[0.5 2 0.5 0.5];

homedir = getenv('HOME');
behr_work_dir = fullfile(homedir,'Documents','MATLAB','BEHR','Workspaces','Wind speed');
%hybrid_dir = 'SE US BEHR Hybrid - No ghost';
wrf_dir = 'SE US WRF Hourly';
hybrid_dir = 'SE US BEHR Hourly - No ghost - lw 14.0 overpass - 18-22 UTC';
%monthly_dir = 'SE US BEHR Monthly - No ghost';
monthly_dir = 'SE US BEHR Monthly - No ghost - lw 14.0 overpass - 18-22 UTC';
%coarse_mn_dir = 'SE US BEHR Monthly - No ghost - Coarse WRF';
coarse_mn_dir = 'SE US BEHR Monthly - No ghost - Coarse WRF - lw 14.0 overpass - 18-22 UTC';


% loads theta and windvel plus city lat and lon
load(fullfile(behr_work_dir,wind_file));

if ischar(wind_crit) && strcmpi(wind_crit,'mean')
    wind_crit = mean(windvel);
elseif ischar(wind_crit)
    E.badinput('wind_crit must be a number or the string ''mean''');
end

if strcmpi(city, 'Atlanta')
    gtcrit = windvel >= wind_crit & (theta < -112.5 | theta > 0);
    ltcrit = windvel < wind_crit & (theta < -112.5 | theta > 0);
else
    gtcrit = windvel >= wind_crit;
    ltcrit = windvel < wind_crit;
end

F = dir(fullfile(behr_work_dir, hybrid_dir, 'OMI_BEHR_*.mat'));

fdnums = nan(size(F));
for a=1:numel(fdnums)
    [s,e] = regexp(F(a).name,'\d\d\d\d\d\d\d\d');
    fdnums(a) = datenum(F(a).name(s:e), 'yyyymmdd');
end
F = F(fdnums >= datenum('2013-06-01'));

fprintf('\tFast hybrid\n')
[S.no2x_hyfast, S.no2ld_hyfast, S.no2ldstd_hyfast, S.lon_hyfast, S.lat_hyfast, S.no2cd_hyfast, ~, S.num_obs_hyfast, S.indiv_swaths_hyfast, S.debug_cell_hyfast] = calc_line_density(fullfile(behr_work_dir, hybrid_dir),F,city_lon,city_lat,theta,'crit_logical',gtcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);
fprintf('\tSlow hybrid\n')
[S.no2x_hyslow, S.no2ld_hyslow, S.no2ldstd_hyslow, S.lon_hyslow, S.lat_hyslow, S.no2cd_hyslow, ~, S.num_obs_hyslow, S.indiv_swaths_hyslow, S.debug_cell_hyslow] = calc_line_density(fullfile(behr_work_dir, hybrid_dir),F,city_lon,city_lat,theta,'crit_logical',ltcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);
fprintf('\tFast monthly\n')
[S.no2x_mnfast, S.no2ld_mnfast, S.no2ldstd_mnfast, S.lon_mnfast, S.lat_mnfast, S.no2cd_mnfast, ~, S.num_obs_mnfast, S.indiv_swaths_mnfast, S.debug_cell_mnfast] = calc_line_density(fullfile(behr_work_dir, monthly_dir),F,city_lon,city_lat,theta,'crit_logical',gtcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);
fprintf('\tSlow monthly\n')
[S.no2x_mnslow, S.no2ld_mnslow, S.no2ldstd_mnslow, S.lon_mnslow, S.lat_mnslow, S.no2cd_mnslow, ~, S.num_obs_mnslow, S.indiv_swaths_mnslow, S.debug_cell_mnslow] = calc_line_density(fullfile(behr_work_dir, monthly_dir),F,city_lon,city_lat,theta,'crit_logical',ltcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);
fprintf('\tFast coarse monthly\n')
[S.no2x_mn108fast, S.no2ld_mn108fast, S.no2ldstd_mn108fast, S.lon_mn108fast, S.lat_mn108fast, S.no2cd_mn108fast, ~, S.num_obs_mn108fast, S.indiv_swaths_mn108fast, S.debug_cell_mn108fast] = calc_line_density(fullfile(behr_work_dir, coarse_mn_dir),F,city_lon,city_lat,theta,'crit_logical',gtcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);
fprintf('\tSlow coarse monthly\n')
[S.no2x_mn108slow, S.no2ld_mn108slow, S.no2ldstd_mn108slow, S.lon_mn108slow, S.lat_mn180slow, S.no2cd_mn108slow, ~, S.num_obs_mn108slow, S.indiv_swaths_mn108slow, S.debug_cell_mn108fast] = calc_line_density(fullfile(behr_work_dir, coarse_mn_dir),F,city_lon,city_lat,theta,'crit_logical',ltcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);

return

F = dir(fullfile(behr_work_dir, wrf_dir, '*.mat'));
fdnums = nan(size(F));
for a=1:numel(fdnums)
    [s,e] = regexp(F(a).name,'\d\d\d\d\d\d\d\d');
    fdnums(a) = datenum(F(a).name(s:e), 'yyyymmdd');
end
F = F(fdnums >= datenum('2013-06-01'));
fprintf('\tFast WRF\n')
[S.no2x_wrffast, S.no2ld_wrffast, S.no2ldstd_wrffast, S.lon_wrffast, S.lat_wrffast, S.no2cd_wrffast, S.num_obs_wrffast] = calc_line_density(fullfile(behr_work_dir, wrf_dir),F,city_lon,city_lat,theta,'crit_logical',gtcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);
fprintf('\tSlow WRF\n')
[S.no2x_wrfslow, S.no2ld_wrfslow, S.no2ldstd_wrfslow, S.lon_wrfslow, S.lat_wrfslow, S.no2cd_wrfslow, S.num_obs_wrfslow] = calc_line_density(fullfile(behr_work_dir, wrf_dir),F,city_lon,city_lat,theta,'crit_logical',ltcrit,'rel_box_corners', box, 'interp', interp_bool, 'DEBUG_LEVEL',0);
end

function save_line_densities(S, wind_crit, save_dir, city_name, box, ismean)
if ~exist('ismean','var')
    ismean = false;
end

% Make the fields of S regular variables

fns = fieldnames(S);
for f=1:numel(fns)
    eval(sprintf('%1$s = S.%1$s;',fns{f}));
end

if ismean
    wind_crit_str = strrep(sprintf('Mean_%.1f',wind_crit),'.','pt');
else
    wind_crit_str = strrep(sprintf('%.1f',wind_crit),'.','pt');
end

%save_name = sprintf('NoInterp-LineDensities-%s-%dkmFor-%dkmBack-%dkmSide-AllAngles-%s.mat',wind_crit_str,box(2)*100,box(1)*100,box(3)*100,city_name);
save_name = sprintf('NoInterp-LineDensities-%s-%dkmFor-%dkmBack-%dkmSide-No0to-112-%s.mat',wind_crit_str,box(2)*100,box(1)*100,box(3)*100,city_name);

fprintf('\n\t Saving %s \n\n',fullfile(save_dir,save_name));
save(fullfile(save_dir,save_name),fns{:},'wind_crit');

end
