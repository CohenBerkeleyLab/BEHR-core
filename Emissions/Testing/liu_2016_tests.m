% Tests for the slow-fast EMG method from Liu et al. 2016

mydir = fileparts(mfilename('fullpath'));
currdir = pwd;
cd(fullfile(mydir,'..'));
addpath(genpath(pwd));
cd(currdir);

harbin_lon = 126.5;
harbin_lat = 45.8;

wind_dir = fullfile(behr_repo_dir, 'Workspaces', 'Wind Files', 'ECMWF');
F = dir(fullfile(wind_dir, '*.nc'));
wind_file = cell(size(F));
for a=1:numel(F)
    wind_file{a} = fullfile(wind_dir, F(a).name);
end

[windvel, theta, wind_dnums] = read_ecmwf(wind_file, harbin_lon, harbin_lat);
calc_line_density_sectors(wind_dnums, harbin_lon, harbin_lat, theta, windvel);