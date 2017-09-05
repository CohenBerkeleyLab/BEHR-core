function [ Data ] = debug_pseudo_amf_changes( start_date, end_date, resolution )
% DATA = DEBUG_PSEUDO_AMF_CHANGES( START_DATE, END_DATE )
%   This is a one-off function intended to see how the AMF from the average
%   of the hybrid profiles in the pseudo-retrieval for the wind effects
%   paper compares with AMFs from the monthly average. This will hopefully
%   help explain why there's a consistent domain wide increase in the AMF
%   and will see if the AMF from the average hybrid profile is different
%   from the average of all the individual AMFs from the hybrid profiles.
%
%   DEBUG_PSEUDO_AMF_CHANGES( START_DATE, END_DATE, 'coarse' ) will use
%   those derived from the coarse simulations instead.
%
%   Josh Laughner <joshlaugh5@gmail.com> 7 Apr 2016

% First concatenate all the pseudo-retrieval files within the time period.
% Need geographic data, SZA, VZA, RAA, Albedo, terrain pressure, to
% recalculate the AMF plus the average profiles. The former should not vary
% across pseudo retrieval days, so in the end we'll get those from the
% existing file.

if exist('resolution','var') && strcmpi(resolution,'coarse')
    hybrid_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hybrid - No clouds - No ghost - Coarse WRF';
else
    hybrid_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hybrid - No clouds - No ghost';
end

[profs, pres] = cat_sat_data(hybrid_path, {'BEHRNO2apriori','BEHRPressureLevels'},...
    'prefix','OMI','newdim',true,'startdate',start_date,'enddate',end_date,'DEBUG_LEVEL',0);

% Load the first file in the time period
D = load(fullfile(hybrid_path, sprintf('OMI_BEHR_%s.mat',datestr(start_date,'yyyymmdd'))));
Data = D.Data;
clear('D');

% Remove the pressure levels not in the default OMNO2 pressure vector
pressure = [1020 1015 1010 1005 1000 990 980 970 960 945 925 900 875 850 825 800 770 740 700 660 610 560 500 450 400 350 280 200];% 100 50 20 5];
sz = size(profs);
profs_cut = nan([numel(pressure), sz(2:end)]);
for a=1:prod(sz(2:end))
    xx = ismember(pres(:,a),pressure);
    if sum(xx) ~= numel(pressure)
        E.callError('bad_indexing','xx did not find every pressure level needed');
    end
    profs_cut(:,a) = profs(xx,a);
end

profs_cut = nanmean(profs_cut,4);

amf_tools_path = BEHR_paths('amf_tools_dir');
fileTmp = fullfile(amf_tools_path,'nmcTmpYr.txt');
fileDamf = fullfile(amf_tools_path,'damf.txt');

dAmfClr = rDamf2(fileDamf, pressure, Data.SolarZenithAngle, Data.ViewingZenithAngle, Data.RelativeAzimuthAngle, Data.MODISAlbedo, Data.GLOBETerpres);
cloudalbedo = 0.8 * ones(size(Data.Longitude));
cldPres = Data.CloudPressure;
dAmfCld = rDamf2(fileDamf, pressure, Data.SolarZenithAngle, Data.ViewingZenithAngle, Data.RelativeAzimuthAngle, cloudalbedo, Data.CloudPressure);

sz = size(Data.Longitude);
% dAmfClr = nan([size(pres,1),sz]);
% dAmfCld = nan([size(pres,1),sz]);
% Get the clear and cloudy scattering weights
% for a=1:numel(Data.Longitude)
%     fprintf('%d of %d\n',a,numel(Data.Longitude));
%     dAmfClr(:,a) = rDamf2(fileDamf, pres(:,a), Data.SolarZenithAngle(a), Data.ViewingZenithAngle(a), Data.RelativeAzimuthAngle(a), Data.MODISAlbedo(a), Data.GLOBETerpres(a)); %JLL 18 Mar 2014: Interpolate the values in dAmf to the albedo and other conditions input
%     cloudalbedo = 0.8; %Assume that any cloud has an albedo of 0.8
%     cldPres = Data.CloudPressure(:); % We're setting cloud fraction to zero anyway, so this won't matter
%     dAmfCld(:,a) = rDamf2(fileDamf, pres(:,a), Data.SolarZenithAngle(a), Data.ViewingZenithAngle(a), Data.RelativeAzimuthAngle(a), cloudalbedo, Data.CloudPressure(a));
% end


% Get the temperature for the cross section correction
if month(start_date) == month(end_date)
    mon = month(start_date);
else
    warning('The start and end of the run are not in the same month and the average will be used for the temperature correction.')
    mon = month(mean([datenum(start_date), datenum(end_date)]));
end
mon = mon*ones(size(Data.Longitude));
% either temperature also needs to go in a for loop or I need to knock out
% the interpolated pressure levels.
temperature = rNmcTmp2(fileTmp, pressure, Data.Longitude, Data.Latitude, mon);

noGhost = 0;
ak = 1;

[amf, ~, ~, scattering_weights, avg_kernels, no2_prof_interp, sw_plevels, ghost_fraction] = omiAmfAK2(Data.GLOBETerpres, cldPres, Data.CloudFraction, Data.CloudRadianceFraction, pressure, dAmfClr, dAmfCld, temperature, profs_cut, profs_cut, noGhost, ak);

Data.BEHRAMFTrop = amf; %JLL 18 Mar 2014: Save the resulting AMF of the pixel
Data.BEHRGhostFraction = ghost_fraction;
Data.BEHRScatteringWeights = scattering_weights;
Data.BEHRAvgKernels = avg_kernels;
Data.BEHRNO2apriori = no2_prof_interp;
Data.BEHRPressureLevels = sw_plevels;
end

