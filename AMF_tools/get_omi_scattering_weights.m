function [ sw ] = get_omi_scattering_weights( pressures, sza, vza, raa, alb, surf_pres, lon, lat, mon )
%get_omi_scattering_weights Returns a vector of scattering weights
%   Inputs:
%       pressures = the pressure levels to return the scattering weights
%       for, range allowed [1020, 0.8]
%
%       sza = solar zenith angle, range [0, 90]
%
%       vza = viewing zenith angle, range [0, 90]
%
%       raa = relative azimuth angle, range [0, 180]
%
%       alb = surface albedo, range [0, 1]
%
%       surf_pres = surface_pressure, range [111.552, 1013]
%
%       lon = longitude used for temperature correction, range [-180, 180]
%
%       lat = latitude used for temperature correction, range [-90, 90]
%
%       mon = month used for the temperature correction, range [1 12]
%
%   Josh Laughner <joshlaugh5@gmail.com> 13 May 2015

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(pressures)
    E.badinput('pressures is empty - did you try to specify it as 1013:200 or similar? Put the smaller number first');
elseif ~isnumeric(pressures) || ~isvector(pressures)
    E.badinput('pressures must be a numeric vector')
end

if ~isscalar(sza) || ~isnumeric(sza) || sza < 0 || sza > 90
    E.badinput('sza must be a scalar between 0 and 90');
end
if ~isscalar(vza) || ~isnumeric(vza) || vza < 0 || vza > 90
    E.badinput('vza must be a scalar between 0 and 90');
end
if ~isscalar(raa) || ~isnumeric(raa) || raa < 0 || raa > 180
    E.badinput('raa must be a scalar between 0 and 180');
end
if ~isscalar(alb) || ~isnumeric(alb) || alb < 0 || alb > 1
    E.badinput('alb must be a scalar between 0 and 1');
end
if ~isscalar(surf_pres) || ~isnumeric(surf_pres) || surf_pres < 111.552 || sza > 1013
    E.badinput('surf_pres must be a scalar between 111.552 and 1013');
end
if ~isscalar(lon) || ~isnumeric(lon) || lon < -180 || lon > 180
    E.badinput('lon must be a scalar between -180 and 180');
end
if ~isscalar(lat) || ~isnumeric(lat) || lat < -90 || lat > 90
    E.badinput('lat must be a scalar between -90 and 90');
end
if ~isscalar(mon) || ~isnumeric(mon) || mon < 1 || mon > 12
    E.badinput('mon must be a scalar between 1 and 12');
end

%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FXN %%%%%
%%%%%%%%%%%%%%%%%%%%

[damf_file, tmp_file] = amf_filepaths;

temperature = rNmcTmp2(tmp_file, pressures, lon, lat, mon);

alpha = 1 - 0.003 * (temperature - 220);   % temperature correction factor vector
% Clip between 0.1 and 10 (following omiAmfAK2)
alpha_i=max(alpha,0.1);
alpha = min(alpha_i,10);

damf = rDamf2(damf_file,pressures,sza,vza,raa,alb,surf_pres);

sw = damf .* alpha;

end

