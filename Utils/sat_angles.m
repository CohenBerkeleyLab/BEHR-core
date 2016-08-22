function [ pix_sza, pix_vza, pix_raa, pix_saa, pix_vaa ] = sat_angles( pix_lon, pix_lat, pix_surfpres, sat_lon, sat_lat, sat_alt, timestr )
%SAT_ANGLES Calculates SZA, VZA, and RAA for satellite pixels.
%   Radiative transfer calculations depend on the sun and satellite
%   position relative to the pixel. This function will calculate those
%   angles. Requires seven inputs:
%
%       pix_lon - a matrix of pixel longitudes.
%
%       pix_lat - a matrix of pixel latitudes. Must be same size as
%       pix_lon.
%
%       pix_surfpres - a matrix of pixel surface pressures in hPa. Must be
%       the same size as pix_lon.
%
%       sat_lat, sat_lon, sat_alt - the coordinates of the satellite
%       (altitude in kilometers). Must either be the same size as the three
%       pixel matrices or be scalars in which case they will be used for
%       all pixels.
%
%       timestr - a string or datenum identifying the date and time (MATLAB
%       must be able to recognize the string format, so yyyy-mm-dd HH:MM:SS
%       is probably best) where the time is given as UTC time OR a
%       structure acceptable for input into sun_position.m. This will not
%       be replicated, so it must correspond to all the pixels input.
%
%   Josh Laughner <joshlaugh5@gmail.com> 21 Jan 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

sz = size(pix_lon);
nd = ndims(pix_lon);

if ndims(pix_lat) ~= nd || any(size(pix_lat) ~= sz) || ndims(pix_surfpres) ~= nd || any(size(pix_surfpres) ~= sz)
    E.badinput('pix_lon, pix_lat, and pix_surfpress must all be the same size');
end
if isscalar(sat_lon)
    sat_lon = repmat(sat_lon,sz);
elseif ndims(sat_lon) ~= nd || any(size(sat_lon) ~= sz)
    E.badinput('sat_lat must be a scalar or the same size as the pixel matrices')
end
if isscalar(sat_lat)
    sat_lat = repmat(sat_lat,sz);
elseif ndims(sat_lat) ~= nd || any(size(sat_lat) ~= sz)
    E.badinput('sat_lat must be a scalar or the same size as the pixel matrices')
end
if isscalar(sat_alt)
    sat_alt = repmat(sat_alt,sz);
elseif ndims(sat_alt) ~= nd || any(size(sat_alt) ~= sz)
    E.badinput('sat_alt must be a scalar or the same size as the pixel matrices')
end
if isstruct(timestr)
    fns = fieldnames(timestr);
    req_fields = {'year','month','day','hour','min','sec','UTC'};
    if any(~ismember(req_fields, fns))
        E.badinput('If giving timestr as a structure, it must contain all the following fields: %s',strjoin(req_fields,', '));
    end
else
    if ischar(timestr)
        try
            timestr = datenum(timestr);
        catch err
            if strcmp(err.identifier,'MATLAB:datestr:ConvertToDateNumber')
                E.badinput('Could not identify date string format')
            else
                rethrow(err);
            end
        end
    end
    time.year = year(timestr);
    time.month = month(timestr);
    time.day = day(timestr);
    time.hour = hour(timestr);
    time.min = minute(timestr);
    time.sec = second(timestr);
    time.UTC = 0;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

r_E = 6375; % earth's radium in km
pix_sza = nan(sz);
pix_vza = nan(sz);
pix_saa = nan(sz); % azimuth angles, will be needed for calculating the rel. azimuth angle.
pix_vaa = nan(sz);

ellip_E = wgs84Ellipsoid('km');

for a=1:numel(pix_lon)
    location.longitude = pix_lon(a);
    location.latitude = pix_lat(a);
    location.altitude = -log(pix_surfpres(a)./1013)*7400;
    sa = sun_position(time,location);
    pix_sza(a) = sa.zenith;
    pix_saa(a) = sa.azimuth;
    
    % The viewing azimuth angle will be easy, as Matlab has built a in
    % function for calculating azimuth angles between points and can
    % generate the WGS84 ellipsoid. The zenith angle will require a bit of
    % vector math: if s_vec is the vector from the earth's center to the
    % satellite and p_vec from earth's center to the pixel, then d_vec will
    % be the difference between them and the zenith angle should be the
    % angle between p_vec and d_vec.
    
    % Zenith first
    s.r = r_E + sat_alt(a);
    s.theta = deg2rad(sat_lon(a));
    s.phi = deg2rad(sat_lat(a));
    
    p.r = r_E + location.altitude/1000;
    p.theta = deg2rad(pix_lon(a));
    p.phi = deg2rad(pix_lat(a));
    
    [s_vec(1), s_vec(2), s_vec(3)] = sph2cart(s.theta, s.phi, s.r);
    [p_vec(1), p_vec(2), p_vec(3)] = sph2cart(p.theta, p.phi, p.r);
    d_vec = s_vec - p_vec;
    
    pix_vza(a) = acosd(dot(p_vec, d_vec)/(norm(p_vec) * norm(d_vec)));
    
    % Azimuth is much easier
    pix_vaa(a) = azimuth(pix_lat(a), pix_lon(a), sat_lat(a), sat_lon(a), ellip_E, 'degrees');
    
end

pix_saa = constrain_angles(pix_saa);
pix_vaa = constrain_angles(pix_vaa);

pix_raa = abs(pix_saa + 180 - pix_vaa);
pix_raa(pix_raa > 180) = 360 - pix_raa(pix_raa > 180);

end

function theta = constrain_angles(theta)
theta(theta>180) = theta(theta>180) - 360;
end
