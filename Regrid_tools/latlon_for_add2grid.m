function [ Longitude, Latitude ] = latlon_for_add2grid( lonbdy, latbdy, reslon, reslat )
%LATLON_FOR_ADD2GRID Creates the latitude/longitude grids for BEHR gridding
%   [ LONGITUDE, LATITUDE ] = LATLONG_FOR_ADD2GRID( LONBDY, LATBDY, RESLON, RESLAT )
%   Returns the grids LONGITUDE and LATITUDE given the min and max
%   longitude (LONBDY) and latitude (LATBDY) and the desired resolution in
%   degrees for each dimension.

E = JLLErrors;

if numel(lonbdy) ~= 2 || numel(latbdy) ~= 2
    E.badinput('LONBDY and LATBDY must be 2 element vectors')
end
if ~isscalar(reslon) || ~isnumeric(reslon) || ~isscalar(reslat) || ~isnumeric(reslon) || reslon <= 0 || reslat <= 0
    E.badinput('RESLON and RESLAT must be scalar, positive numbers')
end

Latitude=(min(latbdy)+reslat/2):reslat:(max(latbdy)-reslat/2); 
Latitude=Latitude'; 
nlat = numel(Latitude);

Longitude=(min(lonbdy)+reslon/2):reslon:(max(lonbdy)-reslon/2); 
nlon = numel(Longitude);

Latitude=repmat(Latitude,1,nlon);
Longitude=repmat(Longitude,nlat,1);

end

