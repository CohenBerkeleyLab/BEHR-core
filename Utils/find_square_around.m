function [xx, yy] = find_square_around(lon, lat, center_lon, center_lat, radius)
% FIND_SQUARE_AROUND Returns indicies in a square around a point.
%   [XX, YY] = FIND_SQUARE_AROUND( LON, LAT, CENTER_LON, CENTER_LAT, RADIUS)
%   Given lon and lat as grids, this will find indicies for a square of
%   points centered on center_lon and center_lat. The square will have
%   sides of length 2*radius + 1. Slicing lon/lat as lon(xx,yy), lat(xx,yy)
%   will give the points.

E = JLLErrors;

% Check input
if ~ismatrix(lon)
    E.badinput('lon must be a 2-D array')
elseif ~ismatrix(lat)
    E.badinput('lat must be a 2-D array')
elseif ~all(size(lon) == size(lat))
    E.badinput('lon and lat must be the same size')
elseif ~isscalar(center_lat) || ~isnumeric(center_lat)
    E.badinput('center_lat must be a numeric scalar')
elseif ~isscalar(center_lon) || ~isnumeric(center_lon)
    E.badinput('center_lon must be a numeric scalar')
elseif ~isscalar(radius) || ~isnumeric(radius) || radius < 0 || mod(radius,1) ~= 0
    E.badinput('radius must be a positive scalar integer')
end

del = abs(lon - center_lon) + abs(lat - center_lat);
[~,I] = min(del(:));
[x,y] = ind2sub(size(del),I);
xx = (x-radius):(x+radius);
yy = (y-radius):(y+radius);
end