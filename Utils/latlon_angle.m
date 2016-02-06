function [ angle ] = latlon_angle( lon1, lat1, lon2, lat2, from_east )
%latlon_angle Calculates angle between two lat/lon coordinates
%   Calculates the angle between two lat/lon coordinates (in degrees). By
%   default, it will return an angle in degrees clockwise from north (which
%   is apparently standard for navigation). However, in some cases it
%   makes more sense to do it in degrees counter-clockwise from each
%   (particularly if comparing to wind vector quantities, since this
%   definition follows the usual trigonometric definition of angle w.r.t. x
%   & y vectors). To use this definition, pass true as the optional fifth
%   argument.

E = JLLErrors;
n = ndims(lon1); 
sz = size(lon1);
if ndims(lat1) ~= n || any(size(lat1) ~= sz)
    E.sizeMismatch('lon1','lat1')
elseif ndims(lon2) ~= n || any(size(lon2) ~= sz)
    E.sizeMismatch('lon1','lon2')
elseif ndims(lat2) ~= n || any(size(lat2) ~= sz)
    E.sizeMismatch('lon1','lat2')
end

angle = nan(sz);

if ~exist('from_east','var')
    from_east = false;
elseif (~islogical(from_east) && ~isnumeric(from_east)) || ~isscalar(from_east)
    E.badinput('from_east, if given, must be a scalar logical or numeric value');
end

x = cosd(lat1).*sind(lat2)-sind(lat1).*cosd(lat2).*cosd(lon2-lon1);
y = sind(lon2-lon1).*cosd(lat2);

if from_east
    angle = mod(atan2d( x , y ), 360);
else
    angle = mod(atan2d( y , x ), 360);
end

end

