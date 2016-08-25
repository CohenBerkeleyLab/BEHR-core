function [ x_box, y_box, lonlim, latlim ] = convert_rel_box_corners( rel_box_corners, center_lon, center_lat )
%CONVERT_REL_BOX_CORNERS Interprets rel_box_corners for rotate_plume
%   [ X_BOX, Y_BOX, LONLIM, LATLIM ] = CONVERT_REL_BOX_CORNERS( REL_BOX_CORNERS, CENTER_LON, CENTER_LAT )
%       This is a utility function for ROTATE_PLUME. ROTATE_PLUME allows
%       the user to define the box to contain the plume as a 4 element
%       vector, REL_BOX_CORNERS, that defines the box as:
%
%       [distance_left, distance_right, distance_up, distance_down]
%
%       in degrees. This function converts this to X_BOX and Y_BOX which
%       give the x- and y- coordinates of the box corners relative to the
%       center and LONLIM and LATLIM which give the longitude and latitude
%       limits to be given to the gridding code. CENTER_LON and CENTER_LAT
%       define the point that REL_BOX_CORNERS should be relative to; they
%       are only needed to return LONLIM and LATLIM.
%
%       This was split into a separate function to allow it to be
%       referenced by calc_line_density to compute the longitude and
%       latitude arrays outside of the parfor loop.


E = JLLErrors;
if ~isvector(rel_box_corners) || numel(rel_box_corners) ~= 4 || ~isnumeric(rel_box_corners) || any(rel_box_corners <= 0) 
    E.badinput('REL_BOX_CORNERS must be a 4-element numeric vector of all positive numbers')
end
if ~isnumeric(center_lon) || ~isscalar(center_lon) || ~isnumeric(center_lat) || ~isscalar(center_lat)
    E.badinput('CENTER_LON and CENTER_LAT must be scalar numbers')
end

x_box = [-rel_box_corners(1), rel_box_corners(2), rel_box_corners(2), -rel_box_corners(1)];
y_box = [-rel_box_corners(3), -rel_box_corners(3), rel_box_corners(4), rel_box_corners(4)];

lonlim = nan(1,2);
latlim = nan(1,2);
lonlim(1) = center_lon + x_box(1);
lonlim(2) = center_lon + x_box(3);
latlim(1) = center_lat + y_box(1);
latlim(2) = center_lat + y_box(3);

end

