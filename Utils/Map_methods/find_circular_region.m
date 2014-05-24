function [ within_circle ] = find_circular_region( center_lon, center_lat, dist_km, lon, lat  )
% Finds all pixels given by the grid matrices with centers "lon" & "lat" within "dist_km" (distance in kilometers) of "center_lon" and "center_lat"
%   Pass a center latitude and longitude and a distance in kilometers as
%   the first three arguments.  The lon and lat matrices can theoretically
%   be gridded or vectors, either will work - however, the pixel indices must 
%   match between the lon & lat matrices and whatever data you will ultimately use this for.   
%   This function returns a logical matrix of pixels that fall within the
%   distance specified of the center pixel.
%
%   Josh Laughner <joshlaugh5@gmail.com> 19 May 2014

within_circle = false(size(lat));
for i = 1:numel(lat)
    if m_lldist([center_lon, lon(i)], [center_lat, lat(i)]) <= dist_km;
        within_circle(i) = 1;
    end
end


end

