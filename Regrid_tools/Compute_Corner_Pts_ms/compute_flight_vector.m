%%compute_flight_vector
%%arr 12/10/2007

%JLL 2-18-2014: The flight vector is the straight line between two points
%on the earth's surface. It is the chord connecting those points, not the
%arc across the surface.
function [flightvector] = compute_flight_vector(lat, lon, lat1, lon1)

earthradius = 6378.5;

n = length(lat); 

flightvector = zeros(n,3);

for i = 1:n;
    
    [x1,y1,z1] = sph2cart(lon(i)*pi/180,lat(i)*pi/180,earthradius);
    [x2,y2,z2] = sph2cart(lon1(i)*pi/180,lat1(i)*pi/180,earthradius);
    
    dx = x2-x1; dy = y2-y1; dz = z2-z1;
    
    flightvector(i,1) = dx;
    flightvector(i,2) = dy;
    flightvector(i,3) = dz;
end


