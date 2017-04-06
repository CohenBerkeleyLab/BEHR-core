%%compute_edge_pts
%%arr 12/10/2007

function [latedge, lonedge] = compute_edge_pts(lat, lon)

n = length(lat); 

latedge = zeros(n+1,1);
lonedge = zeros(n+1,1);

%the extrapolations
%JLL 5-12-2014: Assume that the pixel lat and lon is a linear function of
%"pixel number" and so calculate the edge of the first pixel by assuming
%the distance to the center of the pixel on either side is the same.
xx = find(~isnan(lat) & ~isnan(lon),2,'first');

[x1,y1,z1] = sph2cart(lon(xx(1))*pi/180,lat(xx(1))*pi/180,1); 
[x2,y2,z2] = sph2cart(lon(xx(2))*pi/180,lat(xx(2))*pi/180,1);
    
dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
xi = x1 + dx; yi = y1 + dy; zi = z1 + dz;
[poslon,poslat,Z1] = cart2sph(xi,yi,zi);
lonedge(1) = poslon / pi * 180;
latedge(1) = poslat / pi * 180;


[x1,y1,z1] = sph2cart(lon(n)*pi/180,lat(n)*pi/180,1);
[x2,y2,z2] = sph2cart(lon(n-1)*pi/180,lat(n-1)*pi/180,1);
    
dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
xi = x1 + dx; yi = y1 + dy; zi = z1 + dz;
[poslon,poslat,Z1] = cart2sph(xi,yi,zi);
lonedge(n+1) = poslon / pi * 180;
latedge(n+1) = poslat / pi * 180;
    
    
%the interpolations
%JLL 2-18-2014: Assume that the edge between pixels is halfway in between
%the center points. 
for i = 2:n;
    [x1,y1,z1] = sph2cart(lon(i-1)*pi/180,lat(i-1)*pi/180,1);
    [x2,y2,z2] = sph2cart(lon(i)*pi/180,lat(i)*pi/180,1);
    xi = 0.5*(x1 + x2); yi = 0.5*(y1 + y2); zi = 0.5*(z1 + z2);
    [poslon,poslat,Z1] = cart2sph(xi,yi,zi);
    lonedge(i) = poslon / pi * 180;
    latedge(i) = poslat / pi * 180;
end
  

    