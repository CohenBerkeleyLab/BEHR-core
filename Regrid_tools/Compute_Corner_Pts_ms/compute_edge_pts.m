%%compute_edge_pts
%%arr 12/10/2007

function [latedge, lonedge] = compute_edge_pts(lat, lon)

n = length(lat); 

latedge = zeros(n+1,1);
lonedge = zeros(n+1,1);

%the extrapolations

[x1,y1,z1] = sph2cart(lon(1)*pi/180,lat(1)*pi/180,1);
[x2,y2,z2] = sph2cart(lon(2)*pi/180,lat(2)*pi/180,1);
    
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

for i = 2:n;
    [x1,y1,z1] = sph2cart(lon(i-1)*pi/180,lat(i-1)*pi/180,1);
    [x2,y2,z2] = sph2cart(lon(i)*pi/180,lat(i)*pi/180,1);
    xi = 0.5*(x1 + x2); yi = 0.5*(y1 + y2); zi = 0.5*(z1 + z2);
    [poslon,poslat,Z1] = cart2sph(xi,yi,zi);
    lonedge(i) = poslon / pi * 180;
    latedge(i) = poslat / pi * 180;
end
  

    