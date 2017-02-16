%%compute_edge_pts
%%arr 12/10/2007

function [latedge, lonedge] = compute_edge_pts(lat, lon)

n = length(lat); 

latedge = zeros(n+1,1);
lonedge = zeros(n+1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extrapolate edge pixels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%JLL 14 Feb 2014: Generally, I think NaNs only occur in latitude and
%longitude when the satellite is operating in super-zoom mode and half the
%swath is cut down (I think - could be wrong). This means that half the
%swath has NaNs for lat/lon, so we need to find the beginning and end of
%the valid pixels to do the extrapolations.
%JLL 5-12-2014: Assume that the pixel lat and lon is a linear function of
%"pixel number" and so calculate the edge of the first pixel by assuming
%the distance to the center of the pixel on either side is the same.
xx = find(~isnan(lat) & ~isnan(lon),2,'first');

[x1,y1,z1] = sph2cart(lon(xx(1))*pi/180,lat(xx(1))*pi/180,1); 
[x2,y2,z2] = sph2cart(lon(xx(2))*pi/180,lat(xx(2))*pi/180,1);

% x1, y1, z1 are the center point of the edge pixel, so its corners will be
% 1/2 dx, dy, dz away from it
dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
xi = x1 + 0.5*dx; yi = y1 + 0.5*dy; zi = z1 + 0.5*dz;
[poslon,poslat,~] = cart2sph(xi,yi,zi);
lonedge(xx(1)) = poslon / pi * 180;
latedge(xx(1)) = poslat / pi * 180;

xx2 = find(~isnan(lat) & ~isnan(lon),2,'last');

% x1, y1, z1 here need to the the pixel on the edge, hence why they are
% found as the second element of xx2
[x1,y1,z1] = sph2cart(lon(xx2(2))*pi/180,lat(xx2(2))*pi/180,1);
[x2,y2,z2] = sph2cart(lon(xx2(1))*pi/180,lat(xx2(1))*pi/180,1);
    
dx = x1 - x2; dy = y1 - y2; dz = z1 - z2;
xi = x1 + 0.5*dx; yi = y1 + 0.5*dy; zi = z1 + 0.5*dz;
[poslon,poslat,~] = cart2sph(xi,yi,zi);
lonedge(xx2(2)+1) = poslon / pi * 180;
latedge(xx2(2)+1) = poslat / pi * 180;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolate inner pixels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%JLL 2-18-2014: Assume that the edge between pixels is halfway in between
%the center points. Remember that lonedge and latedge have n+1 elements. 
for i = xx(2):xx2(2);
    [x1,y1,z1] = sph2cart(lon(i-1)*pi/180,lat(i-1)*pi/180,1);
    [x2,y2,z2] = sph2cart(lon(i)*pi/180,lat(i)*pi/180,1);
    xi = 0.5*(x1 + x2); yi = 0.5*(y1 + y2); zi = 0.5*(z1 + z2);
    [poslon,poslat,~] = cart2sph(xi,yi,zi);
    lonedge(i) = poslon / pi * 180;
    latedge(i) = poslat / pi * 180;
end
  

    