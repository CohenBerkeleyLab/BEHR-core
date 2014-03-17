%%compute_corner_pts
%%arr 12/10/2007

function [latcorner, loncorner] = compute_corner_pts(lat, lon, flightvector, fwhm)

earthradius = 6378.5; %km

n = length(lat);

latcorner = zeros(n,2);
loncorner = zeros(n,2);

for i = 1:n;
    [x1,y1,z1] = sph2cart(lon(i)*pi/180,lat(i)*pi/180,earthradius); %JLL 17 Mar 2014: This time we need the earth's radius because we are interested in the actual physical distance between the points, not in units of earth_radius
    FV = flightvector(i,:); %JLL 17 Mar 2014: The 3D vector between adjacent rows
    scale = (0.5 * fwhm(i)) / norm(FV); %JLL 17 Mar 2014: This will give the ratio between the "ideal" pixel size defined by its FWHM and its actual size, determined by the flight vector
    
    posx1 = x1 - scale * FV(1); posy1 = y1 - scale * FV(2); posz1 = z1 - scale * FV(3); %JLL 17 Mar 2014: Find the corners for this side of the pixel
    posx2 = x1 + scale * FV(1); posy2 = y1 + scale * FV(2); posz2 = z1 + scale * FV(3);
    
    [poslatlonx1, poslatlony1, poslatlonz1] = cart2sph(posx1,posy1,posz1); %JLL 17 Mar 2014: Return to spherical coordinates
    [poslatlonx2, poslatlony2, poslatlonz2] = cart2sph(posx2,posy2,posz2);
    
    latcorner(i,1) = poslatlony1 / pi * 180; %JLL 17 Mar 2014: Convert back to degrees
    latcorner(i,2) = poslatlony2 / pi * 180;
    loncorner(i,1) = poslatlonx1 / pi * 180;
    loncorner(i,2) = poslatlonx2 / pi * 180;
end
