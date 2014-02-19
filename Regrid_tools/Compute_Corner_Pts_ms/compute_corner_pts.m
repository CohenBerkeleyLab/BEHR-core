%%compute_corner_pts
%%arr 12/10/2007

function [latcorner, loncorner] = compute_corner_pts(lat, lon, flightvector, fwhm)

earthradius = 6378.5;

n = length(lat);

latcorner = zeros(n,2);
loncorner = zeros(n,2);

for i = 1:n;
    [x1,y1,z1] = sph2cart(lon(i)*pi/180,lat(i)*pi/180,earthradius);
    FV = flightvector(i,:);
    scale = (0.5 * fwhm(i)) / norm(FV);
    
    posx1 = x1 - scale * FV(1); posy1 = y1 - scale * FV(2); posz1 = z1 - scale * FV(3);
    posx2 = x1 + scale * FV(1); posy2 = y1 + scale * FV(2); posz2 = z1 + scale * FV(3);
    
    [poslatlonx1, poslatlony1, poslatlonz1] = cart2sph(posx1,posy1,posz1);
    [poslatlonx2, poslatlony2, poslatlonz2] = cart2sph(posx2,posy2,posz2);
    
    latcorner(i,1) = poslatlony1 / pi * 180;
    latcorner(i,2) = poslatlony2 / pi * 180;
    loncorner(i,1) = poslatlonx1 / pi * 180;
    loncorner(i,2) = poslatlonx2 / pi * 180;
end
