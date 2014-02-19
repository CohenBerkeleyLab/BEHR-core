%%compute_distance
%%arr 12/10/2007

function [distance] = compute_distance(lat, lon, satlat, satlon, satalt)

earthradius = 6378.5;

[earthposx,earthposy,earthposz] = sph2cart(lon.*pi./180, lat.*pi./180, 1);
earthpos = [earthposx', earthposy', earthposz'];

[satposx,satposy,satposz] = sph2cart(satlon.*pi./180, satlat.*pi./180, satalt./earthradius + 1);
satpos = [satposx, satposy, satposz];
satpos=repmat(satpos,length(lat),1);

distance = norm(satpos - earthpos);
distance = distance .* earthradius;

