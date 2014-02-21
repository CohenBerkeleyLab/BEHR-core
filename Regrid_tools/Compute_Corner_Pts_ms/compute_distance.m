%%compute_distance
%%arr 12/10/2007

%JLl 2-19-2014: I believe that lat and lon SHOULD be single numbers but am
%not positive.  This then returns the distance from the satellite to the
%latitude and longitude input.
function [distance] = compute_distance(lat, lon, satlat, satlon, satalt)

earthradius = 6378.5;

[earthposx,earthposy,earthposz] = sph2cart(lon.*pi./180, lat.*pi./180, 1);
earthpos = [earthposx', earthposy', earthposz'];

[satposx,satposy,satposz] = sph2cart(satlon.*pi./180, satlat.*pi./180, satalt./earthradius + 1); %JLL 2-19-2014: satalt is alt. above the surface, the +1 converts that into the proper r-coordinate
satpos = [satposx, satposy, satposz];
satpos=repmat(satpos,length(lat),1);

distance = norm(satpos - earthpos); %JLL 2-19-2014: Assuming that satpos and earthpos are simply 3D vectors, this finds the distance between earth and the satellite.
distance = distance .* earthradius; %JLL 2-19-2014: Since the conversion to cartesian coordinates assumed a unit sphere, this converts back to km.

