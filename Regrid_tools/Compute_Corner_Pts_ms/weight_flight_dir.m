%%weight_flight_dir
%%arr 12/10/2007

%distance is the distance to the satellite

function [wx] = weight_flight_dir(distance, x)

%assume that the fwhm of ifov is 1 deg: (JLL 2-19-2014: ifov =
%intrinsic field of view, c.f. http://www.knmi.nl/omi/research/instrument/characteristics.php)
%This is the ifov along track
fwhm_deg = 1;

%the distance on the ground becomes:
fwhm_km = 2 .* tan(0.5 .* fwhm_deg ./ 180 .* pi) .* distance;

%using flat-topped gaussian with exponent 4:
%f(x) = exp(-c1 * (x - x0)^4)
%compute the constant c1, knowing that f(0.5 * fwhm_km) = 0.5
c1 = log(0.5) ./ (0.5 .* fwhm_km).^4;
%JLL 2-19-2014: This is the definition of ifov, that ifov is the field of
%view that has >0.5 pixel response.  This constant fixes the gaussian s.t.
%at 1/2 fwhm_km from the center, the reponse is 0.5.


%x0 is the center of the pixel, ie the subsatellite point when half of the
%exposure time has passed.
%x0 is defined 0; thus all distances are wrt the center of the pixel.

%x0 is the middle of the ifov. at t=0, x0 is defined to 0. at time t, x0 =
%kt, where k is the ground speed:
%JLL 2-19-2014: k = 2*pi*r / (orbit period in s).  The OMI orbit period is
%~ 100 min (98.83 by Wikipedia).
k = 2 .* pi .* 6378.5 ./ (100 .* 60);

%the OI master clock period is 2 seconds. the integration is done from t=-1
%to t=+1 seconds.
mcp = 2;

%the weight for a specific point x on the earth at time t becomes:
%w(x,t) = exp(c1 * (x1 - x0(t))^4);
%w(x,t) = exp(c1 * (x - kt)^4);

%to compute the weight for position x for a master clock period, we have to
%integrate from t=-1 to t=+1.  this is done numerically:
nsteps = 5;
dt = mcp ./ nsteps;

wx = 0;
for i = 1:nsteps;
    t = (i - 0.5) .* dt - 0.5 .* mcp; %JLL 17 Mar 2014: Set t to be at the middle of each time step
    expo = (c1 .* (x - k .* t).^4); %JLL 17 Mar 2014: Calculate the exponent;
    if (expo < -700); %JLL 17 Mar 2014: If all the elements of the exponent are < -700, go ahead and just set the result to 0
        wxt = 0;
    else
       wxt = exp(expo); %JLL 17 Mar 2014: Otherwise, go ahead and calculate the expression
    end
    wx = wx + wxt;
end
wx = wx ./ nsteps;
%JLL 2-19-2014: Remember that x is a vector (1x101 matrix) so wx is
%actually a vector as well
       