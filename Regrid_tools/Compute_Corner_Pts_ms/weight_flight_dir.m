%%weight_flight_dir
%%arr 12/10/2007

%distance is the distance to the satellite

function [wx] = weight_flight_dir(distance, x)

%assume that the fwhm of ifov is 1 deg:
fwhm_deg = 1;

%the distance on the ground becomes:
fwhm_km = 2 .* tan(0.5 .* fwhm_deg ./ 180 .* pi) .* distance;

%using flat-topped gaussian with exponent 4:
%f(x) = exp(-c1 * (x - x0)^4)
%compute the constant c1, knowing that f(0.5 * fwhm_km) = 0.5
c1 = log(0.5) ./ (0.5 .* fwhm_km).^4;


%x0 is the center of the pixel, ie the subsatellite point when half of the
%exposure time has passed.
%x0 is defined 0; thus all distances are wrt the center of the pixel.

%x0 is the middle of the ifov. at t=0, x0 is defined to 0. at time t, x0 =
%kt, where k is the ground speed:
k = 2 .* pi .* 6378.5 ./ (100 .* 60);

%the OI master clock period is 2 seconds. the integration is done from t=-1
%to t=+1 seconds.
mcp = 2;

%the weight for a specific point x on the earth at time t becomes:
%w(x,t) = exp(c1 * (x1 - x0(t))^4);
%w(x,t) = exp(c1 * (x - kt)^4);

%to compute the wieght for position x for a master clock period, we have to
%integrate from t=-1 to t=+1.  this is done numerically:
nsteps = 5;
dt = mcp ./ nsteps;

wx = 0;
for i = 1:nsteps;
    t = (i - 0.5) .* dt - 0.5 .* mcp;
    expo = (c1 .* (x - k .* t).^4);
    if (expo < -700);
        wxt = 0;
    else
       wxt = exp(expo);
    end
    wx = wx + wxt;
end
wx = wx ./ nsteps;
       