%%weight_distance
%%arr 12/10/2007

%computes the array of weights on the ground as a function of the distance
%to the center of the ground pixel

function [fwhm] = weight_distance(distance) 
%JLL 2-19-2014: These were commented out when I received the file.
%nsteps = 1001;
%start = -50;
%xend = 50;

%dx = (xend - xstart) / (nsteps - 1);
%wx = zeros(nsteps,1);
%x = zeros(nsteps,1);
 
x=linspace(-50,50,101);


wx = weight_flight_dir(distance, x);

%normalize wx
wx = wx / max(wx);

%determine fwhm            
aa = (abs(wx - 0.5));
bb = find(aa==min(aa));

fwhm = abs(2 * x(bb(1)));
