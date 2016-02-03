function [ alb ] = modis_brdf_alb_poly( f_iso, f_vol, f_geo, sza )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

% Convert to radians
sza = sza ./ 180 .* pi;

% The polynomial coefficients given on p. 16 of the MOD43 TBD.
g0_iso = 1.0;
% g1 and g2 iso are 0

g0_vol = -0.007574;
g1_vol = -0.070987;
g2_vol = 0.307588;

g0_geo = -1.284909;
g1_geo = -0.166314;
g2_geo = 0.041840;

% Using the given coefficients and the polynomial for black-sky albedo,
% calculate the albedo

alb = f_iso .* g0_iso +...
    f_vol .* (g0_vol + g1_vol .* sza.^2 + g2_vol .* sza.^3) +...
    f_geo .* (g0_geo + g1_geo .* sza.^2 + g2_geo .* sza.^3);

end

