function [ alb ] = modis_brdf_alb_poly( f_iso, f_vol, f_geo, sza )
% MODIS_BRDF_ALB_POLY Calculate SZA-dependent black-sky albedo
%   ALB = MODIS_BRDF_ALB_POLY( F_ISO, F_VOL, F_GEO, SZA ) Given vectors of
%   BRDF coefficients F_ISO, F_VOL, and F_GEO along with the SZAs,
%   calculate the black sky albedo expected for a surface with those
%   coefficients for that SZA.
%
%   More information is given on pp. 15-16 of the MOD43 theoretical basis
%   document (https://modis.gsfc.nasa.gov/data/atbd/atbd_mod09.pdf).
%   Normally, for a full BRDF representation, analytical or semi-empirical
%   kernel expressions are necessary to compute the exact reflectance for a
%   sun-satellite geometry. These expressions are fairly complicated (c.f.
%   modis_brdf_kernels.m), and if one is only interested in the total
%   albedo of a surface (i.e. the total light reflected in all directions)
%   for a given solar zenith angle, it would be necessary to integrate over
%   all possible viewing angles. However, it has been found that the SZA
%   dependence of that integral can be well approximated with a set of
%   polynomials, contained in this function.

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

