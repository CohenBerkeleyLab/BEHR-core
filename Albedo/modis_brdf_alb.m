function [ alb ] = modis_brdf_alb( f_iso, f_vol, f_geo, sza, vza, raa )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[k_vol, k_geo] = modis_brdf_kernels(sza, vza, raa);
alb = f_iso + f_vol .* k_vol + f_geo .* k_geo;

end

