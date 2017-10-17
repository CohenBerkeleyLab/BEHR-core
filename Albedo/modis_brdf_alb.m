function [ alb ] = modis_brdf_alb( f_iso, f_vol, f_geo, sza, vza, raa )
%MODIS_BRDF_ALB Computes the black sky albedo from MODIS MCD43C1 data
%   [ ALB ] = MODIS_BRDF_ALB( F_ISO, F_VOL, F_GEO, SZA, VZA, RAA ) returns
%   the effective surface reflectance/albedo ALB given the three MODIS BRDF
%   parameters (F_ISO, F_VOL, and F_GEO) and the viewing geometry (SZA,
%   VZA, RAA). Internally, calls MODIF_BRDF_KERNELS with the given viewing
%   geometry to calculate the needed kernels then combines them with the
%   weights F_ISO, F_VOL, and F_GEO.
%
%   A BRDF (bidirectional reflectance distribution function) is a
%   mathematical way of representing the directional dependence of a
%   surface's reflectance. The idea is that scattering is likely
%   asymmetrical (to varying degrees for different surfaces) and so viewing
%   a scene in a forward scattering vs. backward scattering mode can have
%   different effective reflectances.
%
%   See "BRDF explained" on the MODIS BRDF product website
%   (https://www.umb.edu/spectralmass/terra_aqua_modis/modis)

[k_vol, k_geo] = modis_brdf_kernels(sza, vza, raa);
alb = f_iso + f_vol .* k_vol + f_geo .* k_geo;

end

