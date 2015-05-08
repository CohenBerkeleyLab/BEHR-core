function [ cld_rad_frac ] = cld_rad_frac_interp( cld_geo_frac )
%cld_rad_frac_interp Find the OMI radiative cloud fraction corresponding to a given geometric fraction.
%   OMI AMF calculations rely on the use of a radiative cloud fraction that
%   accounts for the fact that, since clouds are usually more reflective
%   than the ground beneath them, the percentage of light they contribute
%   to the scene is usually greater than the geometric fraction of a scene
%   they actually take up.  Directly calculating these is slightly
%   involved, so in lieu of that calculation, this function will use the
%   relationship given in the presentation linked below to estimate the
%   radiative cloud fraction corresponding to a geometric cloud fraction
%   for SZA = 45, VZA = 40, cloud optical thickness = 10 and surface albedo
%   = 0.05.
%
%   Source:
%   http://discover-aq.larc.nasa.gov/pdf/2010STM/Bhartia_discover_aq_pres.pdf)
%   


% This will add the variables cld_geo_frac_vec and cld_rad_frac_vec to the
% workspace. This file is located in BEHR/Utils.  These vectors are the x &
% y values of the solid blue line on the 4th page of the above ppt (NO2 -
% 420 nm, alb = 0.05).
load('cld_frac_vectors.mat');

cld_rad_frac = interp1(cld_geo_frac_vec,cld_rad_frac_vec,cld_geo_frac);

% If the geometric fraction is too low, the interpolation will return a
% NaN. This is undesirable - a cloud fraction of 0 should be a radiance
% fraction of 0 too. Since a cloud radiance fraction of 0.047 corresponds to
% a geometric cloud fraction of 0.02, we'll estimate the radiance fraction
% as twice the geometric fraction. It's very rough, but should do.

if cld_geo_frac < 0.02 && isnan(cld_rad_frac)
    cld_rad_frac = 2*cld_geo_frac;
end


end

