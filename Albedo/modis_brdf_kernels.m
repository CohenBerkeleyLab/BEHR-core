function [ k_vol, k_geo ] = modis_brdf_kernels( sza, vza, raa )
%MODIS_BRDF_KERNELS Calculates k_vol and k_geo for a given geometry
%   The bidirection reflectance distribution function describes the surface
%   reflectivity in a way that should account for variations in reflection
%   depending on viewing geometry. The implementation in the MOD43
%   algorithm uses kernels to account for volumetric and geometric affects
%   on reflectance (i.e. leaves distributed throughout the volume above the
%   surface, or larger obstructions on the surface itself). As far as I can
%   tell, the kernels used for MOD43 do not account for effects like ocean
%   glint.
%
%   This function takes an input viewing geometry (solar zenith, viewing
%   zenith, and relative azimuth angle IN DEGREES) and outputs the kernels
%   necessary to calculate the surface reflectance from the MODIS BRDF
%   algorithm. This would need to be combined with MCD43C1 coefficients to
%   get the actual reflectance.
%
%   IMPORTANT: I believe that these kernel equations expect the traditional
%   description of relative azimuth angle, i.e. |saa - vaa| constrained to
%   [0,180] because on p. 20457 of Roujean et al., it says "...f1 presents
%   a local maximum in the BACKSCATTERING DIRECTION for VZA = SZA, RAA =
%   0..." Therefore, I believe it is important to not use the relative
%   azimuth angle required by the TOMRAD table, which has an extra factor
%   of 180 deg due to the definition of RAA used in that table.
%
%   Further reading: MOD43 Theoretical Basis Document esp. pp. 13-15
%
%   Wanner et al. J. Geophys. Res. (1995), 100, pp. 21077-21089 esp. sec.
%   4.2 for K_geo
%
%   Roujean et al. J. Geophys. Res. (1992), 97, pp. 20455-20468
%
%   Josh Laughner <joshlaugh5@gmail.com> 8 Sept 2015

DEBUG_LEVEL = 0;

% cos(xi) essentially describes the total angle between the incoming direct
% and outgoing viewing paths. From Roujean et al. (1992). Xi needs to be in
% radians for the k_vol eqn.
xi = acos(cosd(sza) .* cosd(vza) + sind(sza) .* sind(vza) .* cosd(raa));

% Calculate k_vol
k_vol = ((pi/2 - xi) .* cos(xi) + sin(xi) ) ./ ( cosd(sza) + cosd(vza) ) - pi/4;

% k_geo is defined as:
%   k_geo = O(sza, vza, raa) - sec(sza') - sec(vza') + 1/2(1 + cos(xi'))*sec(sza')*sec(vza')
% so first we need the quantities O, sza', vza', and xi'.
%
% The primed version of the angles account for the shape of the objects in
% the model - i.e., it will be different than the basic angle if the
% objects are assumed to be ellipsoidal rather than circular/spherical (as
% we do here).
%
% The idea is described in Wanner. Basically, the idea is to consider an
% ellipsoidal object illuminated from the initial angle. The question is,
% what angle would a spherical object need to be illuminated from to
% produce the same shadow as the ellipsoid at the first angle? That new
% angle is the primed version.
%
% The ratio b/r is the ratio of the two axes of the ellipsoid, if one, the
% objects are spherical. h/b is how high off the ground the objects are,
% relative to their size. MODIS assumes this to be 2.
b_r = 1;
h_b = 2;

sza_prime = atan(b_r .* tand(sza));
vza_prime = atan(b_r .* tand(vza));
xi_prime = acos(cos(sza_prime) .* cos(vza_prime) + sin(sza_prime) .* sin(vza_prime) .* cosd(raa));

% The following three quantities deal with the overlap of the objects in
% the model.

if DEBUG_LEVEL > 1
    fprintf('sza_prime = %f\n', sza_prime)
    fprintf('vza_prime = %f\n', vza_prime)
    fprintf('raa = %f\n', raa)
    fprintf('D2 first term = %f\n', tan(sza_prime).^2 + tan(vza_prime).^2)
    fprintf('D2 second term = %f\n', 2 .* tan(sza_prime) .* tan(vza_prime) .* cosd(raa))
end

D2 = tan(sza_prime).^2 + tan(vza_prime).^2 - 2 .* tan(sza_prime) .* tan(vza_prime) .* cosd(raa);
cos_t = h_b .* sqrt( (D2 + (tan(sza_prime) .* tan(vza_prime) .* sind(raa)).^2) ) ./ (sec(sza_prime) + sec(vza_prime));

% The TBD indicates on p. 14 that cos(t) should be constrained to the range
% [-1, 1] as "values outside this range imply no overlap [between solar and
% viewing shadows] and should be disregarded."  I believe the line of
% reasoning goes something like this:
%   If you read through the Wanner paper, you'll find that for a sparse
% leaf canopy (which is what this is trying to model) the goal of the K_geo
% is to account for the effect the canopy (a.k.a. "crowns") as a whole has
% on the reflectivity based on how much of the ground is obscured (while
% K_vol handles the scattering effect of the individual leaves).  In there,
% they point out that there is reflectivity from both the ground and the
% canopies, so they account for both, assuming in this case, that there are
% few enough individual canopies that they do not shadow each other, only
% the ground.
%   Further, they go on to assume that the shadows are completely black, so
% that any ground shadowed either from canopies blocking the incoming
% sunlight or the observer's view of the ground does not contribute to the
% reflectance.  This area is calculated by considering multiplying the
% average area of the hypothetical canopies by the secants of the viewing
% and solar zenith angles to account for the geometric transformation of
% the shadow of a spherical canopy to an ellipsoid shadow on the ground
% (consider shining a flashlight on a baseball from above versus the side -
% from the side, the shadow is elongated).
%   However the last thing that must be accounted for is the overlap of
% shadows of incoming light and from the observer's POV. The term O handles
% that, and t should vary from 0 to pi, i.e. from no overlap to complete
% overlap.

% according to the MOD43 TBD (p 14), this should be contrained between -1
% and 1. Previously I'd found 0 and 1 somewhere... If someone happens to
% find [0,1] again in the future, know that I'd seen [0, 1] too - JLL 22
% Sept 2017
cos_t = min(max(cos_t,-1),1); 
t = acos(cos_t);
O = 1/pi .* (t - sin(t) .* cos(t)) .* (sec(sza_prime) + sec(vza_prime));

% Finally we can compute k_geo
k_geo = O - sec(sza_prime) - sec(vza_prime) + 1/2 .* (1 + cos(xi_prime)) .* sec(sza_prime) .* sec(vza_prime);

if DEBUG_LEVEL > 1
    fprintf('xi_prime = %f\n', xi_prime);
    fprintf('D2 = %f\n', D2);
    fprintf('cos_t = %f\n', cos_t);
    fprintf('t = %f\n', t);
    fprintf('O = %f\n', O);
end

end

