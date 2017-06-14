      SUBROUTINE USER_BRDF_KERNELS(zen_refl, zen_inc, phi_refl, k_vol, 
     $ k_geo)
c     ================
c     Input Parameters
c     ================
      double precision zen_refl, zen_inc, phi_refl

c     ================
c     Output Parameters 
c     ================
      double precision k_vol, k_geo

c     ==================
c     Internal Variables
c     ==================
      double precision, parameter :: Pie = 3.1415927d0
      double precision, parameter :: b_r = 1.0d0
      double precision, parameter :: h_b = 2.0d0
      double precision xi, xi_prime, zen_refl_prime, zen_inc_prime
      double precision D2, cos_t, t, O
      integer, parameter :: dbunit = 10
      logical, parameter :: db_bool = .false.

      if (db_bool) open(dbunit,file='bkdb.txt')

c cos(xi) essentially describes the total angle between the incoming direct
c and outgoing viewing paths. From Roujean et al. (1992).
      xi = dacos(dcos(zen_inc) * dcos(zen_refl) + dsin(zen_inc) * 
     $ dsin(zen_refl) * dcos(phi_refl))

      k_vol = ((Pie/2 - xi) * dcos(xi) + dsin(xi))/(dcos(zen_inc) +
     $ dcos(zen_refl)) - Pie/4.0d0

    
c k_geo is defined as:
c   k_geo = O(sza, vza, raa) - sec(sza') - sec(vza') + 1/2(1 + cos(xi'))*sec(sza')*sec(vza')
c so first we need the quantities O, sza', vza', and xi'.
c
c The primed version of the angles account for the shape of the objects in
c the model - i.e., it will be different than the basic angle if the
c objects are assumed to be ellipsoidal rather than circular/spherical (as
c we do here).
c
c The idea is described in Wanner. Basically, the idea is to consider an
c ellipsoidal object illuminated from the initial angle. The question is,
c what angle would a spherical object need to be illuminated from to
c produce the same shadow as the ellipsoid at the first angle? That new
c angle is the primed version.
c
c The ratio b/r is the ratio of the two axes of the ellipsoid, if one, the
c objects are spherical. h/b is how high off the ground the objects are,
c relative to their size. MODIS assumes this to be 2.

      zen_inc_prime = datan(b_r * dtan(zen_inc))
      zen_refl_prime = datan(b_r * dtan(zen_refl))
      xi_prime = dacos(dcos(zen_inc_prime)*dcos(zen_refl_prime) + 
     $ dsin(zen_inc_prime) * dsin(zen_refl_prime) * dcos(phi_refl))

      if (db_bool) write(dbunit, *) 'xi_prime=', xi_prime

c The following three quantities deal with the overlap of the objects in
c the model.

      if (db_bool) write(dbunit, *) 'zen_inc_prime=', zen_inc_prime
      if (db_bool) write(dbunit, *) 'zen_refl_prime=', zen_refl_prime
      if (db_bool) write(dbunit, *) 'phi_refl=', phi_refl

      if (db_bool) write(dbunit, *) 'D2 first term=',
     $ dtan(zen_inc_prime)**2.0d0 + dtan(zen_refl_prime)**2.0d0
      if (db_bool) write(dbunit, *) 'D2 second term=',
     $ 2.0d0*dtan(zen_inc_prime) * dtan(zen_refl_prime) * dcos(phi_refl)

      D2 = dtan(zen_inc_prime)**2.0d0 + dtan(zen_refl_prime)**2.0d0 - 
     $ 2.0d0*dtan(zen_inc_prime) * dtan(zen_refl_prime) * dcos(phi_refl)

      cos_t = h_b * sqrt( (D2 + (dtan(zen_inc_prime) * 
     $ dtan(zen_refl_prime) * dsin(phi_refl))**2.0d0) ) / 
     $ (1.0d0/dcos(zen_inc_prime)+1.0d0/dcos(zen_refl_prime))

      if (db_bool) write(dbunit, *) 'D2=', D2
      if (db_bool) write(dbunit, *) 'cos_t=', cos_t
      
c The TBD indicates on p. 14 that cos(t) should be constrained to the range
c [-1, 1] as "values outside this range imply no overlap [between solar and
c viewing shadows] and should be disregarded."  I believe the line of
c reasoning goes something like this:
c   If you read through the Wanner paper, you'll find that for a sparse
c leaf canopy (which is what this is trying to model) the goal of the K_geo
c is to account for the effect the canopy (a.k.a. "crowns") as a whole has
c on the reflectivity based on how much of the ground is obscured (while
c K_vol handles the scattering effect of the individual leaves).  In there,
c they point out that there is reflectivity from both the ground and the
c canopies, so they account for both, assuming in this case, that there are
c few enough individual canopies that they do not shadow each other, only
c the ground.
c   Further, they go on to assume that the shadows are completely black, so
c that any ground shadowed either from canopies blocking the incoming
c sunlight or the observer's view of the ground does not contribute to the
c reflectance.  This area is calculated by considering multiplying the
c average area of the hypothetical canopies by the secants of the viewing
c and solar zenith angles to account for the geometric transformation of
c the shadow of a spherical canopy to an ellipsoid shadow on the ground
c (consider shining a flashlight on a baseball from above versus the side -
c from the side, the shadow is elongated).
c   However the last thing that must be accounted for is the overlap of
c shadows of incoming light and from the observer's POV. The term O handles
c that, and t should vary from 0 to pi, i.e. from no overlap to complete
c overlap.
      cos_t = min(max(cos_t,0.0d0),1.0d0)
      t = dacos(cos_t)
      O = 1.0d0/Pie * (t - dsin(t) * dcos(t)) * 
     $ (1.0d0/cos(zen_inc_prime) + 1.0d0/cos(zen_refl_prime))

      if (db_bool) write(dbunit, *) 'cos_t limited=', cos_t
      if (db_bool) write(dbunit, *) 't=', t
      if (db_bool) write(dbunit, *) 'O=', O

      k_geo = O - 1.0d0/dcos(zen_inc_prime) - 1.0d0/cos(zen_refl_prime) 
     $ + 0.5d0 * (1.0d0 + dcos(xi_prime)) * 1.0d0/dcos(zen_inc_prime)
     $ * 1.0d0/cos(zen_refl_prime)

      if (db_bool) close(dbunit)

      return
      END
