      program kernel_test
      implicit none

      double precision, dimension(9), parameter ::
     $ sza = (/ 0., 10., 20., 30., 40., 50., 60., 70., 80. /)
      double precision, dimension(6), parameter ::
     $ vza = (/ 0., 15., 30., 45., 60., 75. /)
      double precision, dimension(5), parameter ::
     $ raa = (/ 0., 45., 90., 135., 180. /)

      integer i, j, k
      double precision kvol, kgeo

      do i=1,size(sza)
        do j=1,size(vza)
          do k=1,size(raa)
            call user_brdf_kernels(vza(j), sza(i), raa(k), kvol, kgeo)
            print *, 'SZA=', sza(i), ' VZA=', vza(j), ' RAA=', raa(k),
     $               ' => k_vol=', kvol, ' k_geo=', kgeo
          enddo
        enddo
      enddo

      end

c #include "user_brdf_kernels.f"
