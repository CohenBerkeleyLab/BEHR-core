#include "fintrf.h"
C user_brdf_kernels - calculates MODIS BRDF kernels k_vol and k_geo
C    [k_vol, k_geo] = user_brdf_kernels(sza, vza, raa) - computes the 
C       volumetric and geometric kernels (k_vol, k_geo) given solar
C       and viewing zenith angles (sza, vza) and the relative azimuth
C       angles (raa) in degrees. The mex interface converts them to
C       radians before calling user_brdf_kernels itself.
C Gateway routine

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)

      implicit none

C mexFunction arguments
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs, mrows, nrows, asize
C Arguments for the brdf test
      double precision sza, vza, raa, kgeo, kvol
      double precision, parameter :: d2r = 3.141592653 / 180.0d0

C Function declarations
      mwPointer mxGetPr
      mwPointer mxCreateDoubleMatrix
      integer mxIsNumeric
      mwPointer mxGetM, mxGetN

C Pointers to input and output arrays
      mwPointer sza_ptr, vza_ptr, raa_ptr, kvol_ptr, kgeo_ptr

C Check for proper number of arguments. 
      if(nrhs .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:user_brdf_kernels:nInput',
     +                           'Three inputs required.')
      elseif(nlhs .gt. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:user_brdf_kernels:nOutput',
     +                           'Too many output arguments.')
      endif

C Check that the inputs are numbers.
      if(mxIsNumeric(prhs(1)) .eq. 0 .or.
     +   mxIsNumeric(prhs(2)) .eq. 0 .or.
     +   mxIsNumeric(prhs(3)) .eq. 0) then
         call mexErrMsgIdAndTxt ('MATLAB:user_brdf_kernels:NonNumeric',
     +                           'All inputs must be numeric.')
      endif

C Get the input pointers
      sza_ptr = mxGetPr(prhs(1))
      vza_ptr = mxGetPr(prhs(2))
      raa_ptr = mxGetPr(prhs(3))

C Get the size of each input
      mrows = mxGetM(prhs(1))
      nrows = mxGetN(prhs(1))
      if ( mrows*nrows .ne. 1 ) then
         call mexErrMsgIdAndTxt ('MATLAB:user_brdf_kernels:NonNumeric',
     +                           'All inputs must be scalar.')
      endif
      asize = mrows*nrows
      
      mrows = mxGetM(prhs(2))
      nrows = mxGetN(prhs(2))
      if ( mrows*nrows .ne. 1 ) then
         call mexErrMsgIdAndTxt ('MATLAB:user_brdf_kernels:NonNumeric',
     +                           'All inputs must be scalar.')
      endif

      mrows = mxGetM(prhs(3))
      nrows = mxGetN(prhs(3))
      if ( mrows*nrows .ne. 1 ) then
         call mexErrMsgIdAndTxt ('MATLAB:user_brdf_kernels:NonNumeric',
     +                           'All inputs must be scalar.')
      endif

C Copy the input pointers
      call mxCopyPtrToReal8(sza_ptr, sza, asize)
      call mxCopyPtrToReal8(vza_ptr, vza, asize)
      call mxCopyPtrToReal8(raa_ptr, raa, asize)

C Prepare the output data
      plhs(1) = mxCreateDoubleMatrix(mrows,nrows,0)
      plhs(2) = mxCreateDoubleMatrix(mrows,nrows,0)
      kvol_ptr = mxGetPr(plhs(1))
      kgeo_ptr = mxGetPr(plhs(2))

C All that just to call the subroutine...
      vza = vza*d2r
      sza = sza*d2r
      raa = raa*d2r
      call user_brdf_kernels(vza, sza, raa, kvol, kgeo)

C Copy the output to the pointers which will be returned to matlab
      call mxCopyReal8ToPtr(kvol,kvol_ptr,asize)
      call mxCopyReal8ToPtr(kgeo,kgeo_ptr,asize)

      end

c #include "user_brdf_kernels.f"
