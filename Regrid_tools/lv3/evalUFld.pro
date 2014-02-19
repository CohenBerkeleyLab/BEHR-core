;written by Eric Bucsela

FUNCTION evalUFld,    longitude,  latitude,  unpolFldCoefs

;; Reconstructs the smooth NO2 field (columnAmountNO2Unpolluted) (e15/cm^2)
;; using the coefficients stored in the Level-2 file.
;; Adapted from OMNO2 FORTRAN routines:  genBasisFunctions  and  evalUnpolFld.
;; EJB 2006-08-31
;;
;; Arguments
;; longitude     = array of longitudes (nPxl x nTime) for OMI pixels (deg)
;; latitude      = array of latitudes  (nPxl x nTime) for OMI pixels (deg)
;; unpolFldCoefs = array of coefficients (nLat x nUnpolFldCoefs) (e15/cm^2)


; Get dimensions and initialize variables and arrays .....................................

; Integers:

dum            = SIZE(longitude)     & nPxl=dum(1) &  nTime=dum(2)          ;60  & nTime
dum            = SIZE(unpolFldCoefs) & nLat=dum(1) &  nUnpolFldCoefs=dum(2) ;180 & 5
nLon           = 360
i              =   0
j              =   0


; Real Float (all geographic coordinates are in degrees):

rd             = !pi/180.
longitudeRad   =  0.0
lonStart       = -180.
lonEnd         =  180.
lonRange       =  0.0
latStart       =  -90.
latEnd         =   90.
latRange       =  0.0
aLonC          = FLTARR( nPxl, nTime )
aLon           = FLTARR( nPxl, nTime )
aLonF          = FLTARR( nPxl, nTime )
aLatC          = FLTARR( nPxl, nTime )
aLat           = FLTARR( nPxl, nTime )
aLatF          = FLTARR( nPxl, nTime )
lon0           = FLTARR( nPxl, nTime )

basisFunctions = FLTARR( nLon, nUnpolFldCoefs )

columnAmountNo2Unpolluted = FLTARR( nPxl, nTime )



; Generate basis functions ...............................................................

basisFunctions(*,0)  =  1.0                       ; first set of functions is a constant

FOR j   = 1,  nUnpolFldCoefs-1,  2   DO BEGIN     ; the rest are cosines and sines
  FOR i = 0, nLon-1 DO BEGIN
    longitudeRad = ( lonStart + (lonEnd-lonStart)*(i+0.5) / nLon ) * rd
    basisFunctions( i, j   )  =  COS( longitudeRad * j/2 )
    basisFunctions( i, j+1 )  =  SIN( longitudeRad * j/2 )
  ENDFOR
ENDFOR



; Compute the longitude and latitude indices .............................................

latRange = latEnd - latStart
lonRange = lonEnd - lonStart
lon0     = longitude  - lonRange * (FLOOR( (longitude  - lonStart) /  lonRange )  -  1)
dum      = WHERE(lon0 GE lonEnd) &  IF MIN(dum) GE 0 THEN  lon0(dum) = lon0(dum) - lonRange
aLon     = nLon   *  (lon0      - lonStart) /  lonRange  - 0.5
alat     = nLat   *  (latitude  - latStart) /  latRange  - 0.5



; Compute the unpolluted field by interpolation using fractional grid indices
;  and then summing the product of the basis functions and the coefficients...............

dum   = WHERE (aLon LT 0)         &  IF MIN(dum) GE 0 THEN aLon(dum) = aLon(dum) + nLon
dum   = WHERE (aLat LT 0)         &  IF MIN(dum) GE 0 THEN aLat(dum) =  0
aLonC = CEIL(aLon)
dum   = WHERE (aLonC GE nLon)     &  IF MIN(dum) GE 0 THEN aLonC(dum)=  0
aLonF = FLOOR(aLon)
aLatC = CEIL(aLat)
dum   = WHERE (aLatC GE nLat)     &  IF MIN(dum) GE 0 THEN aLatC(dum) = nLat - 1
aLatF = FLOOR(aLat)


FOR i = 0, nUnpolFldCoefs-1 DO BEGIN
   columnAmountNo2Unpolluted = columnAmountNo2Unpolluted                           $
 +   (    UnpolFldCoefs (aLatF,i)                                                  $
       + (aLat - aLatF) *(UnpolFldCoefs (aLatC,i)  -  UnpolFldCoefs (aLatF,i)) )   $
   * (    basisFunctions(aLonF,i)                                                  $
       + (aLon - aLonF) *(basisFunctions(aLonC,i)  -  basisFunctions(aLonF,i)) )
ENDFOR

RETURN, columnAmountNO2Unpolluted
END

