;written by Pepijn Veefkind

function fWeightFlightDir, distance, x

; distance is the distance to the satellite

; assume that the FWHM of IFOV is 1 degree:
fwhm_deg = 1.

; the distance on ground becomes:
fwhm_km = 2 * tan( 0.5 * fwhm_deg /180. * !PI ) * distance

; using flat-topped gaussian with exponent 4:
; f(x) = exp( -c1 * ( x - x0)^4 )
; compute the constant c1, knowing that f(0.5*fwhm_km) = 0.5
c1 = alog( 0.5 )/ (0.5 * fwhm_km)^4

; x0 is the center of the pixel, i.e. the subsatellite point when half
; of the exposure time has passed.
; x0 is defined 0, thus all distances are wrt to the center of the pixel.

; x0 is the middle of the IFOV. At t=0 X0 is defined to 0. At time t
; x0 = kt, where k is the ground speed:
k = 2 * !PI * 6378.5 / ( 100. * 60. )


; the OI master clock period is 2 seconds. The integration is done from
; t=-1 to t=+1 seconds.
mcp = 2.


; the weight for a specific point x on earth at time t becomes
; w(x,t) = exp( c1 * ( x -x0(t))^4),
; w(x,t) = exp( c1 * ( x - k t )^4 )

; To compute the weight for position x for a master clock period, we have
; to integrate from t=-1 to +1. This is done numerically.
nsteps = 5
dt = mcp / nsteps

wx = 0.0
for i=0, nsteps-1 do begin
    t = (i + 0.5) * dt - 0.5 * mcp

    expo = double( c1 * ( x - k * t )^4 )

    if ( expo lt -700d) then wxt = 0.0D else  wxt = exp( expo )

    wx = wx + wxt

endfor

wx = wx / double( nsteps )

return, wx

end





function fweightdistance, distance, X = x, WX=wx

; computes the array of weights on the groud,
; as a function of the distance to the center of the ground pixel

nsteps = 1001
xstart = -50.
xend   =  50.


dx = (xend - xstart) / float( nsteps-1)

wx = dblarr( nsteps)
x  = dblarr( nsteps)

for i=0, nsteps-1 do begin

    x[i] = xstart + i * dx
    wx[i] = fWeightFlightDir(distance, x[i] )

endfor

; nofmalize wx
wx = wx / max( wx )

; determine the FWHM
res = min( abs(wx-0.5), index )
fwhm = abs( 2*x[index] )

return, fwhm

end




function fcomputedistance, lat, lon, satLat, satLon, satHeight


earthRadius = 6378.5

earthPosSpherical = [ lon, lat, 1.]
earthPos = CV_COORD( /DEGREES, /DOUBLE, FROM_SPHERE=earthPosSpherical, /TO_RECT)


satPosSpherical = [ satLon, satLat, satHeight/earthRadius + 1.]
satPos = CV_COORD( /DEGREES, /DOUBLE, FROM_SPHERE=satPosSpherical, /TO_RECT)


distance = norm( satPos - earthPos)

distance = distance * earthRadius

return, distance

end



function fcomputeEdgePoints, lat, lon, latEdge, lonEdge


n = n_elements( lat )

latEdge = fltarr( n +1  )
lonEdge = fltarr( n +1  )

; first do the interpolation part
for i=1, n-1 do begin
    pos0 = CV_COORD( /DEGREES, /DOUBLE, FROM_SPHERE=[lon[i-1], lat[i-1], 1.], /TO_RECT)
    pos1 = CV_COORD( /DEGREES, /DOUBLE, FROM_SPHERE=[lon[i], lat[i], 1.], /TO_RECT)

    posi = 0.5* ( pos0 + pos1 )

    posLatLon = CV_COORD( /DEGREES, /DOUBLE, FROM_RECT=posi, /TO_SPHERE)

    latEdge[i] = posLatLon[1]
    lonEdge[i] = posLatLon[0]
endfor


; do the extrapolations
pos0 = CV_COORD( /DEGREES, /DOUBLE, FROM_SPHERE=[lon[0], lat[0], 1.], /TO_RECT)
pos1 = CV_COORD( /DEGREES, /DOUBLE, FROM_SPHERE=[lon[1], lat[1], 1.], /TO_RECT)

difPos = pos0 - pos1
posi = pos0 + difpos

posLatLon = CV_COORD( /DEGREES, /DOUBLE, FROM_RECT=posi, /TO_SPHERE)

latEdge[0] = posLatLon[1]
lonEdge[0] = posLatLon[0]

; do the extrapolations
pos0 = CV_COORD( /DEGREES, /DOUBLE, FROM_SPHERE=[lon[n-1], lat[n-1], 1.], /TO_RECT)
pos1 = CV_COORD( /DEGREES, /DOUBLE, FROM_SPHERE=[lon[n-2], lat[n-2], 1.], /TO_RECT)

difPos = pos0 - pos1
posi = pos0 + difpos

posLatLon = CV_COORD( /DEGREES, /DOUBLE, FROM_RECT=posi, /TO_SPHERE)

latEdge[n] = posLatLon[1]
lonEdge[n] = posLatLon[0]

return, 0

end



function computeFlightVector, lat, lon, lat1, lon1

earthRadius = 6378.5

n = n_elements( lat )

flightVector = dblarr( n, 3)

for i=0, n-1 do begin

    pos0 = CV_COORD( /DEGREES, /DOUBLE, FROM_SPHERE=[lon[i], lat[i], earthRadius], /TO_RECT)
    pos1 = CV_COORD( /DEGREES, /DOUBLE, FROM_SPHERE=[lon1[i], lat1[i], earthRadius], /TO_RECT)

    flightVector[i, *] = (pos1 - pos0)

endfor


return, flightVector

end



function computeCornerPoints, lat, lon, flightvector, fwhm, latCorner, lonCorner

earthRadius = 6378.5

n = n_elements( lat )

latCorner = dblarr(n, 2)
lonCorner = dblarr(n, 2)

for i=0, n-1 do begin

    pos0 = CV_COORD( /DEGREES, /DOUBLE, FROM_SPHERE=[lon[i], lat[i], earthRadius], /TO_RECT)

    FV = reform( flightVector[i, *], 3 )

    scale = ( 0.5* fwhm[i] ) / norm( FV )

    pos1 = pos0 - scale * FV
    pos2 = pos0 + scale * FV

    pos1LatLon = CV_COORD( /DEGREES, /DOUBLE, FROM_RECT=pos1, /TO_SPHERE)
    pos2LatLon = CV_COORD( /DEGREES, /DOUBLE, FROM_RECT=pos2, /TO_SPHERE)

    latCorner[i,0] = pos1LatLon[1]
    latCorner[i,1] = pos2LatLon[1]
    lonCorner[i,0] = pos1LatLon[0]
    lonCorner[i,1] = pos2LatLon[0]

endfor


return, 0

end


function fcomputeFwhm, lat, lon, satLat, satLon, satHeight, DISTANCE=distance


n = n_elements( lat )

fwhm = dblarr(n)
distance = dblarr(n)

for i=0, n-1 do begin

    distance[i] = fcomputedistance( lat[i], lon[i], satLat, satLon, satHeight)

    fwhm[i] = fweightdistance( distance[i])

    ;if (i eq 0) then print, 'NrPixel    Latitude     Longitude     Distance    FWHM'
    ;print, i, lat[i], lon[i], distance[i], fwhm[i]

endfor

return, fwhm

end

pro TestPepijnsPixelsize
; prepare graphics
;simplewin

; compile the routines
;@omipixelcorners
	;set_plot, 'win'
	;Window, 1, Title='Pixelsize', XSize=800,YSize=600
	;Device, Decomposed=0
	set_plot, 'ps'
	DEVICE, FILE='Pepijn Pixelsnew.ps', /COLOR, BITS=8, /portrait


	; extract some real omi data:
	lFileName = "C:\Data\OMI\OMINO2\Lv2\2005\04\15\OMI-Aura_L2-OMNO2_2005m0415t0119-o03987_v002-2005m1221t084732.he5"
	ReadLv2cFiles, lFileName, 'ColumnAmountNO2Trop', latitude, longitude, Data, cc, amf,sza, DataError, satLatf, satLonf, satHeightf
	;data = h5_parse( filename, /READ)

	itimes = 1000
color = 1
	for i = 0, 1 do begin



	    itimes = itimes + 1

	    lat = latitude[*,itimes]
	    lon = longitude[*,itimes]

	    lat1 = latitude[*,itimes+1]
	    lon1 = longitude[*,itimes+1]

	    satLat = satLatf[itimes]
	    satLon = satLonf[itimes]
	    satHeight = satHeightf[itimes]

	    satHeight = 1e-3 * satHeight  ; convert to km


	    ; compute the ground pixel edge lat lon in the across track direction
	    res = fcomputeEdgePoints( lat, lon, latEdge, lonEdge )
	    res = fcomputeEdgePoints( lat1, lon1, latEdge1, lonEdge1 )
			help, latEdge
	    ; compute flight vector
	    flightVector = computeFlightVector(latEdge, lonEdge, latEdge1, lonEdge1)
			help, latEdge1

	    ; compute the FWHM in the flight direction if not already calculated
	    if ( n_elements(fwhm) ne n_elements(latEdge) ) then $
				fwhm = fcomputeFwhm( latEdge, lonEdge, satLat, satLon, satHeight, DISTANCE=distance)$
			ELSE PRINT, 'already calculated'
	    ; compute corner points
	    res = computeCornerPoints( latEdge, lonEdge, flightvector, fwhm, latCorner, lonCorner)
			help, latCorner

	    ; do some visualization
	    ;if i eq 0 then plot, lon, lat, xtitle='Longitude', ytitle='Latitude', title=filename, yrange=[-50,-40], xrange=[-20,20]
	    ;if i eq 0 then plot, lon, lat, xtitle='Longitude', ytitle='Latitude', XTICKLEN=1.0,YTICKLEN=1.0, XGRIDSTYLE=2, yGRIDSTYLE=2

	    color = i + 1

	    ;oplot, lon, lat, psym=1, color = color
	    ;if i eq 0 then plot, loncorner[*,0], latcorner[*,0], xtitle='Longitude', ytitle='Latitude', XTICKLEN=1.0,YTICKLEN=1.0, XGRIDSTYLE=2, yGRIDSTYLE=2,yrange=[-3,2]$
	    if i eq 0 then plot, loncorner[*,0], latcorner[*,0], xtitle='Longitude', ytitle='Latitude', XTICKLEN=1.0,YTICKLEN=1.0, XGRIDSTYLE=2, yGRIDSTYLE=2$
	    else oplot, loncorner[*,0], latcorner[*,0]
	    oplot, loncorner[*,1], latcorner[*,1]
	    for ip=0,60 do oplot,  loncorner[ip,*], latcorner[ip,*]
    	;oplot, longitude[29,*], latitude[29,*]
	endfor
	device, /close

end


function CornerCoordinates, latitude, longitude, satLatf, satLonf, satHeightf
	;latitude=array(60,y), y=number of scanlines in the orbit, 60 scenes per scanline
	;longitude=array(60,y)
	;satLatf=array(60,y)
	;satLonf=array(60,y)
	;satHeightf=array(60,y)
	;returns array(60,y,2,5),2=lon,lat,5=corners(5th=center)
	dims = SIZE(latitude, /DIMENSIONS)
	help, latitude
	print, dims
	;stop
	corners = fltarr(dims(0),dims(1),2,5)

	for y = 0, dims(1)-2 do begin
		lat = latitude[*,y]
		lon = longitude[*,y]

		lat1 = latitude[*,y+1]
		lon1 = longitude[*,y+1]

		satLat = satLatf[y]
		satLon = satLonf[y]
		satHeight = satHeightf[y]

		satHeight = 1e-3 * satHeight  ; convert to km
		; compute the ground pixel edge lat lon in the across track direction
		res = fcomputeEdgePoints( lat, lon, latEdge, lonEdge )
		res = fcomputeEdgePoints( lat1, lon1, latEdge1, lonEdge1 )
		; compute flight vector
		flightVector = computeFlightVector(latEdge, lonEdge, latEdge1, lonEdge1)
		; compute the FWHM in the flight direction if not already calculated
		if ( n_elements(fwhm) ne n_elements(latEdge) ) then $
			fwhm = fcomputeFwhm( latEdge, lonEdge, satLat, satLon, satHeight, DISTANCE=distance);$
		;ELSE PRINT, 'fwhm already calculated'
		; compute corner points
		res = computeCornerPoints( latEdge, lonEdge, flightvector, fwhm, latCorner, lonCorner)
		for x=0,59 do begin
			corners(x,y,0,0)=loncorner[x,0]
			corners(x,y,1,0)=latcorner[x,0]
			corners(x,y,0,1)=loncorner[x+1,0]
			corners(x,y,1,1)=latcorner[x+1,0]
			corners(x,y,0,2)=loncorner[x+1,1]
			corners(x,y,1,2)=latcorner[x+1,1]
			corners(x,y,0,3)=loncorner[x,1]
			corners(x,y,1,3)=latcorner[x,1]
			corners(x,y,0,4)=lon(x)
			corners(x,y,1,4)=lat(x)

			if (abs(loncorner[x+1,0]-loncorner[x,1]) gt 180.0) $
			or (abs(loncorner[x+1,1]-loncorner[x,0]) gt 180.0) $
			or (abs(loncorner[x+1,0]-loncorner[x,0]) gt 180.0) $
			or (abs(loncorner[x+1,1]-loncorner[x,1]) gt 180.0) then begin
				if corners(x,y,0,0) lt 0.0 then corners(x,y,0,0)+=360.0
				if corners(x,y,0,1) lt 0.0 then corners(x,y,0,1)+=360.0
				if corners(x,y,0,2) lt 0.0 then corners(x,y,0,2)+=360.0
				if corners(x,y,0,3) lt 0.0 then corners(x,y,0,3)+=360.0
			end

		end
	end
	return, corners
end