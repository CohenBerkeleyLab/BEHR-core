pro Add2Grid, Grid, GridWeight, Data, DataError, CloudCover, AMF, SZA, Coord, lat1, lon1, lat2, lon2, reslat, reslon
;written by Mark Wenig
	Dimensions = SIZE( Data , /DIMENSIONS)
	avweight=0.0
	avweightc=0.0
	seed = 1001L
	for x = 0,Dimensions(0)-1,1 do begin
	for y = 0,Dimensions(1)-1,1 do begin
			lCoord=fltarr(2,5)
			for c = 0,4 do begin
				lCoord(1,c)=(Coord(x,y,1,c)-lat1)/reslat
				lCoord(0,c)=(Coord(x,y,0,c)-lon1)/reslon
			end
			border=10.0
			dummy=where(lCoord(*,*) lt -border/reslat, count1)
			dummy=where(lCoord(1,*) ge ((lat2-lat1+border)/reslat), count2)
			dummy=where(lCoord(0,*) ge ((lon2-lon1+border)/reslon), count3)
			if (count1 eq 0) and (count2 eq 0) and (count3 eq 0) and (Data(x,y) gt 0.0) then begin
				Error=1.0e-16*(18.5+abs(double(x)-29.7)^3.5*0.00028)*DataError(x,y);/Data(x,y);
				Weight=1.0/(Error*Error)
					avweight = avweight+Weight
					avweightc = avweightc + 1.0
				if SZA(x,y) ge 85 then Weight=0.0
				if y gt 10 then if (Coord(29,y,1,4)-Coord(29,y-10,1,4)) lt 0.0 then Weight=0.0;no descending
				if (Weight gt 0.0) and (DataError(x,y) gt 0.0) then begin
					Quadrangle, Grid, lCoord, Data(x,y)*Weight
					Quadrangle, GridWeight, lCoord, Weight
				end
			end
		end
	end

	;print, 'avweight=', avweight/avweightc
end

pro InitGrid, DataGrid, DataGridWeight, maxx, maxy
	DataGrid = fltarr(maxx,maxy)
	DataGridWeight = fltarr(maxx,maxy)
	DataGridError = fltarr(maxx,maxy)
	DataGridCount = intarr(maxx,maxy)
end

pro NormGrid, DataGrid, DataGridWeight, maxx, maxy
 for x = 0,maxx-1 do begin
		for y = 0,maxy-1 do begin
			if DataGridWeight(x,y) gt 0 then begin
				DataGrid(x,y)/=DataGridWeight(x,y)
			end else DataGrid(x,y)=0
		end
	end
end

pro ShowGrid, DataGrid, FileName, maxx, maxy, scale, latmin, latmax, lonmin, lonmax, DatasetName, minvalue, maxvalue, winnumber
	set_plot, 'ps'
	DEVICE, File=FileName, /COLOR, BITS=8, /landscape, xoff=0.0, yoff=0.0,/helvetica, font_size=12
	LOADCT, 33, NCOLORS=255, BOTTOM=1
	Dimensions = SIZE( DataGrid , /DIMENSIONS)
	print, 'Dimensions=',Dimensions
	print, (lonmin+180)*Dimensions(0)/360,(lonmax+180)*Dimensions(0)/360,(latmin+90)*Dimensions(1)/180,(latmax+90)*Dimensions(1)/180
	image = BYTSCL( DataGrid[(lonmin+180)*Dimensions(0)/360:((lonmax+180)*Dimensions(0)/360-1),(latmin+90)*Dimensions(1)/180:((latmax+90)*Dimensions(1)/180-1)] , MAX=maxvalue , MIN=minvalue )
	MAP_SET, 0.0, (lonmin+lonmax)/2.0, /CYLIN, /ISOTROPIC, LIMIT=[latmin, lonmin, latmax, lonmax]
	result = MAP_IMAGE(image, Startx, Starty, COMPRESS=1, LATMIN=latmin, LONMIN=lonmin, LATMAX=latmax, LONMAX=lonmax)
	LOADCT, 33, NCOLORS=254, BOTTOM=1
	TV, result, Startx, Starty, XSize=maxx*scale,YSize=maxy*scale
	LOADCT, 0, NCOLORS=256,BOTTOM=0
	MAP_GRID, /LABEL, /HORIZON, COLOR=255
	MAP_CONTINENTS, /coasts, COLOR=255, /COUNTRIES, /USA, /hires
	device, /close
end

pro Lv2toLv3, PathLv2, Lv3Path, Day, Month, Year, DatasetName, minvalue, maxvalue, resx, resy, lon1, lat1, lon2, lat2, winnumber, scale, show_plot_flag
	maxx=round((lon2-lon1)/resx)
	maxy=round((lat2-lat1)/resy)
	GenerateFilename, filename, Lv3Path, Day, Month, Year, Day, Month, Year, DatasetName, 'hdf5p', resx, resy, lon1, lat1, lon2, lat2, 0, ok
	gotdata=0
	if ok eq 1 then begin
		print, filename, ' already exists'
		return ;comment out this line if you want to recreate jpeg, ascii and hdf5 files
		ReadLv3plus, filename, DatasetName, DataGrid, DataGridWeight
		gotdata=1
	end
	if gotdata eq 0 then begin
		print, 'creating ',filename
		filename = PathLv2+"/"+string(Year,"$(i4.4)")+"/"+string(Month,"$(i2.2)")+"/"+string(Day,"$(i2.2)")+'/*.he5'
		print, 'searching for ',filename
		FileList = FILE_SEARCH(filename)
		NFiles = SIZE( FileList , /DIMENSIONS)
		if (NFiles(0) lt 15) and (NFiles(0) gt 0) then begin
			if (Year mod 4) eq 0 then DayOfTheMonth=[31,29,31,30,31,30,31,31,30,31,30,31] else DayOfTheMonth=[31,28,31,30,31,30,31,31,30,31,30,31]
			lYear=Year
			lMonth=Month
			lDay=Day
			if lDay eq DayOfTheMonth(lMonth-1) then begin
				lDay=1
				lMonth+=1
			end else lDay+=1
			if lMonth eq 13 then begin
				lMonth=1
				lYear+=1
			end
			filename2 = PathLv2+"/"+string(lYear,"$(i4.4)")+"/"+string(lMonth,"$(i2.2)")+"/"+string(lDay,"$(i2.2)")+'/*.he5'
			FileList2 = FILE_SEARCH(filename2)
			NFiles2 = SIZE( FileList2 , /DIMENSIONS)
			if NFiles2(0) gt 0 then FileList = [FileList,FileList2(0)]
			NFiles = SIZE( FileList , /DIMENSIONS)
		end
		if NFiles(0) gt 0 then begin
			for x = 0,NFiles(0)-1,1 do begin
				ReadLv2cFiles, FileList(x), DatasetName, Latitude, Longitude, Data, CloudCover, AMF, SZA, DataError, satLat, satLon, satHeight, ok
				if ok eq 1 then begin
					Dimensions = SIZE( Longitude , /DIMENSIONS)
					if (Dimensions(0) eq 60) and (Longitude(Dimensions(0)-1,Dimensions(1)/2) gt -1.0e10) then begin
						Coord = CornerCoordinates(Latitude,Longitude, satLat, satLon, satHeight)
						Add2Grid, DataGrid, DataGridWeight, Data, DataError, CloudCover, AMF, SZA, Coord, lat1, lon1, lat2, lon2, resy, resx
					end
				end
			end
			NormGrid, DataGrid, DataGridWeight, maxx, maxy
			gotdata=1
		end else gotdata=0
	end
	if (gotdata eq 1) and (max(DataGrid) gt 0.0) then begin
		GenerateFilename, filename, Lv3Path, Day, Month, Year, Day, Month, Year, DatasetName, 'hdf5p', resx, resy, lon1, lat1, lon2, lat2, 1, ok
		WriteLv3Plus, filename, DatasetName, DataGrid, DataGridWeight
		if show_plot_flag eq 1 then ShowGrid, DataGrid, maxx, maxy, scale, lat1, lat2, lon1, lon2, DatasetName, minvalue, maxvalue, winnumber
		GenerateFilename, filename, Lv3Path, Day, Month, Year, Day, Month, Year, DatasetName, 'hdf5', resx, resy, lon1, lat1, lon2, lat2, 1, ok
		WriteLv3, filename, DatasetName, DataGrid
		if (lon1 eq -180.0) and (lon2 eq 180.0) and (lat1 eq -90.0) and (lat2 eq 90.0) and (resx eq 0.25) and (resy eq 0.25) then begin
			GenerateFilename, filename, Lv3Path, Day, Month, Year, Day, Month, Year, DatasetName, 'ASCII', resx, resy, lon1, lat1, lon2, lat2, 1, ok
			WriteASCII, filename, Year, Month, Day, DataGrid, DatasetName, 0, 2.7e16
		end
		DataGrid*=DataGridWeight
		SmoothGrid, DataGrid, 3, 0, 1
		SmoothGrid, DataGridWeight, 3, 0, 1
		NormGrid, DataGrid, DataGridWeight, maxx, maxy
		GenerateFilename, filename, Lv3Path, Day, Month, Year, Day, Month, Year, DatasetName, 'Jpeg', resx, resy, lon1, lat1, lon2, lat2, 1, ok
		if (lon1 eq -180.0) and (lon2 eq 180.0) and (lat1 eq -90.0) and (lat2 eq 90.0) and (resx eq 0.25) and (resy eq 0.25) then begin
			WriteJpeg, DataGrid, filename, DatasetName, Lv3Path,  minvalue, maxvalue, 1, 75
		end
		GenerateFilename, filename, Lv3Path, Day, Month, Year, Day, Month, Year, DatasetName, 'Jpeg_plain', resx, resy, lon1, lat1, lon2, lat2, 1, ok
		DataGrid*=DataGridWeight
		SmoothGrid, DataGrid, 2, 1, 1
		SmoothGrid, DataGridWeight, 2, 1, 1
		NormGrid, DataGrid, DataGridWeight, maxx, maxy
		WriteJpeg, DataGrid, filename, DatasetName, Lv3Path,  minvalue, maxvalue, 0, 75
	end
end