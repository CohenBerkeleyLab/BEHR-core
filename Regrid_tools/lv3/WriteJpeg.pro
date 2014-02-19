;written by Mark Wenig
pro AddContinentLines, colorimage, Lv3Path
	Dimensions = SIZE( colorimage , /DIMENSIONS)
	if (Dimensions(1) eq 1440) then Worldmap = READ_TIFF(Lv3Path+'/LUTs/Worldmap1440(continents).tif')
	colorimagesmooth=colorimage
	for y = 0,Dimensions(2)-1,1 do begin
		for x = 0,Dimensions(1)-1,1 do begin
			if (colorimage(0,x,y) eq 255) and (colorimage(1,x,y) eq 255) and (colorimage(2,x,y) eq 255) then begin
				colorimagesmooth(0,x,y)=128
				colorimagesmooth(1,x,y)=128
				colorimagesmooth(2,x,y)=128
			end
		end
	end
	SmoothGrid, colorimagesmooth, 10, 10, 2
	for y = 0,Dimensions(2)-1,1 do begin
		for x = 0,Dimensions(1)-1,1 do begin
			if (colorimage(0,x,y) eq 255) and (colorimage(1,x,y) eq 255) and (colorimage(2,x,y) eq 255) then begin
				colorimagesmooth(0,x,y)=255
				colorimagesmooth(1,x,y)=255
				colorimagesmooth(2,x,y)=255
			end
			if (Worldmap(x,y) ne 0) then begin
				br=colorimagesmooth(0,x,y)*3/8+colorimagesmooth(1,x,y)/2+colorimagesmooth(2,x,y)/8
				if br gt 80 then br=0 else br=255
				colorimage(0,x,y)=br
				colorimage(1,x,y)=br
				colorimage(2,x,y)=br
			end
		end
	end
end

pro ConvertDataToImage, Data, colorimage, ColorLUT, Lv3Path,  minvalue, maxvalue
	Dimensions = SIZE( Data , /DIMENSIONS)
	colorimage=BYTARR(3, Dimensions(0), Dimensions(1))
	image = BYTSCL( Data , MAX=maxvalue , MIN=minvalue )
	for y = 0,Dimensions(1)-1,1 do begin
		for x = 0,Dimensions(0)-1,1 do begin
			if Data(x,y) ne 0 then begin
				colorimage(0,x,Dimensions(1)-1-y)=ColorLUT(0,image(x,y))
				colorimage(1,x,Dimensions(1)-1-y)=ColorLUT(1,image(x,y))
				colorimage(2,x,Dimensions(1)-1-y)=ColorLUT(2,image(x,y))
			end else begin
				colorimage(0,x,Dimensions(1)-1-y)=255
				colorimage(1,x,Dimensions(1)-1-y)=255
				colorimage(2,x,Dimensions(1)-1-y)=255
			end
		end
	end
end

pro WriteJpeg, Data, FileName, DatasetName, Lv3Path,  minvalue, maxvalue, WithContours, Qual
	Dimensions = SIZE( Data , /DIMENSIONS)
	image = BYTSCL( Data , MAX=maxvalue , MIN=minvalue )
	ColorLUT = READ_TIFF(Lv3Path+'/LUTs/OMI2.tif')
	Worldmap = READ_TIFF(Lv3Path+'/LUTs/Worldmap1440(continents).tif')
	colorimagesmooth=BYTARR(3, Dimensions(0), Dimensions(1))
	ConvertDataToImage, Data, colorimage, ColorLUT, Lv3Path,  minvalue, maxvalue
	if (WithContours eq 1) then AddContinentLines, colorimage, Lv3Path
	print, 'writing ', FileName
	write_jpeg, FileName, colorimage, QUALITY=Qual,TRUE=1, /ORDER
end

pro WriteTransparentPng, Data, FileName, DatasetName, Lv3Path,  minvalue, maxvalue, LutName
	Dimensions = SIZE( Data , /DIMENSIONS)
	image = BYTSCL( Data , MAX=maxvalue , MIN=minvalue )
	ColorLUT = READ_TIFF(Lv3Path+'/LUTs/'+LutName)
	colorimagesmooth=BYTARR(3, Dimensions(0), Dimensions(1))
	ConvertDataToImage, Data, colorimage, ColorLUT, Lv3Path,  minvalue, maxvalue
	print, 'writing ', FileName
	transpmin=100
	transp=INDGEN(transpmin)
	for y = 0,Dimensions(1)-1,1 do begin
		for x = 0,Dimensions(0)-1,1 do begin
			if (image(x,y) lt transpmin) then image(x,y)=0
		end
	end
	WRITE_PNG, FileName, image, ColorLUT(0,*), ColorLUT(1,*), ColorLUT(2,*), /VERBOSE, TRANSPARENT=transp
end

pro WriteJpegWithClouds, Data, Clouds, FileName, DatasetName, Lv3Path,  minvalue, maxvalue, WithContours, Qual
	Dimensions = SIZE( Data , /DIMENSIONS)
	print, 'converting'
	image = BYTSCL( Data , MAX=maxvalue , MIN=minvalue )
	ColorLUT = READ_TIFF(Lv3Path+'/LUTs/OMI2.tif')
	GrayLUTT = READ_TIFF(Lv3Path+'/LUTs/grau.tif')
	Worldmap = READ_TIFF(Lv3Path+'/LUTs/Worldmap1440(continents).tif')
	ConvertDataToImage, Data, colorimage, ColorLUT, Lv3Path,  minvalue, maxvalue
	ConvertDataToImage, Clouds, cloudimage, GrayLUTT, Lv3Path,  0.0, 1.0
	th=0.0
	for y = 0,Dimensions(1)-1,1 do begin
		for x = 0,Dimensions(0)-1,1 do begin
			if Clouds(x,y) gt th then begin
				w=(Clouds(x,y)-th)/(1.0-th)
				colorimage(0,x,Dimensions(1)-1-y)=w*cloudimage(0,x,Dimensions(1)-1-y)+(1.0-w)*colorimage(0,x,Dimensions(1)-1-y)
				colorimage(1,x,Dimensions(1)-1-y)=w*cloudimage(1,x,Dimensions(1)-1-y)+(1.0-w)*colorimage(1,x,Dimensions(1)-1-y)
				colorimage(2,x,Dimensions(1)-1-y)=w*cloudimage(2,x,Dimensions(1)-1-y)+(1.0-w)*colorimage(2,x,Dimensions(1)-1-y)
			end
		end
	end
	if (WithContours eq 1) then AddContinentLines, colorimage, Lv3Path
	print, 'writing ', FileName
	write_jpeg, FileName, colorimage, QUALITY=Qual,TRUE=1, /ORDER
end

