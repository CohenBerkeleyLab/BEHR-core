;written by Mark Wenig
function ReadFromLv2, File_ID, Path, name
	ID = H5D_OPEN(File_ID, Path+name)
	data = H5D_READ(ID)
	H5D_CLOSE, ID
	return, data
end

pro ReadLv2cFiles, lFileName, lDatasetName, lat, lon, Data, cc, amfunpol ,sza, DataError, SpaceLat, SpaceLon, SpaceHigh, ok
	ok=1
	print, 'reading ',lDatasetName,' from ',lFileName
	if (lDatasetName eq 'SO205KM') then begin
		Path='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Data Fields/'
		PathGeo='/HDFEOS/SWATHS/OMI Total Column Amount SO2/Geolocation Fields/'
	end else begin
		Path='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/'
		PathGeo='/HDFEOS/SWATHS/ColumnAmountNO2/Geolocation Fields/'
	end
	File_ID = H5F_OPEN(lFileName)
	lat = ReadFromLv2(File_ID, PathGeo, 'Latitude')
	lon = ReadFromLv2(File_ID, PathGeo, 'Longitude')
	sza = ReadFromLv2(File_ID, PathGeo, 'SolarZenithAngle')
	SpaceLat = ReadFromLv2(File_ID, PathGeo, 'SpacecraftLatitude')
	SpaceLon = ReadFromLv2(File_ID, PathGeo, 'SpacecraftLongitude')
	SpaceHigh = ReadFromLv2(File_ID, PathGeo, 'SpacecraftAltitude')
	Dimensions = SIZE( lat , /DIMENSIONS)
	if (check_if_new_lv1_file(File_ID) eq 0) then begin
		ok=0
		return
	end
	if (lDatasetName eq 'SO205KM') then begin
		cc = ReadFromLv2(File_ID, Path, 'CloudFraction')
		cc*=1000
		Flag=INTARR(Dimensions)
	end else begin
		amfunpol = ReadFromLv2(File_ID, Path, 'AMFUnpolluted')
		cc = ReadFromLv2(File_ID, Path, 'CloudRadianceFraction')
		Flag = ReadFromLv2(File_ID, Path, 'vcdQualityFlags')
	end
	CASE lDatasetName OF
	'TerrainReflectivity': begin
		Data = ReadFromLv2(File_ID, Path, lDatasetName)
		DataError=fltarr(Dimensions)
		DataError+=1.0
  end
  'TerrainReflectivitySmooth': begin
		Data = ReadFromLv2(File_ID, Path, 'TerrainReflectivity')
		SmoothGrid, Data, 1, 1, 3
		DataError=fltarr(Dimensions)
		DataError+=1.0
  end
	'TerrainHeight': begin
		Data = ReadFromLv2(File_ID, Path, lDatasetName)
		DataError=fltarr(Dimensions)
		DataError+=1.0
  end
	'TerrainPressure': begin
		Data = ReadFromLv2(File_ID, Path, lDatasetName)
		DataError=fltarr(Dimensions)
		DataError+=1.0
  end
	'CloudRadianceFraction': begin
		Data = cc/1000.0
		DataError=fltarr(Dimensions)
		for x = 0,Dimensions(0)-1,1 do begin
			for y = 0,Dimensions(1)-1,1 do begin
				DataError(x,y)=(1.0+cc(x,y)/1000.0*3.0)
			end
		end
  end
	'SlantColumnAmountNO2Std': begin
		Data = ReadFromLv2(File_ID, Path, lDatasetName)
		DataError=fltarr(Dimensions)
		DataError+=1.0
  end
  'SlantColumnAmountNO2Col2': begin
		Data = ReadFromLv2(File_ID, Path, 'SlantColumnAmountNO2')
		DataError=fltarr(Dimensions)
		DataError+=1.0
  end
	'SlantColumnAmountNO2Col3': begin
		Data = ReadFromLv2(File_ID, Path, 'SlantColumnAmountNO2')
		DataError=fltarr(Dimensions)
		DataError+=1.0
  end
	'SlantColumnAmountNO2RelErr': begin
		Data = ReadFromLv2(File_ID, Path, 'SlantColumnAmountNO2')
		DataError = ReadFromLv2(File_ID, Path, 'SlantColumnAmountNO2Std')
		for x = 0,Dimensions(0)-1,1 do begin
			for y = 0,Dimensions(1)-1,1 do begin
				Data(x,y)=DataError(x,y)/Data(x,y)
				DataError(x,y)=1.0
			end
		end
	end

 'NO2Total': begin;
    print, 'calculating NO2Total'
	Data = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2')
	Data2 = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2BelowCloud')

	DataError = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2Std')
	DataError2 = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2BelowCloudStd')
	unpolFldCoefs=ReadFromLv2(File_ID, Path, 'UnpolFldCoefficients')
	background = evalUFld(lon,  lat,  unpolFldCoefs)
	cloudFraction = ReadFromLv2(File_ID, Path, 'CloudFraction')
	cloudFractionStd=ReadFromLv2(File_ID, Path, 'CloudFractionStd')
	terrainReflectivity =	ReadFromLv2(File_ID, Path, 'TerrainReflectivity')
	amfPollutedClear = ReadFromLv2(File_ID, Path, 'AMFPollutedClear')
	amfPollutedCloudy= ReadFromLv2(File_ID, Path, 'AMFPollutedCloudy')
	amfPollutedToGround= ReadFromLv2(File_ID, Path, 'AMFPollutedToGround')
	amfPolluted= ReadFromLv2(File_ID, Path, 'AMFPolluted')
	columnAmountNo2= ReadFromLv2(File_ID, Path, 'ColumnAmountNO2')
	columnAmountNo2Polluted= ReadFromLv2(File_ID, Path, 'ColumnAmountNO2Polluted')
	SlantColumnAmountNO2Std= ReadFromLv2(File_ID, Path, 'SlantColumnAmountNO2Std')
	SlantColumnAmountNO2Std=SlantColumnAmountNO2Std/3.0
	rms=ReadFromLv2(File_ID, Path, 'RootMeanSquareErrorOfFit')

	OmiStd, terrainReflectivity,  cloudFraction,   cc,  amfunpol,                 $  ;;Inputs
        amfPollutedCloudy,    amfPollutedClear,   amfPolluted,      amfPollutedToGround,              $
        background, columnAmountNo2Polluted, columnAmountNo2, SlantColumnAmountNO2Std, $
        AMFUnpollutedStd, AMFPollutedStd, AMFPollutedToGroundStd, columnAmountNO2UnpollutedStd,       $  ;;Outputs
        columnAmountNO2PollutedStd,   columnAmountNO2Std,     columnAmountNO2TropStd,                 $
        columnAmountNO2BelowCloudStd, columnAmountNO2InitStd, columnAmountNO2Init, columnAmountNO2BelowCloud
	DataError=	columnAmountNO2TropStd
	DataError2= columnAmountNO2BelowCloudStd

	for x = 0,Dimensions(0)-1,1 do begin
	for y = 0,Dimensions(1)-1,1 do begin
		DataError(x,y)=1.5e15*(1.0+cc(x,y)/1000.0*3.0);

		if (Data2(x,y) gt 0) and (DataError2(x,y) gt 0) then begin
	 		Data(x,y)+=Data2(x,y)
	 	end
		if (rms(x,y) gt 0.0003) then begin
			Data(x,y)=0.0
			DataError(x,y)=0.0
		end
	 end
	end
  end
 'NO2TotalCS30': begin;
    print, 'calculating NO2Total'
		Data = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2')
		Data2 = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2BelowCloud')
		DataError = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2Std')
		DataError2 = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2BelowCloudStd')
		unpolFldCoefs=ReadFromLv2(File_ID, Path, 'UnpolFldCoefficients')
		background = evalUFld(lon,  lat,  unpolFldCoefs)
		cloudFraction = ReadFromLv2(File_ID, Path, 'CloudFraction')
		cloudFractionStd=ReadFromLv2(File_ID, Path, 'CloudFractionStd')
		terrainReflectivity =	ReadFromLv2(File_ID, Path, 'TerrainReflectivity')
		amfPollutedClear = ReadFromLv2(File_ID, Path, 'AMFPollutedClear')
		amfPollutedCloudy= ReadFromLv2(File_ID, Path, 'AMFPollutedCloudy')
		amfPollutedToGround= ReadFromLv2(File_ID, Path, 'AMFPollutedToGround')
		amfPolluted= ReadFromLv2(File_ID, Path, 'AMFPolluted')
		columnAmountNo2= ReadFromLv2(File_ID, Path, 'ColumnAmountNO2')
		columnAmountNo2Polluted= ReadFromLv2(File_ID, Path, 'ColumnAmountNO2Polluted')
		SlantColumnAmountNO2Std= ReadFromLv2(File_ID, Path, 'SlantColumnAmountNO2Std')
		SlantColumnAmountNO2Std=SlantColumnAmountNO2Std/3.0
		rms=ReadFromLv2(File_ID, Path, 'RootMeanSquareErrorOfFit')

		OmiStd, terrainReflectivity,  cloudFraction,   cc,  amfunpol,                 $  ;;Inputs
            amfPollutedCloudy,    amfPollutedClear,   amfPolluted,      amfPollutedToGround,              $
            background, columnAmountNo2Polluted, columnAmountNo2, SlantColumnAmountNO2Std, $
            AMFUnpollutedStd, AMFPollutedStd, AMFPollutedToGroundStd, columnAmountNO2UnpollutedStd,       $  ;;Outputs
            columnAmountNO2PollutedStd,   columnAmountNO2Std,     columnAmountNO2TropStd,                 $
            columnAmountNO2BelowCloudStd, columnAmountNO2InitStd, columnAmountNO2Init, columnAmountNO2BelowCloud

		DataError=	columnAmountNO2TropStd
		DataError2= columnAmountNO2BelowCloudStd

		for x = 0,Dimensions(0)-1,1 do begin
			for y = 0,Dimensions(1)-1,1 do begin
				DataError(x,y)=1.5e15*(1.0+cc(x,y)/1000.0*3.0);
				if (Data2(x,y) gt 0) and (DataError2(x,y) gt 0) then begin
			 		Data(x,y)+=Data2(x,y)
			 	end
				if (rms(x,y) gt 0.0003) or (cc(x,y) gt 300) then begin
					Data(x,y)=0.0
					DataError(x,y)=0.0
				end
			end
		end
  end
'NO2Trop': begin
    print, 'calculating NO2Trop_New'
		Data = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2Trop')
		Data2 = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2BelowCloud')
		DataError = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2TropStd')
		DataError2 = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2BelowCloudStd')
		unpolFldCoefs=ReadFromLv2(File_ID, Path, 'UnpolFldCoefficients')
		background = evalUFld(lon,  lat,  unpolFldCoefs)
		cloudFraction = ReadFromLv2(File_ID, Path, 'CloudFraction')
		cloudFractionStd=ReadFromLv2(File_ID, Path, 'CloudFractionStd')
		terrainReflectivity =	ReadFromLv2(File_ID, Path, 'TerrainReflectivity')
		amfPollutedClear = ReadFromLv2(File_ID, Path, 'AMFPollutedClear')
		amfPollutedCloudy= ReadFromLv2(File_ID, Path, 'AMFPollutedCloudy')
		amfPollutedToGround= ReadFromLv2(File_ID, Path, 'AMFPollutedToGround')
		amfPolluted= ReadFromLv2(File_ID, Path, 'AMFPolluted')
		columnAmountNo2= ReadFromLv2(File_ID, Path, 'ColumnAmountNO2')
		columnAmountNo2Polluted= ReadFromLv2(File_ID, Path, 'ColumnAmountNO2Polluted')
		SlantColumnAmountNO2Std= ReadFromLv2(File_ID, Path, 'SlantColumnAmountNO2Std')
		SlantColumnAmountNO2Std=SlantColumnAmountNO2Std/3.0

		OmiStd, terrainReflectivity,  cloudFraction,   cc,  amfunpol,                 $  ;;Inputs
            amfPollutedCloudy,    amfPollutedClear,   amfPolluted,      amfPollutedToGround,              $
            background, columnAmountNo2Polluted, columnAmountNo2, SlantColumnAmountNO2Std, $

            AMFUnpollutedStd, AMFPollutedStd, AMFPollutedToGroundStd, columnAmountNO2UnpollutedStd,       $  ;;Outputs
            columnAmountNO2PollutedStd,   columnAmountNO2Std,     columnAmountNO2TropStd,                 $
            columnAmountNO2BelowCloudStd, columnAmountNO2InitStd, columnAmountNO2Init, columnAmountNO2BelowCloud

		DataError=	columnAmountNO2TropStd
		DataError2= columnAmountNO2BelowCloudStd
		rms=ReadFromLv2(File_ID, Path, 'RootMeanSquareErrorOfFit')
		for x = 0,Dimensions(0)-1,1 do begin
			for y = 0,Dimensions(1)-1,1 do begin
				DataError(x,y)=1.5e15*(1.0+cc(x,y)/1000.0*3.0);
				if (Data2(x,y) gt 0) and (DataError2(x,y) gt 0) then begin
			 		Data(x,y)+=Data2(x,y)
			 	end
				if (rms(x,y) gt 0.0003) then begin
					Data(x,y)=0.0
					DataError(x,y)=0.0
				end
			 end
		end
  end

  'NO2TropCS30': begin
    print, 'calculating NO2TropCS30'
		Data = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2Trop')
		Data2 = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2BelowCloud')
		DataError = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2TropStd')
		DataError2 = ReadFromLv2(File_ID, Path, 'ColumnAmountNO2BelowCloudStd')
		unpolFldCoefs=ReadFromLv2(File_ID, Path, 'UnpolFldCoefficients')
		background = evalUFld(lon,  lat,  unpolFldCoefs)
		cloudFraction = ReadFromLv2(File_ID, Path, 'CloudFraction')
		cloudFractionStd=ReadFromLv2(File_ID, Path, 'CloudFractionStd')
		terrainReflectivity =	ReadFromLv2(File_ID, Path, 'TerrainReflectivity')
		amfPollutedClear = ReadFromLv2(File_ID, Path, 'AMFPollutedClear')
		amfPollutedCloudy= ReadFromLv2(File_ID, Path, 'AMFPollutedCloudy')
		amfPollutedToGround= ReadFromLv2(File_ID, Path, 'AMFPollutedToGround')
		amfPolluted= ReadFromLv2(File_ID, Path, 'AMFPolluted')
		columnAmountNo2= ReadFromLv2(File_ID, Path, 'ColumnAmountNO2')
		columnAmountNo2Polluted= ReadFromLv2(File_ID, Path, 'ColumnAmountNO2Polluted')
		SlantColumnAmountNO2Std= ReadFromLv2(File_ID, Path, 'SlantColumnAmountNO2Std')
		SlantColumnAmountNO2Std=SlantColumnAmountNO2Std/3.0

		OmiStd, terrainReflectivity,  cloudFraction,   cc,  amfunpol,                 $  ;;Inputs
            amfPollutedCloudy,    amfPollutedClear,   amfPolluted,      amfPollutedToGround,              $
            background, columnAmountNo2Polluted, columnAmountNo2, SlantColumnAmountNO2Std, $

            AMFUnpollutedStd, AMFPollutedStd, AMFPollutedToGroundStd, columnAmountNO2UnpollutedStd,       $  ;;Outputs
            columnAmountNO2PollutedStd,   columnAmountNO2Std,     columnAmountNO2TropStd,                 $
            columnAmountNO2BelowCloudStd, columnAmountNO2InitStd, columnAmountNO2Init, columnAmountNO2BelowCloud

		DataError=	columnAmountNO2TropStd
		DataError2= columnAmountNO2BelowCloudStd
		rms=ReadFromLv2(File_ID, Path, 'RootMeanSquareErrorOfFit')
		for x = 0,Dimensions(0)-1,1 do begin
			for y = 0,Dimensions(1)-1,1 do begin
				DataError(x,y)=1.5e15*(1.0+cc(x,y)/1000.0*3.0);
				if (Data2(x,y) gt 0) and (DataError2(x,y) gt 0) then begin
			 		Data(x,y)+=Data2(x,y)
			 	end
				if (rms(x,y) gt 0.0003) or (cc(x,y) gt 300) then begin
					Data(x,y)=0.0
					DataError(x,y)=0.0
				end
			end
		end
  end

 'SO205KM': begin
    print, 'calculating SO2ColumnAmount05KM'
		Data = ReadFromLv2(File_ID, Path, 'SO2ColumnAmount05KM')
		DataError = ReadFromLv2(File_ID, Path, 'SO2ColumnAmount05KMbrd')
		cloudFraction = ReadFromLv2(File_ID, Path, 'CloudFraction')
		for x = 0,Dimensions(0)-1,1 do begin
			for y = 0,Dimensions(1)-1,1 do begin
				if (Data(x,y) lt -10.0) then begin
					DataError(x,y)=0.
				end else begin
		 			DataError(x,y)=1.5e15*(1.0+cc(x,y)/1000.0*3.0);
		 		end
			end
		end
  end
ELSE: BEGIN
	Data = ReadFromLv2(File_ID, Path, lDatasetName)
	DataError = ReadFromLv2(File_ID, Path, lDatasetName+'Std')
	for x = 0,Dimensions(0)-1,1 do begin
		for y = 0,Dimensions(1)-1,1 do begin
			if (DataError(x,y)/Data(x,y) lt 0.05) then begin
				DataError(x,y)=Data(x,y)*0.05
			end
		 	if (cc(x,y) lt 0) and ((Flag(x,y) mod 2) eq 0) then begin
		 		PRINT, 'negative cloud cover, but flag not set at pixel ',x,',',y, 'cc(x,y)=', cc(x,y), 'Data(x,y)=',Data(x,y)
		 	end


		end
	end
END
ENDCASE

;checking flags
for x = 0,Dimensions(0)-1,1 do begin
	for y = 0,Dimensions(1)-1,1 do begin
	 	if ((Flag(x,y) mod 2) ne 0) then begin
	 		Data(x,y)=0.0
	 		DataError(x,y)=0.0
		end
	end
end
H5F_CLOSE, File_ID
end
