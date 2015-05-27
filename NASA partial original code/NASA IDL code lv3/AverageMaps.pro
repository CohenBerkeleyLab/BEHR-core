pro AverageMaps, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName, resx, resy, lon1, lat1, lon2, lat2, minvalue, maxvalue
	GenerateFilename, FileName, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName, 'hdf5p', resx, resy, lon1, lat1, lon2, lat2, 0, ok
	if ok eq 1 then begin
		print, Filename,' already exists'
		;return
		ReadLv3plus, FileName, DatasetName, AverageGrid, AverageWeight

	end else begin
		maxx=round((lon2-lon1)/resx)
		maxy=round((lat2-lat1)/resy)
		InitGrid, AverageGrid, AverageWeight, maxx, maxy
		;NBadDays=1
		;BadDays=INTARR(3,NBadDays)

		;openr,lun,Lv3Path+'/LUTs/baddayslist.txt',error=err, /GET_LUN
		;if err(0) eq 0 then begin
		;	readf,lun,NBadDays
		;	BadDays=INTARR(3,NBadDays)
		;	readf,lun,BadDays
		;	CLOSE, lun
		;	FREE_LUN, lun
		;end
		;BadDay=0
		for Year = StartYear ,EndYear, 1 do begin

			if Year eq StartYear then lStartMonth=StartMonth else  lStartMonth=1
			if Year eq EndYear then lEndMonth=EndMonth else  lEndMonth=12
			for Month = lStartMonth ,lEndMonth, 1 do begin

				if (Year eq StartYear) and (Month eq StartMonth) then lStartDay=StartDay else  lStartDay=1
				if (Year eq EndYear) and (Month eq EndMonth) then lEndDay=EndDay else lEndDay=DaysInMonth(Month, Year)


				for Day = lStartDay, lEndDay, 1 do begin
				print, 'date=',Day, Month, Year
					GenerateFilename, FileName, Lv3Path, Day, Month, Year, Day, Month, Year, DatasetName, 'hdf5p', resx, resy, lon1, lat1, lon2, lat2, 0, ok
					if ok ne 1 then print, FileName, ' not found'

					;while (JULDAY(BadDays(1,BadDay), BadDays(2,BadDay), BadDays(0,BadDay), 0,0,0)-JULDAY(Month, Day, Year, 0,0,0) lt 0.0) and (BadDay lt NBadDays-1) do begin
					;	BadDay+=1
					;end
					;if JULDAY(BadDays(1,BadDay), BadDays(2,BadDay), BadDays(0,BadDay), 0,0,0)-JULDAY(Month, Day, Year, 0,0,0) eq 0 then ok=0
					;if ok eq 0 then print, Day, Month, Year, ' is a bad day!'
					if ok eq 1 then begin
						ReadLv3plus, FileName, DatasetName, Data, DataWeight
						for x = 0,maxx-1 do begin
							for y = 0,maxy-1 do begin
								AverageGrid(x,y)+=Data(x,y)*DataWeight(x,y)
								AverageWeight(x,y)+=DataWeight(x,y)
							end
						end
					end
				end
			end
		end
		print, 'done averaging'
		NormGrid, AverageGrid, AverageWeight, maxx, maxy
	end
	;if (StartDay eq 1) and (StartMonth eq 1) and (EndDay eq 31) and (EndMonth eq 12)  then time='annual'$
	;else time=string(1+EndDay-StartDay+(EndMonth-StartMonth)*30,"$(i2.2)")+'days'
	print, 'writing data'
	GenerateFilename, FileName, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName, 'hdf5p', resx, resy, lon1, lat1, lon2, lat2, 1, ok
	WriteLv3Plus, FileName, DatasetName, AverageGrid, AverageWeight

	GenerateFilename, FileName, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName, 'hdf5', resx, resy, lon1, lat1, lon2, lat2, 1, ok
	;print, FileName
	WriteLv3, FileName, DatasetName, AverageGrid


	GenerateFilename, FileName, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName, 'ASCII', resx, resy, lon1, lat1, lon2, lat2, 1, ok
	WriteASCII, FileName, StartYear, StartMonth, 0, AverageGrid, DatasetName, 0, 2.7e16

	GenerateFilename, FileName, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName, 'Jpeg_plain', resx, resy, lon1, lat1, lon2, lat2, 1, ok
	WriteJpeg, AverageGrid, FileName, DatasetName, Lv3Path,  minvalue, maxvalue, 0, 75

	if resx eq 0.25 then begin
		GenerateFilename, FileName, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName, 'Jpeg', resx, resy, lon1, lat1, lon2, lat2, 1, ok
		WriteJpeg, AverageGrid, FileName, DatasetName, Lv3Path,  minvalue, maxvalue, 1, 75
	end

end



pro Lv3toMpeg, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName, scale, minvalue, maxvalue, winnumber
	maxx=1440
	maxy=720
	lon1= -180.0
	lat1= -90.0
	lon2= 180.0
	lat2= 90.0
	res=0.25

	mpegID = MPEG_OPEN([1440,720]);,QUALITY=100)
	t=0;
 	for Year = StartYear ,EndYear, 1 do begin

		if (Year mod 4) eq 0 then DayOfTheMonth=[31,29,31,30,31,30,31,31,30,31,30,31] else DayOfTheMonth=[31,28,31,30,31,30,31,31,30,31,30,31]

		if Year eq StartYear then lStartMonth=StartMonth else  lStartMonth=1
		if Year eq EndYear then lEndMonth=EndMonth else  lEndMonth=12
		for Month = lStartMonth ,lEndMonth, 1 do begin

			if (Year eq StartYear) and (Month eq StartMonth) then lStartDay=StartDay else  lStartDay=1
			if (Year eq EndYear) and (Month eq EndMonth) then lEndDay=EndDay else lEndDay=DayOfTheMonth(Month-1)


			for Day = lStartDay, lEndDay, 1 do begin

				GenerateFilename, FileName, Lv3Path, Day, Month, Year, Day, Month, Year, DatasetName, 'hdf5p', res, res, lon1, lat1, lon2, lat2, 0, ok
				if ok eq 1 then begin

					ReadLv3plus, FileName, DatasetName, Data, DataWeight
				;	if t lt 2 then begin
				;		Data1=Data2
				;		Data2=Data3
				;	end else if t eq 1 then begin
				;
				;	end

				;	AvData=Data1*0.25+Data2*0.5+Data3*0.25

					datestr=string(day,"$(i2.2)")+'.'+string(month,"$(i2.2)")+'.'+string(year,"$(i4.4)")
					GridToImage, Data, maxx, maxy, scale, -90.0, 90.0, -180.0, 180.0, DatasetName, minvalue, maxvalue, winnumber, im, datestr

					MPEG_PUT, mpegID, IMAGE=im, FRAME=t
					MPEG_PUT, mpegID, IMAGE=im, FRAME=t+1
					MPEG_PUT, mpegID, IMAGE=im, FRAME=t+2
					MPEG_PUT, mpegID, IMAGE=im, FRAME=t+3
					MPEG_PUT, mpegID, IMAGE=im, FRAME=t+4
					t=t+5
				end
			end
		end
	end

	GenerateFilename, FileName, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName, 'Mpeg', res, res, lon1, lat1, lon2, lat2, 1, ok
	print, 'writing ',FileName
	MPEG_SAVE, mpegID, FILENAME=FileName

	; Close the MPEG sequence:
	MPEG_CLOSE, mpegID
end

pro Convert2GMIres, Lv3Path, Day, Month, Year, Day2, Month2, Year2, DatasetName, minvalue, maxvalue
	;if (Year mod 4) eq 0 then DayOfTheMonth=[31,29,31,30,31,30,31,31,30,31,30,31] else DayOfTheMonth=[31,28,31,30,31,30,31,31,30,31,30,31]
	res=0.25
	GenerateFilename, FileName, Lv3Path, Day, Month, Year, Day2, Month2, Year2, DatasetName, 'hdf5p', res, res, -180.0, -90.0, 180.0, 90.0, 0, ok
	;FileName="C:\Data\OMI\OMINO2\Lv3\AGU06\res025x025\2005_04_"+DatasetName+".hdf5"
	ok=1
	if ok eq 1 then begin

		ReadLv3plus, FileName, DatasetName, Data, DataWeight
		maxx=360.0/res
		maxy=180.0/res

		newresx=2.5
		newresy=2.0
		newmaxx=360.0/newresx
		newmaxy=180.0/newresy+1
		NewData=fltarr(newmaxx,newmaxy)
		NewDataWeight=fltarr(newmaxx,newmaxy)
		NewDataNewWeight=fltarr(newmaxx,newmaxy)


		for x = 0,maxx-1 do begin
			for y = 0,maxy-1 do begin
				NewData((res/newresx*(x+0.5)),(res/newresy*y))+=Data(x,y)
				NewDataWeight((res/newresx*(x+0.5)),(res/newresy*y))+=DataWeight(x,y)
				NewDataNewWeight((res/newresx*(x+0.5)),(res/newresy*y))+=1.0
			end
		end

		for x = 0,newmaxx-1 do begin
			for y = 0,newmaxy-1 do begin
				if NewDataNewWeight(x,y) ne 0 then NewData(x,y)/=NewDataNewWeight(x,y)
				if NewDataNewWeight(x,y) ne 0 then NewDataWeight(x,y)/=NewDataNewWeight(x,y)
			end
		end
		GenerateFilename, FileName, Lv3Path, Day, Month, Year, Day2, Month2, Year2, DatasetName, 'hdf5p', newresx, newresy, -180.0, -90.0, 180.0, 90.0, 1, ok
		;FileName="C:\Data\OMI\OMINO2\Lv3\AGU06\res25x2\2005_04_"+DatasetName+".hdf5"
		WriteLv3Plus, filename, DatasetName, NewData, NewDataWeight

		GenerateFilename, filename, Lv3Path, Day, Month, Year, Day2, Month2, Year2, DatasetName, 'Jpeg_plain', newresx, newresy, -180.0, -90.0, 180.0, 90.0, 1, ok
		;FileName="C:\Data\OMI\OMINO2\Lv3\AGU06\res25x2\2005_04_"+DatasetName+".jpg"
		WriteJpeg, NewData, filename, DatasetName, Lv3Path,  minvalue, maxvalue, 0, 75


	end else begin
		print, 'error reading file'
	end
end

pro AverageWeekDayMaps, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName, resx, resy, lon1, lat1, lon2, lat2, minvalue, maxvalue, weekday
	WeekDayName=['Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat','Sun']
	maxx=round((lon2-lon1)/resx)
	maxy=round((lat2-lat1)/resy)
	InitGrid, AverageGrid, AverageWeight, maxx, maxy
	NBadDays=1
	BadDays=INTARR(3,NBadDays)

	openr,lun,Lv3Path+'/LUTs/baddayslist.txt',error=err, /GET_LUN
	if err(0) eq 0 then begin
		readf,lun,NBadDays
		BadDays=INTARR(3,NBadDays)
		readf,lun,BadDays
		CLOSE, lun
		FREE_LUN, lun
	end
	BadDay=0
	Day=StartDay
	Month=StartMonth
	Year=StartYear

	while (round(JULDAY(Month, Day, Year, 12,0,0) mod 7) ne weekday) do begin
		NextDay, Day, Month, Year
		print, Day, Month, Year, weekday,round(JULDAY(Month, Day, Year, 12,0,0) mod 7)
	end

	WHILE date1_before_date2(Day, Month, Year, EndDay, EndMonth, EndYear) DO BEGIN
  	print, Day, Month, Year

		GenerateFilename, FileName, Lv3Path, Day, Month, Year, Day, Month, Year, DatasetName, 'hdf5p', resx, resy, lon1, lat1, lon2, lat2, 0, ok
		if ok ne 1 then print, FileName, ' not found'

		while (JULDAY(BadDays(1,BadDay), BadDays(2,BadDay), BadDays(0,BadDay), 0,0,0)-JULDAY(Month, Day, Year, 0,0,0) lt 0.0) and (BadDay lt NBadDays-1) do begin
			BadDay+=1
		end
		if JULDAY(BadDays(1,BadDay), BadDays(2,BadDay), BadDays(0,BadDay), 0,0,0)-JULDAY(Month, Day, Year, 0,0,0) eq 0 then ok=0
		if ok eq 0 then print, Day, Month, Year, ' is a bad day!'
		if ok eq 1 then begin
			ReadLv3plus, FileName, DatasetName, Data, DataWeight
			for x = 0,maxx-1 do begin
				for y = 0,maxy-1 do begin
					AverageGrid(x,y)+=Data(x,y)*DataWeight(x,y)
					AverageWeight(x,y)+=DataWeight(x,y)
				end
			end
		end
		for x = 1,7 do NextDay, Day, Month, Year
	ENDWHILE

;return
	print, 'done averaging'
	NormGrid, AverageGrid, AverageWeight, maxx, maxy

	print, 'writing data'
	GenerateFilename, FileName, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName+STRTRIM(weekday,1), 'hdf5p', resx, resy, lon1, lat1, lon2, lat2, 1, ok
	WriteLv3Plus, FileName, DatasetName, AverageGrid, AverageWeight

	;GenerateFilename, FileName, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName+STRTRIM(weekday,1), 'hdf5', resx, resy, lon1, lat1, lon2, lat2, 1, ok
	;print, FileName
	;WriteLv3, FileName, DatasetName, AverageGrid


	;GenerateFilename, FileName, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName+STRTRIM(weekday,1), 'ASCII', resx, resy, lon1, lat1, lon2, lat2, 1, ok
	;WriteASCII, FileName, StartYear, StartMonth, 0, AverageGrid, DatasetName, 0, 2.7e16

	GenerateFilename, FileName, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName+STRTRIM(weekday,1), 'Jpeg_plain', resx, resy, lon1, lat1, lon2, lat2, 1, ok
	WriteJpeg, AverageGrid, filename, DatasetName, Lv3Path,  minvalue, maxvalue, 0, 75

	if resx eq 0.25 then begin
		GenerateFilename, FileName, Lv3Path, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DatasetName+STRTRIM(weekday,1), 'Jpeg', resx, resy, lon1, lat1, lon2, lat2, 1, ok
		WriteJpeg, AverageGrid, filename, DatasetName, Lv3Path,  minvalue, maxvalue, 1, 75
	end

end