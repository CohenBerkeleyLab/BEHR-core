;written by Mark Wenig
pro WriteASCII, FileName, Year, Month, Day, Data, DatasetName, minvalue, maxvalue
	if (Year mod 4) eq 0 then SumDayOfTheMonth=[0,31,60,91,121,152,182,213,244,274,305,335] else SumDayOfTheMonth=[0,31,59,90,120,151,181,212,243,273,304,334]
	NameOfMonths=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
	Dimensions = SIZE( Data , /DIMENSIONS)
	Julian = SYSTIME( /JULIAN)
	CALDAT, Julian, GenMonth , GenDay , GenYear
	print, GenYear, GenMonth, GenDay
	OPENW, lun, FileName, /GET_LUN
	if Day ne 0 then begin
		PRINTF, lun, ' Day:',string((SumDayOfTheMonth(Month-1)+Day),"$(i4)"),' ',NameOfMonths(Month-1), string(Day,"$(i3)"), ', ',$
			string(Year,"$(i4.4)"),'    OMI/AURA   ',DatasetName, '     GEN:',string((GenYear mod 100),"$(i2.2)"),':',string((SumDayOfTheMonth(GenMonth-1)+GenDay),"$(i3.3)")," Asc LECT:01:45"
	end	else begin
		PRINTF, lun, 'Startday:',string((SumDayOfTheMonth(Month-1)),"$(i4)"),' ',NameOfMonths(Month-1), ', ',$
			string(Year,"$(i4.4)"),'    OMI/AURA   ',DatasetName, '     GEN:',string((GenYear mod 100),"$(i2.2)"),':',string((SumDayOfTheMonth(GenMonth-1)+GenDay),"$(i3.3)")," Asc LECT:01:45"
	end
	PRINTF, lun, 'Longitudes: ', STRTRIM(Dimensions(0), 1),' bins centered on ', STRTRIM(180.0-180.0/Dimensions(0), 1), ' W to ',STRTRIM(180.0-180.0/Dimensions(0), 1),' E (',string(FORMAT='(1F0.2)',360.0/Dimensions(0)), ' degrees steps)'
	PRINTF, lun, 'Latitudes : ', STRTRIM(Dimensions(1), 1),' bins centered on ', STRTRIM(90.0-90.0/Dimensions(1), 1), ' S to ',STRTRIM(90.0-90.0/Dimensions(1), 1),' N (',string(FORMAT='(1F0.2)',180.0/Dimensions(1)), ' degrees steps)'
	Data = Data < maxvalue*999.0/1000.0
	line=STRARR(Dimensions(0)+1)
	for y = 0, Dimensions(1)-1,1 do begin
		for x = 0,Dimensions(0)-1,1 do begin
			line(x)=string(round(Data(x,y)/maxvalue*1000.0),"$(i3)")
		end
		line(Dimensions(0))=' lat = '+ string(FORMAT='(1F10.6)',90-180.0/Dimensions(1)*(double(y)+0.5))
		PRINTF,lun, FORMAT = '(25(A))', line
	end
	CLOSE, lun
	FREE_LUN, lun
end

pro ReadLv3ASCII, data, FileName
	data=fltarr(25*57+15,720)
	openr,lun,FileName,error=err, /GET_LUN
	if err(0) eq 0 then begin
		readf,lun,data,format='(49x/,78x/,76x,720(57(/25i3),/15i3,18x))'
		close,lun
	end else print, 'error reading ', FileName
	data*=2.7e13
	FREE_LUN, lun
end