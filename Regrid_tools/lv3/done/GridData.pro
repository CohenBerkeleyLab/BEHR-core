pro GridAllData, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DataSetName, Lv2Path, Lv3Path, res, maxvalue
;written by Mark Wenig
;grids data and creates global maps
	lon1= -180.0
	lat1= -90.0
	lon2= 180.0
	lat2= 90.0
	for Year = StartYear ,EndYear, 1 do begin
		if Year eq StartYear then lStartMonth=StartMonth else  lStartMonth=1
		if Year eq EndYear then lEndMonth=EndMonth else  lEndMonth=12
		for Month = lStartMonth ,lEndMonth, 1 do begin
			if (Year eq StartYear) and (Month eq StartMonth) then lStartDay=StartDay else  lStartDay=1
			if (Year eq EndYear) and (Month eq EndMonth) then lEndDay=EndDay else lEndDay=DaysInMonth(Month,Year)
			for Day = lStartDay, lEndDay, 1 do begin
				Lv2toLv3, Lv2Path,Lv3Path, Day, Month, Year, DataSetName, 0, maxvalue, res, res, lon1,lat1, lon2, lat2, 0, 1.0, 0
			end
		end
	end
end

pro GridArea, StartDay, StartMonth, StartYear, EndDay, EndMonth, EndYear, DataSetName, Lv2Path, Lv3Path, lon1, lat1, lon2, lat2, res, minvalue, maxvalue
;written by Mark Wenig
;grids data and creates maps
	for Year = StartYear ,EndYear, 1 do begin
		if Year eq StartYear then lStartMonth=StartMonth else  lStartMonth=1
		if Year eq EndYear then lEndMonth=EndMonth else  lEndMonth=12
		for Month = lStartMonth ,lEndMonth, 1 do begin
			if (Year eq StartYear) and (Month eq StartMonth) then lStartDay=StartDay else  lStartDay=1
			if (Year eq EndYear) and (Month eq EndMonth) then lEndDay=EndDay else lEndDay=DaysInMonth(Month, Year)
			for Day = lStartDay, lEndDay, 1 do begin
				Lv2toLv3, Lv2Path,Lv3Path, Day, Month, Year, DataSetName, 0, 1.0e16, res, res, lon1,lat1, lon2, lat2, 0, 1.0, 0
			end
		end
	end
end


