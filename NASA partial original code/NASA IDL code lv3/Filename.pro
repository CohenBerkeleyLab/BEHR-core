pro GenerateFilename, filename, Path, Day, Month, Year,  Day2, Month2, Year2, DatasetName, type, resx, resy, lon1, lat1, lon2, lat2, write, ok
;written by Mark Wenig
;generates a path and file name for a map from covered time and space, data type and resolution
	maxx=(lon2-lon1)/resx
	maxy=(lat2-lat1)/resy
	filename = Path+"/"+DatasetName+'/'
	if (lon1 eq -180.0) and (lon2 eq 180.0) and (lat1 eq -90.0) and (lat2 eq 90.0) then begin
		filename=filename+'global_'+STRTRIM(long(360.0/resx),1)+"x"+STRTRIM(long(180.0/resy),1)
	end else begin
		filename=filename+'regional_'
		if lat1 lt 0.0 then filename=filename+STRTRIM(long(-lat1),1)+"S"$
		else filename=filename+STRTRIM(long(lat1),1)+"N"

		if lon1 lt 0.0 then filename=filename+STRTRIM(long(-lon1),1)+"W-"$
		else filename=filename+STRTRIM(long(lon1),1)+"E-"

		if lat2 lt 0.0 then filename=filename+STRTRIM(long(-lat2),1)+"S"$
		else filename=filename+STRTRIM(long(lat2),1)+"N"

		if lon2 lt 0.0 then filename=filename+STRTRIM(long(-lon2),1)+"W"$
		else filename=filename+STRTRIM(long(lon2),1)+"E"
		filename=filename+'_'+STRTRIM(round(maxx),1)+"x"+STRTRIM(round(maxy),1)
	end
	if (Day eq Day2) and (Month eq Month2) and (Year eq Year2) then time='daily'$
	else if (Day eq 1) and (Day2 eq DaysInMonth(Month, Year)) and (Month eq Month2) and (Year eq Year2) then time='monthly'$
	else if (Day eq 1) and (Month eq 1) and (Day2 eq 31) and (Month2 eq 12) and (Year eq Year2) then time='annual'$
	else time='average'
	filename=filename+"/"+time+'/'+type+'/'
	if (time ne 'annual') and (time ne 'monthly') then filename=filename+string(Year,"$(i4.4)")+'/'
	if write eq 1 then FILE_MKDIR, filename
	if time eq 'daily' then filename=filename+string(Year,"$(i4.4)")+"_"+string(Month,"$(i2.2)")+"_"+string(Day,"$(i2.2)")$
	else if time eq 'monthly' then  filename=filename+string(Year,"$(i4.4)")+"_"+string(Month,"$(i2.2)")$
	else if time eq 'annual' then  filename=filename+string(Year,"$(i4.4)")$
	else filename=filename+string(Year,"$(i4.4)")+"_"+string(Month,"$(i2.2)")+"_"+string(Day,"$(i2.2)")+'-'+string(Year2,"$(i4.4)")+"_"+string(Month2,"$(i2.2)")+"_"+string(Day2,"$(i2.2)")
	filename=filename+'_'+DatasetName
	if type eq 'ASCII' then filename=filename+'.txt'$
	else if type eq 'hdf5' then filename=filename+'.hdf5'$
	else if type eq 'hdf5p' then filename=filename+'.hdf5'$
	else if type eq 'Jpeg' then filename=filename+'.jpg'$
	else if type eq 'Jpeg_plain' then filename=filename+'.jpg'$
	else if type eq 'Jpeg_plain_Weight' then filename=filename+'.jpg'$
	else if type eq 'Mpeg' then filename=filename+'.mpg'$
	else if type eq 'kml' then filename=filename+'.kml'$
	else if type eq 'kmz' then filename=filename+'.kmz'$
	else if type eq 'html' then filename=filename+'.html'$
	else if type eq 'Png' then filename=filename+'.png'$
	else filename=filename+'.'+type
	ok=0
	if write eq 0 then begin
		ok = FILE_TEST( filename , /READ)
	end
end
