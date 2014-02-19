;written by Mark Wenig
function check_if_new_lv1_file, File_ID
	Path='/HDFEOS/SWATHS/ColumnAmountNO2/Data Fields/'
	qflagnotfound=0
	num=H5G_GET_NMEMBERS(File_ID, Path)
	for i=0, num-1, 1 do begin
		objname=H5G_GET_MEMBER_NAME(File_ID, Path, i)
		if objname eq 'vcdQualityFlags' then qflagnotfound=1
	end
	return, qflagnotfound
end

function DaysInMonth, month, year
	if (year mod 4) eq 0 then DaysInMonthArr=[31,29,31,30,31,30,31,31,30,31,30,31] else DaysInMonthArr=[31,28,31,30,31,30,31,31,30,31,30,31]
	return, 	DaysInMonthArr(month-1)
end

pro SmoothGrid, Data, nx, ny, type
	if type eq 1 then begin
		kernelx=fltarr(3,1)
		kernelx(0,0)=1
		kernelx(1,0)=100
		kernelx(2,0)=1
		sx=102
		kernely=fltarr(1,3)
		kernely(0,0)=1
		kernely(0,1)=100
		kernely(0,2)=1
		sy=102
	end else if type eq 2 then begin
		kernelx=fltarr(1,3,3)
		kernelx(0,0,1)=1
		kernelx(0,1,1)=2
		kernelx(0,2,1)=1
		kernelx(0,2,2)=0
		sx=4
		kernely=fltarr(1,1,3)
		kernely(0,0,0)=1
		kernely(0,0,1)=2
		kernely(0,0,2)=1
		sy=4
	end else if type eq 3 then begin
		kernelx=fltarr(3,1)
		kernelx(0,0)=1
		kernelx(1,0)=2
		kernelx(2,0)=1
		sx=4
		kernely=fltarr(1,3)
		kernely(0,0)=1
		kernely(0,1)=2
		kernely(0,2)=1
		sy=4
	end
	for x = 1,nx do Data=CONVOL(Data,kernelx,sx, /EDGE_TRUNCATE)
	for y = 1,ny do Data=CONVOL(Data,kernely,sy, /EDGE_TRUNCATE)
end