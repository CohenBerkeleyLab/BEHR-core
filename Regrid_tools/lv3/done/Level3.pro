pro Level3
;written by Mark Wenig

print
print, "---start------------------------------------------------------------------------------------------------------------------------------------"
print

!EXCEPT=0


On_OMI_Cluster=0



if On_OMI_Cluster eq 1 then begin
	show_plot_flag=0
	Lv2Path='/omi/live/ecs2/OMNO2'
	Lv3Path='/omi/home/wenig/Lv2toLv3/Data'
	LUTFileName=Lv3Path+'/LUTs/OMI2.tif'
	JpegPath='/omi/home/wenig/Lv2toLv3/Data/Jpegs'
	AsciiPath='/omi/home/wenig/Lv2toLv3/Data/ASCII'
end else if On_OMI_Cluster eq 0 then begin
	show_plot_flag=0
	Lv2Path='C:/Data/OMI/OMINO2/Lv2'
	Lv3Path='C:/Data/OMI/OMINO2/Lv3'
	LUTFileName=Lv3Path+'/LUTs/OMI2.tif'
	AsciiPath=Lv3Path+'/ASCII'
	JpegPath=Lv3Path+'/Jpegs'
end

GridAllData, 29,9,2004, 30, 9, 2007, 'NO2TropCS30', Lv2Path, Lv3Path, 0.25, 1.0e16
GridAllData, 29,9,2004, 30, 9, 2007, 'NO2TotalCS30', Lv2Path, Lv3Path, 0.25, 1.0e16

print
print, "----end--------------------------------------------------------------------------------------------------------------------------------"
print
end


