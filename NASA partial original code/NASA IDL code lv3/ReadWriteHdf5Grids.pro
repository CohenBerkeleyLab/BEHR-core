;written by Mark Wenig
pro WriteLv3, lFileName, lDatasetName, Data
	print, 'writing ',lDatasetName,' to ',lFileName
	File_ID = H5F_CREATE(lFileName)
	datatype_id = H5T_IDL_CREATE(Data)
	dataspace_id = H5S_CREATE_SIMPLE(size(Data,/DIMENSIONS))
	dataset_id = H5D_CREATE(File_ID,lDatasetName,datatype_id,dataspace_id)
	H5D_WRITE,dataset_id,Data
	H5D_CLOSE,dataset_id
	H5S_CLOSE,dataspace_id
	H5T_CLOSE,datatype_id
	H5F_CLOSE,File_ID
end

pro WriteLv3plus, lFileName, lDatasetName, Data, DataWeight
	print, 'writing ',lDatasetName,' to ',lFileName
	File_ID = H5F_CREATE(lFileName)
	datatype_id = H5T_IDL_CREATE(Data)
	dataspace_id = H5S_CREATE_SIMPLE(size(Data,/DIMENSIONS))
	dataset_id = H5D_CREATE(File_ID,lDatasetName,datatype_id,dataspace_id)
	H5D_WRITE,dataset_id,Data
	datatype_id = H5T_IDL_CREATE(DataWeight)
	dataspace_id = H5S_CREATE_SIMPLE(size(DataWeight,/DIMENSIONS))
	dataset_id = H5D_CREATE(File_ID,lDatasetName+'Weight',datatype_id,dataspace_id)
	H5D_WRITE,dataset_id,DataWeight
	H5D_CLOSE,dataset_id
	H5S_CLOSE,dataspace_id
	H5T_CLOSE,datatype_id
	H5F_CLOSE,File_ID
end

pro ReadLv3, FileName, DatasetName, Data
	print, 'reading ',FileName
	File_ID = H5F_OPEN(FileName)
 	Dataset_ID = H5D_OPEN(File_ID, DatasetName)
 	Data = H5D_READ(Dataset_ID)
	H5D_CLOSE, Dataset_ID
	H5F_CLOSE, File_ID
end

pro ReadLv3plusDate, Day, Month, Year, DatasetName, Data, DataWeight, Lv3Path, resx, resy, lon1, lat1, lon2, lat2
	GenerateFilename, FileName, Lv3Path, Day, Month, Year, Day, Month, Year, DatasetName, 'hdf5p', resx, resy, lon1, lat1, lon2, lat2, 0, ok
	if ok eq 1 then ReadLv3plus, FileName, DatasetName, Data, DataWeight else print, 'file ', FileName,' not found'
end

pro ReadLv3plus, FileName, DatasetName, Data, DataWeight
	print, 'reading ',FileName
	File_ID = H5F_OPEN(FileName)
 	Dataset_ID = H5D_OPEN(File_ID, DatasetName)
 	Data = H5D_READ(Dataset_ID)
	Dataset_ID = H5D_OPEN(File_ID, DatasetName+'Weight')
 	DataWeight = H5D_READ(Dataset_ID)
	H5F_CLOSE, File_ID
end