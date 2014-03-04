% readhe5_omno2_v_aug2012
%Reads omno2 he5 files as of the Aug 2012 version; saves the resulting .mat
%file as <satellite>_<retrieval>_<year><month><day>. Based on
%readhe5_neus_wcld by Ashley Russel.
%
%Josh Laughner <joshlaugh5@gmail.com> 27 Feb 2014

%Add the 'Utils' folder and all subfolders to MATLAB's search path. Within
%the Git repository for BEHR, this is the /Utils folder.
addpath(genpath('/Users/Josh/Documents/MATLAB/BEHR/Utils'))

%Specify the longitude and latitude ranges of interest for this retrieval.
%****************************%
lonmin = -125;  lonmax = -95;
latmin = 37.5;    latmax = 50;
%****************************%
if lonmin > lonmax
    error('read_omno2:maxmin','Lonmin is greater than lonmax')
elseif latmin > latmax
    error('read_omno2:maxmin', 'Latmin is greater than latmax')
end
    
%These will be included in the file name
%****************************%
satellite='OMI';
retrieval='SP';
%****************************%

%This is the directory where the final .mat file will be saved. This will
%need to be changed to match your machine and the files' location. Do not
%include a trailing separator, i.e. '/my/favorite/directory' not
%'my/favorite/directory/
mat_dir = '/Volumes/share/GROUP/SAT/BEHR/Test_SP_files';

%This is the directory where the he5 files are saved. Do not include a
%trailing separator.
he5_dir = '/Volumes/share/GROUP/SAT/OMI/OMNO2_32';

%This is the directory where the MODIS myd06_L2*.hdf files are saved. It should include subfolders organized by year.
%Do not include a trailing separator.
modis_myd06_dir = '/Volumes/share/GROUP/SAT/MODIS/MYD06_L2';

%Process all files between these dates, in yyyy/mm/d format
%****************************%
date_start='2007/02/01';
date_end='2007/02/28';
%****************************%

tic %Start the timer

%Initialize matrices to hold the OMI data
Latitude=zeros(60,2000);
Longitude=zeros(60,300);
SpacecraftAltitude=zeros(300,1);
SpacecraftLatitude=zeros(300,1);
SpacecraftLongitude=zeros(300,1);
Time=zeros(300,1);
ViewingZenithAngle=zeros(60,300);
SolarZenithAngle=zeros(60,300);
ViewingAzimuthAngle=zeros(60,300);
SolarAzimuthAngle=zeros(60,300);
AMFStrat=zeros(60,300);
AMFTrop=zeros(60,300);
CloudFraction=zeros(60,300);
CloudRadianceFraction=zeros(60,300);
ColumnAmountNO2=zeros(60,300);
SlantColumnAmountNO2=zeros(60,300);
TerrainHeight=zeros(60,300);
TerrainPressure=zeros(60,300);
TerrainReflectivity=zeros(60,300);
vcdQualityFlags=zeros(60,300);
CloudPressure=zeros(60,300);
ColumnAmountNO2Trop=zeros(60,300);


%File names will be prefixed with "<satellite>_<retrieval>_", e.g. for OMI
%satellite SP retrieval, the prefix will be "OMI_SP_" and then the date in
%year, month, date order.  This section checks to see if the last file in
%the mat directory has the expected prefix.  If so, that date is taken as
%the last date completed, otherwise it is assumed that the retrieval will
%need to start from the specified start date. This allows he5 reading to be
%stopped and restarted with minimal intervention.
last_file=dir(fullfile(mat_dir,'*.mat'));
file_prefix = [satellite,'_',retrieval,'_']; l = length(file_prefix);

if isempty(last_file) || ~strcmp(last_file(end).name(1:l),file_prefix);
    last_date=date_start;
else
    last_date=last_file(end).name((l+1):(l+8));
    last_date=([last_date(1:4),'/',last_date(5:6),'/',last_date(7:8)]);
end

%For loop over all days from the starting or last finished date to the end
%date. We will give the absolute paths to files rather than changing the
%active directory, as MATLAB seems to run slightly slower if the current
%working directory is on the server.
total_days=datenum(date_end)-datenum(last_date)+1;
for j=1:total_days;
    %Read the desired year, month, and day
    R=addtodate(datenum(last_date), j, 'day');
    date=datestr(R,26);
    year=date(1:4);
    month=date(6:7);
    day=date(9:10);
    
    %Prepare a data structure to receive the final data.
    Data=struct('Longitude',0,'Latitude',0,'Loncorn',0,'Latcorn',0,'Time',0,'ViewingZenithAngle',0,'SolarZenithAngle',0,'ViewingAzimuthAngle',0,'SolarAzimuthAngle',0,'AMFStrat',0,'AMFTrop',0,'CloudFraction',0,'CloudRadianceFraction',0,'TerrainHeight',0,'TerrainPressure',0,'TerrainReflectivity',0,'vcdQualityFlags',0,'CloudPressure',0,'ColumnAmountNO2',0,'SlantColumnAmountNO2',0,'ColumnAmountNO2Trop',0,'MODISCloud',0);
    
    %Set the file path and name, assuming that the file structure is
    %<he5_directory>/<year>/<month>/...files...  Then figure out how many
    %files there are
    short_filename=['OMI-Aura_L2-OMNO2_',year,'m',month,day,'*.he5'];
    file_dir = fullfile(he5_dir,year,month); %Used both here to find all he5 files and in the swath for loop to identify each file.
    file=fullfile(file_dir,short_filename);
    sp_files = dir(file);
    n = length(sp_files);
    E=0;
    if isempty(sp_files);
        disp(['No Data Available For ',month,' ',day,' ',year])
    else
        for e=1:n %For loop over all the swaths in a given day.
            if e==1 || mod(e,10)==0; fprintf('Swath %u of %s/%s/%s \n',e,month,day,year); end
            %Read in each file, saving the hierarchy as 'hinfo'
            filename= sp_files(e).name;
            hinfo = h5info(fullfile(file_dir,filename));
            
            %Read in the full latitude data set; this will be used to determine
            %which pixels to read in later.
            Latitude = h5read(fullfile(file_dir,filename), h5dsetname(hinfo,1,2,1,2,'Latitude')); %h5dsetname takes 1) the object returned by h5info, 2) The indicies of the group tree 3) The last argument may be the index or name of the dataset of interest
            Row=0:59; Row=Row'; Row=repmat(Row,1,size(Latitude,2));
            Swath=filename(35:39); Swath=str2double(Swath).*ones(size(Latitude));
            
            %Restrict latitude to those that fall within the bounds specified
            %at the begininning of the file. Also pivot the dataset so that
            %each row is a swath.
            lat=Latitude';
            lat_i=[latmin, latmax];
            [i_i, j_j]=find(lat > lat_i(1) - 0.25 & lat < lat_i(2) + 0.25);
            cut_y=min(i_i):max(i_i);
            cut_x = 1:60;
            lat=double(lat(cut_y,cut_x));
            Latitude=Latitude(cut_x,cut_y)'; Latitude=double(Latitude);
            Row=Row(cut_x,cut_y)';
            Swath=Swath(cut_x,cut_y)';
            
            %Set up to use low-level HDF5 functions to read in only the parts
            %of the data set that fall within latitude boundaries (to save
            %memory).
            stride = [];
            blocksize = [];
            offset = [(min(i_i)-1),0];
            slabsize = [length(cut_y),60];
            memspaceID = H5S.create_simple(length(slabsize), slabsize, slabsize);
            fileID = H5F.open(fullfile(file_dir,filename), 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
            
            %This will handle each of the variables that are 60x(number of
            %swaths).  It also converts all data from single precision to
            %double precision and pivots the matrix the the convention of row =
            %swath. These are the values needed to compute the pixel corner
            %points.
            
            %Longitude
            datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'Longitude')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); Longitude = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); Longitude=double(Longitude); Longitude=Longitude';
            %ViewingAzimuthAngle
            datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'ViewingAzimuthAngle')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ViewingAzimuthAngle = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ViewingAzimuthAngle=double(ViewingAzimuthAngle); ViewingAzimuthAngle=ViewingAzimuthAngle';
            %ViewingZenithAngle
            datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'ViewingZenithAngle')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ViewingZenithAngle = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ViewingZenithAngle=double(ViewingZenithAngle); ViewingZenithAngle=ViewingZenithAngle';
            %SolarAzimuthAngle
            datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'SolarAzimuthAngle')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SolarAzimuthAngle = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SolarAzimuthAngle=double(SolarAzimuthAngle); SolarAzimuthAngle=SolarAzimuthAngle';
            %SolarZenithAngle
            datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'SolarZenithAngle')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SolarZenithAngle = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SolarZenithAngle=double(SolarZenithAngle); SolarZenithAngle=SolarZenithAngle';
            
            
            %This will handle values that only have a single value per swath.
            %They are still converted to double precision numbers and pivoted.
            offset = [(min(i_i)-1)];
            slabsize = [length(cut_y)];
            memspaceID = H5S.create_simple(length(slabsize), slabsize, slabsize);
            
            %SpacecraftAltitude
            datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'SpacecraftAltitude')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SpacecraftAltitude = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SpacecraftAltitude=double(SpacecraftAltitude); SpacecraftAltitude=SpacecraftAltitude';
            %SpacecraftLatitude
            datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'SpacecraftLatitude')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SpacecraftLatitude = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SpacecraftLatitude=double(SpacecraftLatitude); SpacecraftLatitude=SpacecraftLatitude';
            %SpacecraftLongitude
            datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'SpacecraftLongitude')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SpacecraftLongitude = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SpacecraftLongitude=double(SpacecraftLongitude); SpacecraftLongitude=SpacecraftLongitude';
            %Time
            datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'Time')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); Time = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); Time=double(Time); Time=Time';
            Time=repmat(Time,1,60);
            
            
            
            %Deletes any point that falls outside of the boundaries specified.
            Lat=Latitude; Lon=Longitude;
            x=find(Lon>lonmax | Lon<lonmin);
            y=find(Lat>latmax | Lat<latmin);
            Lon(x)=NaN;     Lon(y)=NaN;     Lon(isnan(Lon))=[];
            Lat(x)=NaN;     Lat(y)=NaN;     Lat(isnan(Lat))=[];
            
            if isempty(Lon)==1 || isempty(Lat)==1 || length(Lat)==1;
                disp('No points within lat/lon boundaries')
                continue
            else
                disp('Founds points within lat/lon boundaries')
                corners = fxn_corner_coordinates(Latitude, Longitude, SpacecraftLatitude, SpacecraftLongitude, SpacecraftAltitude);
                E=E+1;
                lat = corners(:,:,2,5); %Assign the center of each pixel to lat and lon
                lon = corners(:,:,1,5);
                latcorn = corners(:,:,2,1:4); latcorn = squeeze(latcorn);
                a = latcorn(:,:,1); a = a(:); b = latcorn(:,:,2); b = b(:); c = latcorn(:,:,3); c = c(:); d = latcorn(:,:,4); d = d(:);
                latcorn = [a,b,c,d]; latcorn = latcorn';
                loncorn = corners(:,:,1,1:4); loncorn = squeeze(loncorn);
                a = loncorn(:,:,1); a = a(:); b = loncorn(:,:,2); b = b(:); c = loncorn(:,:,3); c = c(:); d = loncorn(:,:,4); d = d(:);
                loncorn = [a,b,c,d]; loncorn = loncorn';
                
                %Import the FoV75 corner lat and lons.  These will be
                %ordered the same as the BEHR-calculated corners, i.e. corner x
                %along track x across track
                slabsize = [length(cut_y),60,4];
                memspaceID = H5S.create_simple(length(slabsize),slabsize,slabsize);
                offset = [(min(i_i)-1),0,0];
                
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'FoV75CornerLatitude')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); FoV75CornerLatitude = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); FoV75CornerLatitude = double(FoV75CornerLatitude); FoV75CornerLatitude = permute(FoV75CornerLatitude, [2 3 1]);
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'FoV75CornerLongitude')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); FoV75CornerLongitude = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); FoV75CornerLongitude = double(FoV75CornerLongitude); FoV75CornerLongitude = permute(FoV75CornerLongitude, [2 3 1]);
                
                %Import all remaining pieces of information from the standard
                %product.
                offset = [(min(i_i)-1),0];
                slabsize = [length(cut_y),60];
                memspaceID = H5S.create_simple(length(slabsize), slabsize, slabsize);
                
                %AMFStratsphere
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'AmfStrat')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); AMFStrat = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); AMFStrat=double(AMFStrat); AMFStrat=AMFStrat';
                %AMFTroposphere
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'AmfTrop')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); AMFTrop = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); AMFTrop=double(AMFTrop); AMFTrop=AMFTrop';
                %CloudFraction
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'CloudFraction')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudFraction = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudFraction=double(CloudFraction); CloudFraction=CloudFraction';
                %CloudFractionError
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'CloudFractionStd')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudFractionError = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudFractionError=double(CloudFractionError); CloudFractionError=CloudFractionError';
                %CloudPressure
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'CloudPressure')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudPressure = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudPressure=double(CloudPressure); CloudPressure=CloudPressure';
                %CloudPressureError
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'CloudPressureStd')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudPressureError = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudPressureError=double(CloudPressureError); CloudPressureError=CloudPressureError';
                %CloudRadianceFraction
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1, 'CloudRadianceFraction')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudRadianceFraction = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudRadianceFraction=double(CloudRadianceFraction); CloudRadianceFraction=CloudRadianceFraction';
                %ColumnAmountNO2
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'ColumnAmountNO2')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ColumnAmountNO2 = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ColumnAmountNO2=double(ColumnAmountNO2); ColumnAmountNO2=ColumnAmountNO2';
                %ColumnAmountNO2Trop
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'ColumnAmountNO2Trop')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ColumnAmountNO2Trop = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ColumnAmountNO2Trop=double(ColumnAmountNO2Trop); ColumnAmountNO2Trop=ColumnAmountNO2Trop';
                %SlantColumnAmountNO2
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'SlantColumnAmountH2O')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SlantColumnAmountNO2 = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SlantColumnAmountNO2=double(SlantColumnAmountNO2); SlantColumnAmountNO2=SlantColumnAmountNO2';
                %TerrainHeight
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'TerrainHeight')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); TerrainHeight = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); TerrainHeight=double(TerrainHeight); TerrainHeight=TerrainHeight';
                %TerrainPressure
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'TerrainPressure')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); TerrainPressure = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); TerrainPressure=double(TerrainPressure); TerrainPressure=TerrainPressure';
                %TerrainReflectivity
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'TerrainReflectivity')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); TerrainReflectivity = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); TerrainReflectivity=double(TerrainReflectivity); TerrainReflectivity=TerrainReflectivity';
                %vcdQualityFlags
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'VcdQualityFlags')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); vcdQualityFlags = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); vcdQualityFlags=double(vcdQualityFlags); vcdQualityFlags=vcdQualityFlags';
                
                
                H5F.close(fileID); %close omi file to free up space
                
                s=size(SolarAzimuthAngle);
                latcorn=reshape(latcorn,4,s(1),s(2));
                loncorn=reshape(loncorn,4,s(1),s(2));
                Time=reshape(Time,s(1),s(2));
                
                RelativeAzimuthAngle=abs(SolarAzimuthAngle+180-ViewingAzimuthAngle);
                RelativeAzimuthAngle(RelativeAzimuthAngle > 180)=360-RelativeAzimuthAngle(RelativeAzimuthAngle > 180);
                
                loncorn=loncorn(1:4,:);
                latcorn=latcorn(1:4,:);
                
                
                x=find(lon<lonmin | lon>lonmax);
                y=find(lat<latmin | lat>latmax);
                lon(x)=NaN;                         lon(y)=NaN;                         lon(isnan(lon))=[];
                lat(x)=NaN;                         lat(y)=NaN;                         lat(isnan(lat))=[];
                loncorn(:,x)=NaN;                   loncorn(:,y)=NaN;                   loncorn=loncorn';                   loncorn(any(isnan(loncorn)'),:) = [];               loncorn=loncorn';
                latcorn(:,x)=NaN;                   latcorn(:,y)=NaN;                   latcorn=latcorn';                   latcorn(any(isnan(latcorn)'),:) = [];               latcorn=latcorn';
                SolarAzimuthAngle(x)=NaN;           SolarAzimuthAngle(y)=NaN;           SolarAzimuthAngle(isnan(SolarAzimuthAngle))=[];
                SolarZenithAngle(x)=NaN;            SolarZenithAngle(y)=NaN;            SolarZenithAngle(isnan(SolarZenithAngle))=[];
                ViewingAzimuthAngle(x)=NaN;         ViewingAzimuthAngle(y)=NaN;         ViewingAzimuthAngle(isnan(ViewingAzimuthAngle))=[];
                ViewingZenithAngle(x)=NaN;          ViewingZenithAngle(y)=NaN;          ViewingZenithAngle(isnan(ViewingZenithAngle))=[];
                Time(x)=NaN;                        Time(y)=NaN;                        Time(isnan(Time))=[];
                AMFStrat(x)=NaN;                    AMFStrat(y)=NaN;                    AMFStrat(isnan(AMFStrat))=[];
                AMFTrop(x)=NaN;                     AMFTrop(y)=NaN;                     AMFTrop(isnan(AMFTrop))=[];
                CloudFraction(x)=NaN;               CloudFraction(y)=NaN;               CloudFraction(isnan(CloudFraction))=[];
                CloudPressure(x)=NaN;               CloudPressure(y)=NaN;               CloudPressure(isnan(CloudPressure))=[];
                CloudRadianceFraction(x)=NaN;       CloudRadianceFraction(y)=NaN;       CloudRadianceFraction(isnan(CloudRadianceFraction))=[];
                ColumnAmountNO2(x)=NaN;             ColumnAmountNO2(y)=NaN;             ColumnAmountNO2(isnan(ColumnAmountNO2))=[];
                SlantColumnAmountNO2(x)=NaN;        SlantColumnAmountNO2(y)=NaN;        SlantColumnAmountNO2(isnan(SlantColumnAmountNO2))=[];
                TerrainHeight(x)=NaN;               TerrainHeight(y)=NaN;               TerrainHeight(isnan(TerrainHeight))=[];
                TerrainPressure(x)=NaN;             TerrainPressure(y)=NaN;             TerrainPressure(isnan(TerrainPressure))=[];
                TerrainReflectivity(x)=NaN;         TerrainReflectivity(y)=NaN;         TerrainReflectivity(isnan(TerrainReflectivity))=[];
                vcdQualityFlags(x)=NaN;             vcdQualityFlags(y)=NaN;             vcdQualityFlags(isnan(vcdQualityFlags))=[];
                RelativeAzimuthAngle(x)=NaN;        RelativeAzimuthAngle(y)=NaN;        RelativeAzimuthAngle(isnan(RelativeAzimuthAngle))=[];
                ColumnAmountNO2Trop(x)=NaN;         ColumnAmountNO2Trop(y)=NaN;         ColumnAmountNO2Trop(isnan(ColumnAmountNO2Trop))=[];
                Row(x)=NaN;                         Row(y)=NaN;                         Row(isnan(Row))=[];
                Swath(x)=NaN;                       Swath(y)=NaN;                       Swath(isnan(Swath))=[];
                
                %Save the imported items to the structure 'Data'.  As is,
                %these structures will be 1 x n, where n is the number of
                %valid pixels.
                Data(E).Latitude = lat(:);
                Data(E).Longitude = lon(:);
                Data(E).Loncorn = loncorn(1:4,:);
                Data(E).Latcorn = latcorn(1:4,:);
                Data(E).SolarAzimuthAngle = SolarAzimuthAngle(:);           Data(E).AMFTrop = AMFTrop(:);
                Data(E).SolarZenithAngle = SolarZenithAngle(:);             Data(E).AMFStrat = AMFStrat(:);
                Data(E).ViewingAzimuthAngle = ViewingAzimuthAngle(:);       Data(E).TerrainHeight = TerrainHeight(:);
                Data(E).ViewingZenithAngle = ViewingZenithAngle(:);         Data(E).TerrainPressure = TerrainPressure(:);
                Data(E).Time = Time(:);                                     Data(E).TerrainReflectivity = TerrainReflectivity(:);
                Data(E).ColumnAmountNO2 = ColumnAmountNO2(:);               Data(E).vcdQualityFlags = vcdQualityFlags(:);
                Data(E).ColumnAmountNO2Trop = ColumnAmountNO2Trop(:);       Data(E).SlantColumnAmountNO2 = SlantColumnAmountNO2(:);
                Data(E).CloudRadianceFraction = CloudRadianceFraction(:);   Data(E).CloudPressure = CloudPressure(:);
                Data(E).RelativeAzimuthAngle = RelativeAzimuthAngle(:);     Data(E).CloudFraction = CloudFraction(:);
                Data(E).Row = Row(:);
                Data(E).Swath = Swath(:);
                
                %Add MODIS cloud info to the files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Convert the OMI date to a Julian calendar day
                d2=1+datenum(str2double(year),str2double(month),str2double(day))-datenum(str2double(year),1,1);
                x=numel(num2str(d2));
                if x==1;
                    day2=(['00',num2str(d2)]);
                elseif x==2;
                    day2=(['0',num2str(d2)]);
                elseif x==3;
                    day2=num2str(d2);
                end
                
                %Find all MODIS files that occur after the current OMI file
                %but before the next OMI file.
                modis_file=(['MYD06_L2.A',year,day2,'*.hdf']);
                modis_files=dir(fullfile(modis_myd06_dir,year,modis_file));
                n=length(modis_files);
                for ii=1:n;
                    mod_filename=modis_files(ii).name;
                    if str2double(mod_filename(19:22))<str2double(sp_files(e).name(29:32));
                        continue
                    elseif str2double(mod_filename(19:22))>str2double(sp_files(e+1).name(29:32));
                        continue
                    else
                        %For each file that fits the criteria mentioned
                        %above, import its latitude, longitude, and cloud
                        %fraction.
                        mod_filename=fullfile(modis_myd06_dir,year,modis_files(ii).name); %Redefine the filename to have the full path to the file
                        mod_fileinfo=hdfinfo(mod_filename);
                        Latitude=hdfread(hdf_dsetID(mod_fileinfo,1,1,'Latitude')); Latitude=double(Latitude); Latitude=Latitude(:);
                        Longitude=hdfread(hdf_dsetID(mod_fileinfo,1,1,'Longitude')); Longitude=double(Longitude); Longitude=Longitude(:);
                        CloudFraction=hdfread(hdf_dsetID(mod_fileinfo,1,2,'Cloud_Fraction')); CloudFraction=double(CloudFraction); CloudFraction=CloudFraction(:);
                        CloudFraction(CloudFraction==127)=100; CloudFraction=CloudFraction*0.009999999776482582;
                        
                        x=find(Longitude>lonmax | Longitude<lonmin);
                        y=find(Latitude>latmax | Latitude<latmin);
                        Longitude(x)=NaN;           Longitude(y)=NaN;               Longitude(isnan(Longitude))=[];
                        Latitude(x)=NaN;            Latitude(y)=NaN;                Latitude(isnan(Latitude))=[];
                        CloudFraction(x)=NaN;       CloudFraction(y)=NaN;           CloudFraction(isnan(CloudFraction))=[];
                        
                        if isempty(Longitude)||isempty(Latitude);
                        else
                            if exist('mod_Data','var')==0;
                                mod_Data(1).Longitude=Longitude;              mod_Data(1).Latitude=Latitude;
                                mod_Data(1).CloudFraction=CloudFraction;
                            elseif exist('mod_Data','var')==1;
                                mod_Data(1).Longitude=[mod_Data(1).Longitude;Longitude];
                                mod_Data(1).Latitude=[mod_Data(1).Latitude;Latitude];
                                mod_Data(1).CloudFraction=[mod_Data(1).CloudFraction;CloudFraction];
                            end
                        end
                    end
                end
                
                %If there is no "mod_Data" variable, fill the regular Data
                %field with -127. Otherwise, find all the MODIS cloud
                %pixels in each OMI pixel and average them together.
                if exist('mod_Data','var')==0;
                    Data(E).MODISCloud=-127*ones(length(Data(E).Latitude),1);
                else
                    for jj=1:length(Data(E).Latitude);
                        x = [];
                        x1 = Data(E).Loncorn(1,jj);   y1 = Data(E).Latcorn(1,jj);
                        x2 = Data(E).Loncorn(2,jj);   y2 = Data(E).Latcorn(2,jj);
                        x3 = Data(E).Loncorn(3,jj);   y3 = Data(E).Latcorn(3,jj);
                        x4 = Data(E).Loncorn(4,jj);   y4 = Data(E).Latcorn(4,jj);
                        
                        xall=[x1;x2;x3;x4;x1];
                        yall=[y1;y2;y3;y4;y1];
                        xx = inpolygon(mod_Data.Latitude,mod_Data.Longitude,yall,xall);
                        
                        cld_vals=mod_Data.CloudFraction(xx);
                        cld_vals(isnan(cld_vals))=[];
                        Data(E).MODISCloud(jj,1)=mean(cld_vals);
                        
                        clear lat1 lat2 lat3 lat4 xx
                    end
                end
                clear mod_Data
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end %End the section carried out only if there are OMI pixels in the area of interest
            
        end %End the loop over all swaths in a day
        savename=[satellite,'_',retrieval,'_',year,month,day];
        save(fullfile(mat_dir,savename), 'Data')
        clear Data 
        toc
        t=toc;
        if t>1200
            disp('Time exceeded 20 min. Stopping')
            quit
        end
    end %End the section checking if there are OMI files for the given time period
end %End the loop over all days
