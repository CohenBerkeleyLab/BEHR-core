% readhe5_omno2_v_aug2012
%Reads omno2 he5 files as of the Aug 2012 version; saves the resulting .mat
%file as <satellite>_<retrieval>_<year><month><day>. Based on
%readhe5_neus_wcld by Ashley Russel.
%
%Josh Laughner <joshlaugh5@gmail.com> 27 Feb 2014

%****************************%
% CONSOLE OUTPUT LEVEL - 0 = none, 1 = minimal, 2 = all messages, 3 = times %
%   4 = save certain debugging variables in the final data structure, including:
%        MODISAlb_Ocean: A matrix that has a value of 1 anywhere the ocean
%        albedo lookup table is used.
% Allows for quick control over the amount of output to the console.
% Choose a higher level to keep track of what the script is doing.
% 3 or less recommended for final products, as 4 will store debugging
% variables in the output file, increasing its size.
DEBUG_LEVEL = 1;
%****************************%

%Add the 'Utils' folder and all subfolders to MATLAB's search path. Within
%the Git repository for BEHR, this is the /Utils folder.
addpath(genpath('/Users/Josh/Documents/MATLAB/BEHR/Utils'))

%Specify the longitude and latitude ranges of interest for this retrieval.
%****************************%
lonmin = -125;    lonmax = -65;
latmin = 25;    latmax = 50;
%****************************%
if lonmin > lonmax %Just in case I enter something backwards...
    error('read_omno2:maxmin','Lonmin is greater than lonmax')
elseif latmin > latmax
    error('read_omno2:maxmin', 'Latmin is greater than latmax')
end

%Process all files between these dates, in yyyy/mm/dd format
%****************************%
date_start='2013/09/01';
date_end='2013/09/30';
%****************************%

%These will be included in the file name
%****************************%
satellite='OMI';
retrieval='SP';
%****************************%

%This will help reject descending node swaths, which occasionally creep in.
%Set to 'US' to implement that feature, set to anything else to disable.
region = 'US';

%This is the directory where the final .mat file will be saved. This will
%need to be changed to match your machine and the files' location. Do not
%include a trailing separator, i.e. '/my/favorite/directory' not
%'my/favorite/directory/'
mat_dir = '/Volumes/share-sat/SAT/BEHR/SP_Files_2014';

%This is the directory where the he5 files are saved. Do not include a
%trailing separator.
he5_dir = '/Volumes/share/GROUP/SAT/OMI/OMNO2_32';

%This is the directory where the MODIS myd06_L2*.hdf files are saved. It should include subfolders organized by year.
%Do not include a trailing separator.
modis_myd06_dir = '/Volumes/share-sat/SAT/MODIS/MYD06_L2';

%This is the directory where the MODIS MCD43C3*.hdf files are saved. It should include subfolders organized by year.
%Do not include a trailing separator.
modis_mcd43_dir = '/Volumes/share/GROUP/SAT/MODIS/MCD43C3';


%This is the directory where the GLOBE data files and their headers (.hdr files) are saved.
%Do not include a trailing separator.
globe_dir = '/Volumes/share/GROUP/SAT/BEHR/GLOBE_files';


%This number will be used

tic %Start the timer

%Initialize matrices to hold the OMI data
Latitude=zeros(60,2000);
Longitude=zeros(60,300);
FoV75CornerLatitude=zeros(300,60,4);
FoV75CornerLongitude=zeros(300,60,4);
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

% [terpres, refvec] = globedem(globe_dir,1,[latmin, latmax],[lonmin, lonmax]);
% %refvec will contain (1) number of cells per degree, (2)
% %northwest corner latitude, (3) NW corner longitude.
% %(2) & (3) might differ from the input latmin & lonmin
% %because of where the globe cell edges fall
% cell_size = refvec(1);
% 
% %GLOBE matrices are arrange s.t. terpres(1,1) is in the SW
% %corner and terpres(end, end) is in the NE corner.
% 
% if DEBUG_LEVEL > 0; fprintf('\n Creating lon/lat matrices for GLOBE data \n'); end
% globe_latmax = refvec(2); globe_latmin = globe_latmax - size(terpres,1)*(1/cell_size);
% globe_lat_matrix = (globe_latmin + 1/(2*cell_size)):(1/cell_size):globe_latmax;
% globe_lat_matrix = globe_lat_matrix';
% globe_lat_matrix = repmat(globe_lat_matrix,1,size(terpres,2));
% 
% globe_lonmin = refvec(3); globe_lonmax = globe_lonmin + size(terpres,2)*(1/cell_size);
% globe_lon_matrix = globe_lonmin + 1/(2*cell_size):(1/cell_size):globe_lonmax;
% globe_lon_matrix = repmat(globe_lon_matrix,size(terpres,1),1);

%File names will be prefixed with "<satellite>_<retrieval>_", e.g. for OMI
%satellite SP retrieval, the prefix will be "OMI_SP_" and then the date in
%year, month, date order.  This section checks to see if the last file in
%the mat directory has the expected prefix.  If so, that date is taken as
%the last date completed, otherwise it is assumed that the retrieval will
%need to start from the specified start date. This allows he5 reading to be
%stopped and restarted with minimal intervention.
file_prefix = [satellite,'_',retrieval,'_']; l = length(file_prefix);
last_file=dir(fullfile(mat_dir,[file_prefix,'*.mat']));

if ~isempty(last_file)
    last_datenum = datenum(last_file(end).name(l+1:l+8),'yyyymmdd')+1;
else
    last_datenum = 0;
end

if last_datenum >= datenum(date_start) && last_datenum <= datenum(date_end)
    datenums = last_datenum:datenum(date_end);
else
    datenums = datenum(date_start):datenum(date_end);
end

% if isempty(last_file) || ~strcmp(last_file(end).name(1:l),file_prefix);
%     last_date=datestr(datenum(date_start)-1,26);
% else
%     last_date=last_file(end).name((l+1):(l+8));
%     last_date=([last_date(1:4),'/',last_date(5:6),'/',last_date(7:8)]);
% end

%Go ahead and load the terrain pressure data - only need to do this once
[terpres, refvec] = globedem(globe_dir,1,[latmin, latmax],[lonmin, lonmax]);
    %refvec will contain (1) number of cells per degree, (2)
    %northwest corner latitude, (3) NW corner longitude.
    %(2) & (3) might differ from the input latmin & lonmin
    %because of where the globe cell edges fall
if DEBUG_LEVEL > 0; fprintf('\n Creating lon/lat matrices for GLOBE data \n'); end
cell_count = refvec(1);
globe_latmax = refvec(2); globe_latmin = globe_latmax - size(terpres,1)*(1/cell_count);
globe_lat_matrix = (globe_latmin + 1/(2*cell_count)):(1/cell_count):globe_latmax;
globe_lat_matrix = globe_lat_matrix';
globe_lat_matrix = repmat(globe_lat_matrix,1,size(terpres,2));

globe_lonmin = refvec(3); globe_lonmax = globe_lonmin + size(terpres,2)*(1/cell_count);
globe_lon_matrix = globe_lonmin + 1/(2*cell_count):(1/cell_count):globe_lonmax;
globe_lon_matrix = repmat(globe_lon_matrix,size(terpres,1),1); 

terpres(isnan(terpres)) = -500;

%For loop over all days from the starting or last finished date to the end
%date. We will give the absolute paths to files rather than changing the
%active directory, as MATLAB seems to run slightly slower if the current
%working directory is on the server.
% total_days=datenum(date_end)-datenum(last_date)+1;
% for j=1:total_days;
for j=1:length(datenums)
    %Read the desired year, month, and day
    R=datenums(j);
    date=datestr(R,26);
    year=date(1:4);
    month=date(6:7);
    day=date(9:10);
    
    %Prepare a data structure to receive the final data.
    Data=struct('Date',0,'Longitude',0,'Latitude',0,'LatBdy',[],'LonBdy',[],'Loncorn',0,'Latcorn',0,'Time',0,'ViewingZenithAngle',0,'SolarZenithAngle',0,'ViewingAzimuthAngle',0,'SolarAzimuthAngle',0,'AMFStrat',0,'AMFTrop',0,'CloudFraction',0,'CloudRadianceFraction',0,'TerrainHeight',0,'TerrainPressure',0,'TerrainReflectivity',0,'vcdQualityFlags',0,'CloudPressure',0,'ColumnAmountNO2',0,'SlantColumnAmountNO2',0,'ColumnAmountNO2Trop',0,'MODISCloud',0,'MODIS_Cloud_File','','MODISAlbedo',0,'MODIS_Albedo_File','','GLOBETerpres',0,'XTrackQualityFlags',0);
    
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
            if DEBUG_LEVEL > 0
                if e==1 || mod(e,10)==0; fprintf('Swath %u of %s/%s/%s \n',e,month,day,year); end
            end
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
            
            
            time_ind = regexp(filename,'t\d\d\d\d-o');
            omi_starttime = str2double(filename(time_ind+1:time_ind+4));
            
            if isempty(Lon)==1 || isempty(Lat)==1 || length(Lat)==1;
                if DEBUG_LEVEL > 1; disp('No points within lat/lon boundaries'); end
                continue
            elseif omi_starttime < 1500 && strcmp(region,'US')
                %If start time is < 1500 and we want to look at the US, reject
                %the file, as it is probably descending nodes only.
                if DEBUG_LEVEL > 0; fprintf(' Swath %d: Nighttime granule skipped\n',e); end
                continue
            else
                if DEBUG_LEVEL > 1; disp('Founds points within lat/lon boundaries'); end
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
                
                if DEBUG_LEVEL > 0; fprintf('\n Importing OMI data fields \n'); end
                
                %Import the FoV75 corner lat and lons.  These will be
                %ordered the same as the BEHR-calculated corners, i.e. corner x
                %along track x across track
                slabsize = [length(cut_y),60,4];
                memspaceID = H5S.create_simple(length(slabsize),slabsize,slabsize);
                offset = [(min(i_i)-1),0,0];
                
                %Occasionally there are problems where the corner lat/lon
                %fields in the OMI files don't have 4 points.  If that is
                %the case, fill those fields with NaNs.  If some other
                %error occurs, rethrow it.
                try
                    datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'FoV75CornerLatitude')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); FoV75CornerLatitude = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); FoV75CornerLatitude = double(FoV75CornerLatitude); FoV75CornerLatitude = permute(FoV75CornerLatitude, [1 3 2]);
                    datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,2,'FoV75CornerLongitude')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); FoV75CornerLongitude = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); FoV75CornerLongitude = double(FoV75CornerLongitude); FoV75CornerLongitude = permute(FoV75CornerLongitude, [1 3 2]);
                catch err
                    if strcmp(err.identifier,'MATLAB:imagesci:hdf5lib:libraryError') 
                        FoV75CornerLatitude = nan(4,length(cut_y),60);
                        FoV75CornerLongitude = nan(4,length(cut_y),60);
                    else
                        rethrow(err);
                    end
                end
                    
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
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'CloudFraction')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudFraction = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudFraction=double(CloudFraction); CloudFraction=CloudFraction'; CloudFraction = CloudFraction/1000;
                %CloudFractionError
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'CloudFractionStd')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudFractionError = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudFractionError=double(CloudFractionError); CloudFractionError=CloudFractionError';
                %CloudPressure
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'CloudPressure')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudPressure = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudPressure=double(CloudPressure); CloudPressure=CloudPressure';
                %CloudPressureError
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'CloudPressureStd')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudPressureError = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudPressureError=double(CloudPressureError); CloudPressureError=CloudPressureError';
                %CloudRadianceFraction
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1, 'CloudRadianceFraction')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); CloudRadianceFraction = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); CloudRadianceFraction=double(CloudRadianceFraction); CloudRadianceFraction=CloudRadianceFraction'; CloudRadianceFraction = CloudRadianceFraction/1000;
                %ColumnAmountNO2
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'ColumnAmountNO2')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ColumnAmountNO2 = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ColumnAmountNO2=double(ColumnAmountNO2); ColumnAmountNO2=ColumnAmountNO2';
                %ColumnAmountNO2Trop
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'ColumnAmountNO2Trop')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ColumnAmountNO2Trop = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ColumnAmountNO2Trop=double(ColumnAmountNO2Trop); ColumnAmountNO2Trop=ColumnAmountNO2Trop';
                %SlantColumnAmountNO2
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'SlantColumnAmountNO2')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SlantColumnAmountNO2 = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SlantColumnAmountNO2=double(SlantColumnAmountNO2); SlantColumnAmountNO2=SlantColumnAmountNO2';
                %TerrainHeight
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'TerrainHeight')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); TerrainHeight = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); TerrainHeight=double(TerrainHeight); TerrainHeight=TerrainHeight';
                %TerrainPressure
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'TerrainPressure')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); TerrainPressure = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); TerrainPressure=double(TerrainPressure); TerrainPressure=TerrainPressure';
                %TerrainReflectivity
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'TerrainReflectivity')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); TerrainReflectivity = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); TerrainReflectivity=double(TerrainReflectivity); TerrainReflectivity=TerrainReflectivity';
                %TropopausePressure
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'TropopausePressure')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); TropopausePressure = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); TropopausePressure=double(TropopausePressure); TropopausePressure=TropopausePressure';
                %vcdQualityFlags
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'VcdQualityFlags')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); vcdQualityFlags = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); vcdQualityFlags=double(vcdQualityFlags); vcdQualityFlags=vcdQualityFlags';
                %XTrackQualityFlags
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'XTrackQualityFlags')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); XTrackQualityFlags = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); XTrackQualityFlags=double(XTrackQualityFlags); XTrackQualityFlags=XTrackQualityFlags';
                
                H5F.close(fileID); %close omi file to free up space
                
                RelativeAzimuthAngle=abs(SolarAzimuthAngle+180-ViewingAzimuthAngle);
                RelativeAzimuthAngle(RelativeAzimuthAngle > 180)=360-RelativeAzimuthAngle(RelativeAzimuthAngle > 180);
                
                % Use Inf instead of NaN here because some days seems to
                % have NaN values populated in NO2/NO2 Trop column density,
                % with the result that this section causes unequal removal
                % of matrix elements if NaNs used as markers.
                %x=find(lon<lonmin | lon>lonmax);
                %y=find(lat<latmin | lat>latmax);
                x = find((lon < lonmin | lon > lonmax) | (lat < latmin | lat > latmax));
                lon(x)=[];                         %lon(y)=Inf;                         lon(isinf(lon))=[];
                lat(x)=[];                         %lat(y)=Inf;                         lat(isinf(lat))=[];
                loncorn(:,x)=[];                   %loncorn(:,y)=Inf;                   loncorn(:,any(isinf(loncorn))) = [];
                latcorn(:,x)=[];                   %latcorn(:,y)=Inf;                   latcorn(:,any(isinf(latcorn))) = [];
                FoV75CornerLatitude(:,x)=[];       %FoV75CornerLatitude(:,y)=Inf;       FoV75CornerLatitude(:,any(isinf(FoV75CornerLatitude))) = [];
                FoV75CornerLongitude(:,x)=[];      %FoV75CornerLongitude(:,y)=Inf;      FoV75CornerLongitude(:,any(isinf(FoV75CornerLongitude))) = [];
                SolarAzimuthAngle(x)=[];           %SolarAzimuthAngle(y)=Inf;           SolarAzimuthAngle(isinf(SolarAzimuthAngle))=[];
                SolarZenithAngle(x)=[];            %SolarZenithAngle(y)=Inf;            SolarZenithAngle(isinf(SolarZenithAngle))=[];
                ViewingAzimuthAngle(x)=[];         %ViewingAzimuthAngle(y)=Inf;         ViewingAzimuthAngle(isinf(ViewingAzimuthAngle))=[];
                ViewingZenithAngle(x)=[];          %ViewingZenithAngle(y)=Inf;          ViewingZenithAngle(isinf(ViewingZenithAngle))=[];
                Time(x)=[];                        %Time(y)=Inf;                        Time(isinf(Time))=[];
                AMFStrat(x)=[];                    %AMFStrat(y)=Inf;                    AMFStrat(isinf(AMFStrat))=[];
                AMFTrop(x)=[];                     %AMFTrop(y)=Inf;                     AMFTrop(isinf(AMFTrop))=[];
                CloudFraction(x)=[];               %CloudFraction(y)=Inf;               CloudFraction(isinf(CloudFraction))=[];
                CloudPressure(x)=[];               %CloudPressure(y)=Inf;               CloudPressure(isinf(CloudPressure))=[];
                CloudRadianceFraction(x)=[];       %CloudRadianceFraction(y)=Inf;       CloudRadianceFraction(isinf(CloudRadianceFraction))=[];
                ColumnAmountNO2(x)=[];             %ColumnAmountNO2(y)=Inf;             ColumnAmountNO2(isinf(ColumnAmountNO2))=[];
                SlantColumnAmountNO2(x)=[];        %SlantColumnAmountNO2(y)=Inf;        SlantColumnAmountNO2(isinf(SlantColumnAmountNO2))=[];
                TerrainHeight(x)=[];               %TerrainHeight(y)=Inf;               TerrainHeight(isinf(TerrainHeight))=[];
                TerrainPressure(x)=[];             %TerrainPressure(y)=Inf;             TerrainPressure(isinf(TerrainPressure))=[];
                TerrainReflectivity(x)=[];         %TerrainReflectivity(y)=Inf;         TerrainReflectivity(isinf(TerrainReflectivity))=[];
                TropopausePressure(x) = [];
                vcdQualityFlags(x)=[];             %vcdQualityFlags(y)=Inf;             vcdQualityFlags(isinf(vcdQualityFlags))=[];
                XTrackQualityFlags(x)=[];          %XTrackQualityFlags(y)=Inf;          XTrackQualityFlags(isinf(XTrackQualityFlags))=[];
                RelativeAzimuthAngle(x)=[];        %RelativeAzimuthAngle(y)=Inf;        RelativeAzimuthAngle(isinf(RelativeAzimuthAngle))=[];
                ColumnAmountNO2Trop(x)=[];         %ColumnAmountNO2Trop(y)=Inf;         ColumnAmountNO2Trop(isinf(ColumnAmountNO2Trop))=[];
                Row(x)=[];                         %Row(y)=Inf;                         Row(isinf(Row))=[];
                Swath(x)=[];                       %Swath(y)=Inf;                       Swath(isinf(Swath))=[];
                
                if DEBUG_LEVEL > 0; disp(' Saving imported OMI fields to "Data"'); end
                %Save the imported items to the structure 'Data'.  As is,
                %these structures will be vectors.
                Data(E).Latitude = lat(:);                                  Data(E).LatBdy = [latmin latmax];
                Data(E).Longitude = lon(:);                                 Data(E).LonBdy = [lonmin lonmax];
                Data(E).Loncorn = loncorn(1:4,:);                           Data(E).FoV75CornerLongitude = FoV75CornerLongitude(1:4,:);
                Data(E).Latcorn = latcorn(1:4,:);                           Data(E).FoV75CornerLatitude = FoV75CornerLatitude(1:4,:);
                Data(E).SolarAzimuthAngle = SolarAzimuthAngle(:);           Data(E).AMFTrop = AMFTrop(:);
                Data(E).SolarZenithAngle = SolarZenithAngle(:);             Data(E).AMFStrat = AMFStrat(:);
                Data(E).ViewingAzimuthAngle = ViewingAzimuthAngle(:);       Data(E).TerrainHeight = TerrainHeight(:);
                Data(E).ViewingZenithAngle = ViewingZenithAngle(:);         Data(E).TerrainPressure = TerrainPressure(:);
                Data(E).Time = Time(:);                                     Data(E).TerrainReflectivity = TerrainReflectivity(:);
                Data(E).ColumnAmountNO2 = ColumnAmountNO2(:);               Data(E).vcdQualityFlags = vcdQualityFlags(:);
                Data(E).ColumnAmountNO2Trop = ColumnAmountNO2Trop(:);       Data(E).SlantColumnAmountNO2 = SlantColumnAmountNO2(:);
                Data(E).CloudRadianceFraction = CloudRadianceFraction(:);   Data(E).CloudPressure = CloudPressure(:);
                Data(E).RelativeAzimuthAngle = RelativeAzimuthAngle(:);     Data(E).CloudFraction = CloudFraction(:);
                Data(E).Row = Row(:);                                       Data(E).XTrackQualityFlags = XTrackQualityFlags(:);
                Data(E).Swath = Swath(:);                                   Data(E).Date=date;
                Data(E).TropopausePressure = TropopausePressure(:);
                
                %Add MODIS cloud info to the files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if DEBUG_LEVEL > 0; fprintf('\n Adding MODIS cloud data \n'); end
                
                %Convert the OMI date to a Julian calendar day
                d2=1+datenum(str2double(year),str2double(month),str2double(day))-datenum(str2double(year),1,1);
                x=numel(num2str(d2));
                if x==1;
                    julian_day=(['00',num2str(d2)]);
                elseif x==2;
                    julian_day=(['0',num2str(d2)]);
                elseif x==3;
                    julian_day=num2str(d2);
                end
                
                %Find all MODIS files that occur after the current OMI file
                %but before the next OMI file.
                modis_file=(['MYD06_L2.A',year,julian_day,'*.hdf']);
                modis_files=dir(fullfile(modis_myd06_dir,year,modis_file));
                
                    % Calculate the start time for the next OMI swath.
                    % Usually there is at least one swath starting after
                    % the one overflying the US west coast for that day,
                    % but occasionally that swath is not present in the
                    % OMNO2 data.  In those cases, we need to calculate the
                    % time it should have started, knowing that the Aura
                    % orbit period is ~99 min.  If for some reason the
                    % calculated start time should end up being in the next
                    % day, error out so the user is aware of that.
                if e < n
                    next_omi_swath_time = str2double(sp_files(e+1).name(29:32));
                else
                    omi_hr = str2double(sp_files(e).name(29:30));
                    omi_min = str2double(sp_files(e).name(31:32));
                    next_omi_swath_time = 100*(floor((omi_min + 99)/60)+omi_hr) + mod(omi_min + 99,60);
                    if next_omi_swath_time > 2359; error('read_omno2:modis_cloud','Next OMI swath time for MODIS cloud binning calculated to be in the next day. \nManual attention recommended'); end
                end
                
                for ii=1:length(modis_files);
                    mod_filename=modis_files(ii).name;
                    if str2double(mod_filename(19:22))<str2double(sp_files(e).name(29:32));
                        continue
                    elseif e < n && str2double(mod_filename(19:22))>next_omi_swath_time;
                        continue
                    else
                        %For each file that fits the criteria mentioned
                        %above, import its latitude, longitude, and cloud
                        %fraction.
                        if DEBUG_LEVEL > 0; fprintf('  Averaging MODIS cloud file %s\n',mod_filename); end
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
                        Data(E).MODIS_Cloud_File=mod_filename;
                        
                        clear xx
                    end
                end
                clear mod_Data
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Add MODIS albedo info to the files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if DEBUG_LEVEL>0; fprintf('\n Adding MODIS albedo information \n'); end
                alb_dir = fullfile(modis_mcd43_dir,year);
                
                %Find the closest MCD file
                in=[0 1 -1 2 -2 3 -3 4 -4 5 -5 6 -6 7 -7 8 -8 9 -9 10 -10 11 -11 12 -12 13 -13 14 -14 15 -15 16 -16 17 -17 18 -18 19 -19 20 -20 21 -21];
                for ii=1:length(in);
                    mcd_date = num2str(str2double(julian_day) + in(ii),'%03g');
                    alb_filename = fullfile(alb_dir,['MCD43C3.A',year,mcd_date,'*.hdf']);
                    alb_files = dir(alb_filename);
                    if DEBUG_LEVEL > 1; fprintf('Looking for %s \n', alb_filename); end
                    %if exist('alb_filename','file') == 2
                    if ~isempty(alb_files)
                        if DEBUG_LEVEL > 0; fprintf(' Found mcd43 file %s \n', alb_filename); end
                        break
                    elseif ii==length(in)
                        error('read_omno2:mcd43c3_find','Could not find an MCD43C3 (Albedo) file within 21 days');
                    end
                end
                
                mcd43_info = hdfinfo(fullfile(alb_dir,alb_files(1).name));
                band3 = hdfread(hdf_dsetID(mcd43_info,1,1,'Albedo_BSA_Band3'));
                band3 = double(band3);
                band3 = flipud(band3); 
                
                %MODIS albedo is given in 0.05 degree cells and a single file covers the
                %full globe, so figure out the lat/lon of the middle of the grid cells as:
                band3_lat=-90+0.05/2:0.05:90-0.05/2; band3_lats=band3_lat'; band3_lats=repmat(band3_lats,1,7200);
                band3_lon=-180+0.05/2:0.05:180-0.05/2; band3_lons=repmat(band3_lon,3600,1);
                
                %To speed up processing, restrict the MODIS albedo data to
                %only the area we need to worry about.  This will
                %significantly speed up the search for albedo values within
                %each pixel.
                lat_min=Data(E).Latcorn(:); lat_min(lat_min==0)=[]; lat_min=floor(min(lat_min));
                lat_max=Data(E).Latcorn(:); lat_max(lat_max==0)=[]; lat_max=ceil(max(lat_max));
                lon_min=Data(E).Loncorn(:); lon_min(lon_min==0)=[]; lon_min=floor(min(lon_min));
                lon_max=Data(E).Loncorn(:); lon_max(lon_max==0)=[]; lon_max=ceil(max(lon_max));
                
                in_lats = find(band3_lat>=lat_min & band3_lat<=lat_max);
                in_lons = find(band3_lon>=lon_min & band3_lon<=lon_max);
                band3=band3(in_lats,in_lons); 
                band3(band3==32767)=NaN; %JLL 11 Apr 2014: 32767 is the fill value for this data set; we will remove the NaNs further down
                band3 = band3 * 1e-3; %JLL 11-Apr-2014 Albedo matrix needs to have the scale factor applied 
                band3_lats=band3_lats(in_lats,in_lons);
                band3_lons=band3_lons(in_lats,in_lons);
                s=size(Data(E).Latitude);
                c=numel(Data(E).Latitude);
                MODISAlbedo=zeros(s);
                
                if DEBUG_LEVEL > 3; MODISAlb_Ocn = zeros(s); end %JLL 
                
                %Now actually average the MODIS albedo for each OMI pixel
                if DEBUG_LEVEL > 0; disp(' Averaging MODIS albedo to OMI pixels'); end
                for k=1:c;
                    if DEBUG_LEVEL > 2; tic; end
                    x = [];
                    x1 = Data(E).Loncorn(1,k);   y1 = Data(E).Latcorn(1,k);
                    x2 = Data(E).Loncorn(2,k);   y2 = Data(E).Latcorn(2,k);
                    x3 = Data(E).Loncorn(3,k);   y3 = Data(E).Latcorn(3,k);
                    x4 = Data(E).Loncorn(4,k);   y4 = Data(E).Latcorn(4,k);
                    
                    
                    xall=[x1;x2;x3;x4;x1];
                    yall=[y1;y2;y3;y4;y1];
                    
                    xx = inpolygon(band3_lats,band3_lons,yall,xall);
                    
                    band3_vals=band3(xx);  band3_zeros=find(band3_vals==0);
                    band3_vals(band3_zeros)=NaN; band3_vals(isnan(band3_vals))=[];
                    band3_avg=mean(band3_vals);
                    
                    %put in ocean surface albedo from LUT
                    if isnan(band3_avg)==1;
                        sza_vec = [5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 89];
                        alb_vec = [0.038 0.038 0.039 0.039 0.040 0.042 0.044 0.046 0.051 0.058 0.068 0.082 0.101 0.125 0.149 0.158 0.123 0.073];
                        alb = interp1(sza_vec,alb_vec,Data(E).SolarZenithAngle(k));
                        band3_avg = alb;
                        if DEBUG_LEVEL > 3; MODISAlb_Ocn(k) = 1; end
                    end
                    
                    MODISAlbedo(k) = band3_avg;
                    if DEBUG_LEVEL > 2; telap = toc; fprintf(' Time for MODIS alb --> pixel %u/%u = %g sec \n',k,c,telap); end
                end
                
                Data(E).MODISAlbedo = MODISAlbedo;
                Data(E).MODIS_Albedo_File = fullfile(alb_dir,alb_files(1).name);
                if DEBUG_LEVEL > 3; Data(E).MODISAlb_Ocean = MODISAlb_Ocn; end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %Add GLOBE terrain pressure to the files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if DEBUG_LEVEL > 0; fprintf('\n Adding GLOBE terrain data \n'); end
                
                GLOBETerpres = zeros(size(Data(E).Latitude));
                %GLOBETerHeight = zeros(size(Data(E).Latitude));
                
                %GLOBE matrices are arrange s.t. terpres(1,1) is in the SW
                %corner and terpres(end, end) is in the NE corner.

                for k=1:c
                    
                    if DEBUG_LEVEL > 1; fprintf('Averaging GLOBE data to pixel %u of %u \n',k,c); end
                    if DEBUG_LEVEL > 2; tic; end
                    x1 = Data(E).Loncorn(1,k);   y1 = Data(E).Latcorn(1,k);
                    x2 = Data(E).Loncorn(2,k);   y2 = Data(E).Latcorn(2,k);
                    x3 = Data(E).Loncorn(3,k);   y3 = Data(E).Latcorn(3,k);
                    x4 = Data(E).Loncorn(4,k);   y4 = Data(E).Latcorn(4,k);
                    
                    
                    xall=[x1;x2;x3;x4;x1];
                    yall=[y1;y2;y3;y4;y1];
                    
                    %%%%SPEED IT UP%%%%
                    ai=find(globe_lat_matrix(:,1)>=min(yall) & globe_lat_matrix(:,1)<=max(yall));
                    bi=find(globe_lon_matrix(1,:)>=min(xall) & globe_lon_matrix(1,:)<=max(xall));
                    pressurex=terpres(ai,bi);
                    pressure_latx=globe_lat_matrix(ai,bi);
                    pressure_lonx=globe_lon_matrix(ai,bi);
                    %%%%%%%%%%%%%%%%%%%
                    
                    xx = inpolygon(pressure_latx,pressure_lonx,yall,xall);
                    
                    pres_vals=pressurex(xx);  pres_zeros=find(pres_vals==0);
                    pres_vals(pres_zeros)=NaN; pres_vals(isnan(pres_vals))=[];
                    %GLOBETerHeight(k) = mean(pres_vals);
                    GLOBETerpres(k)=1013.25 .* exp(-mean(pres_vals) / 7400 ); %Originally divided by 7640 m
                    
                    
                    if DEBUG_LEVEL > 2; telap = toc; fprintf('Time for GLOBE --> pixel %u/%u = %g sec \n',k,c,telap); end
                end
                
                Data(E).GLOBETerpres = GLOBETerpres;
                %Data(E).GLOBEHeight = GLOBETerHeight;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end %End the section carried out only if there are OMI pixels in the area of interest
        
        end %End the loop over all swaths in a day
        savename=[satellite,'_',retrieval,'_',year,month,day];
        save(fullfile(mat_dir,savename), 'Data')
        clear Data
        toc
        t=toc;
        if t>1200
            %error('Time exceeded 20 min. Stopping')
        end
    end %End the section checking if there are OMI files for the given time period
end %End the loop over all days
