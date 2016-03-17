function read_omno2_v_aug2012(date_start, date_end)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PARALLELIZATION OPTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specifies whether the script is executing on a cluster; this must be set
% (globally) in the calling script.  This allows for the execution of code
% needed on the cluster (i.e. adding necessary folders to the Matlab path,
% opening a parallel pool) without running them on the local machine.  If
% onCluster hasn't been defined yet, set it to false.
global onCluster;
if isempty(onCluster); 
    if DEBUG_LEVEL > 0; fprintf('Assuming onCluster is false\n'); end
    onCluster = false; 
end

% Defined the number of threads to run, this will be used to open a
% parallel pool. numThreads should be set in the calling run script,
% otherwise it will default to 1.
global numThreads;
if isempty(numThreads)
    numThreads = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DEPENDENCIES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%Add the 'Utils' folder and all subfolders to MATLAB's search path. Within
%the Git repository for BEHR, this is the /Utils folder.
addpath(genpath('~/Documents/MATLAB/BEHR/Utils'))


% Add the paths needed to run on the cluster
if onCluster;
    addpath(genpath('~/MATLAB/Classes'));
    addpath(genpath('~/MATLAB/Utils'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%
%%%%% USER INPUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%

%Specify the longitude and latitude ranges of interest for this retrieval.
%****************************%
lonmin = -125;    lonmax = -65;
latmin = 25;    latmax = 50;
%****************************%
if lonmin > lonmax %Just in case I enter something backwards...
    error(E.badinput('Lonmin is greater than lonmax'))
elseif latmin > latmax
    error(E.badinput('Latmin is greater than latmax'))
end

%Process all files between these dates, in yyyy/mm/dd format unless
%overriding dates are passed into the function.
%****************************%
if nargin < 2
    date_start='2013/08/01';
    date_end='2013/08/06';
end

% Set to true to resume from last completed file,
% set to false to overwrite any existing files.
restart = false;
%****************************%

%These will be included in the file name
%****************************%
satellite='OMI';
retrieval='SP';
%****************************%

%This will help reject descending node swaths, which occasionally creep in.
%Set to 'US' to implement that feature, set to anything else to disable.
region = 'US';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DATA DIRECTORIES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The location of the directories to read or save data to.  If onCluster is
% true, these will need to be set in the runscript - I figured this would
% be easier than setting them as shell environmental variables and using
% getenv - JLL 15 Jan 2015

if onCluster
    global sp_mat_dir;
    global omi_he5_dir;
    global modis_myd06_dir;
    global modis_mcd43_dir;
    global globe_dir;
    
    % Verify the paths integrity.
    nonexistant = {};
    
    if ~exist(sp_mat_dir,'dir')
        nonexistant{end+1} = 'sp_mat_dir';
    end
    if ~exist(omi_he5_dir,'dir')
        nonexistant{end+1} = 'he5_dir';
    end
    if ~exist(modis_myd06_dir,'dir')
        nonexistant{end+1} = 'modis_myd06_dir';
    end
    if ~exist(modis_mcd43_dir,'dir')
        nonexistant{end+1} = 'modis_mcd43_dir';
    end
    if ~exist(globe_dir,'dir')
        nonexistant{end+1} = 'globe_dir';
    end
    
    if numel(nonexistant)>0
        string_spec = [repmat('\n\t%s',1,numel(nonexistant)),'\n\n'];
        msg = sprintf('The following paths are not valid: %s Please double check them in the run file',string_spec);
        error(E.callError('bad_cluster_path',sprintf(msg,nonexistant{:})));
    end
    
    
else
    %This is the directory where the final .mat file will be saved. This will
    %need to be changed to match your machine and the files' location.
    sp_mat_dir = '/Volumes/share-sat/SAT/BEHR/SP_Files_2014';
    
    %This is the directory where the he5 files are saved.
    omi_he5_dir = '/Volumes/share-sat/SAT/OMI/OMNO2';
    
    %This is the directory where the MODIS myd06_L2*.hdf files are saved. It should include subfolders organized by year.
    modis_myd06_dir = '/Volumes/share-sat/SAT/MODIS/MYD06_L2';
    
    %This is the directory where the MODIS MCD43C3*.hdf files are saved. It should include subfolders organized by year.
    modis_mcd43_dir = '/Volumes/share-sat/SAT/MODIS/MCD43C3';
    
    %This is the directory where the GLOBE data files and their headers (.hdr files) are saved.
    %Do not include a trailing separator.
    %globe_dir = '/Volumes/share/GROUP/SAT/BEHR/GLOBE_files';
    globe_dir = globepath; % a function that returns the GLOBE directory path
end



%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN BODY %%%%%
%%%%%%%%%%%%%%%%%%%%%

%This number will be used

tic %Start the timer

%Initialize matrices to hold the OMI data
% Latitude=zeros(60,2000);
% Longitude=zeros(60,300);
% FoV75CornerLatitude=zeros(300,60,4);
% FoV75CornerLongitude=zeros(300,60,4);
% SpacecraftAltitude=zeros(300,1);
% SpacecraftLatitude=zeros(300,1);
% SpacecraftLongitude=zeros(300,1);
% Time=zeros(300,1);
% ViewingZenithAngle=zeros(60,300);
% SolarZenithAngle=zeros(60,300);
% ViewingAzimuthAngle=zeros(60,300);
% SolarAzimuthAngle=zeros(60,300);
% AMFStrat=zeros(60,300);
% AMFTrop=zeros(60,300);
% CloudFraction=zeros(60,300);
% CloudRadianceFraction=zeros(60,300);
% ColumnAmountNO2=zeros(60,300);
% SlantColumnAmountNO2=zeros(60,300);
% TerrainHeight=zeros(60,300);
% TerrainPressure=zeros(60,300);
% TerrainReflectivity=zeros(60,300);
% vcdQualityFlags=zeros(60,300);
% CloudPressure=zeros(60,300);
% ColumnAmountNO2Trop=zeros(60,300);

%File names will be prefixed with "<satellite>_<retrieval>_", e.g. for OMI
%satellite SP retrieval, the prefix will be "OMI_SP_" and then the date in
%year, month, date order.  This section checks to see if the last file in
%the mat directory has the expected prefix.  If so, that date is taken as
%the last date completed, otherwise it is assumed that the retrieval will
%need to start from the specified start date. This allows he5 reading to be
%stopped and restarted with minimal intervention.
file_prefix = [satellite,'_',retrieval,'_']; l = length(file_prefix);
last_file=dir(fullfile(sp_mat_dir,[file_prefix,'*.mat']));

if ~isempty(last_file) && restart == true
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

if onCluster && isempty(gcp('nocreate'))
    parpool(numThreads);
end

parfor j=1:length(datenums)
    %Read the desired year, month, and day
    R=datenums(j);
    date=datestr(R,26);
    year=date(1:4);
    month=date(6:7);
    day=date(9:10);
    
    %Prepare a data structure to receive the final data. 
    % As I understand this currently (15 Jan 2015) "Data" should be a local
    % variable on each worker for parallelization, and this line will reset
    % its value within each loop. Thus we shouldn't need to worry about
    % cross-communication between workers.
    Data=struct('Date',0,'Longitude',0,'Latitude',0,'LatBdy',[],'LonBdy',[],'Loncorn',0,'Latcorn',0,'Time',0,'ViewingZenithAngle',0,'SolarZenithAngle',0,'ViewingAzimuthAngle',0,'SolarAzimuthAngle',0,'AMFStrat',0,'AMFTrop',0,'CloudFraction',0,'CloudRadianceFraction',0,'TerrainHeight',0,'TerrainPressure',0,'TerrainReflectivity',0,'vcdQualityFlags',0,'CloudPressure',0,'ColumnAmountNO2',0,'SlantColumnAmountNO2',0,'ColumnAmountNO2Trop',0,'MODISCloud',0,'MODIS_Cloud_File','','MODISAlbedo',0,'MODIS_Albedo_File','','GLOBETerpres',0,'XTrackQualityFlags',0);
    
    %Set the file path and name, assuming that the file structure is
    %<he5_directory>/<year>/<month>/...files...  Then figure out how many
    %files there are
    short_filename=['OMI-Aura_L2-OMNO2_',year,'m',month,day,'*.he5'];
    file_dir = fullfile(omi_he5_dir,year,month); %Used both here to find all he5 files and in the swath for loop to identify each file.
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
            % Row will keep track of the pixel's location in the
            % across-track direction; Pixel will track it's position in the
            % along-track direction.  These will allow the spatial
            % arrangement of the pixels to be reconstructed if desired (as
            % a sparse matrix); Row is also used for row anomaly rejection.
            Row=0:59; Row=Row'; Row=repmat(Row,1,size(Latitude,2));
            Pixel=0:length(Latitude)-1; Pixel=repmat(Pixel,60,1);
            Swath=filename(35:39); Swath=str2double(Swath).*ones(size(Latitude));
            
            %Restrict latitude to those that fall within the bounds specified
            %at the begininning of the file. Also pivot the dataset so that
            %the matrix is along track x across track.
            lat=Latitude';
            lat_i=[latmin, latmax];
            [i_i, j_j]=find(lat > lat_i(1) - 0.25 & lat < lat_i(2) + 0.25);
            cut_y=min(i_i):max(i_i);
            cut_x = 1:60;
            lat=double(lat(cut_y,cut_x));
            Latitude=Latitude(cut_x,cut_y)'; Latitude=double(Latitude);
            Row=Row(cut_x,cut_y)';
            Pixel=Pixel(cut_x,cut_y)';
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
            %lines of pixels along track).  It also converts all data from
            %single precision to double precision and pivots the matrix the
            %the convention of row = swath. These are the values needed to
            %compute the pixel corner points.
            
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
            
            
            %This will handle values that only have a single value per
            %across track line of pixels. They are still converted to
            %double precision numbers and pivoted.
            offset = [(min(i_i)-1)]; % need to change to 0-based indexing for HDF files
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
            Time=repmat(Time',1,60);
            
            
            
            %Deletes any rows that have no points in the domain specified
            %Lon and Lat will be used
            Lat=Latitude; Lon=Longitude;
            x = Lon < lonmax & Lon > lonmin;
            y = Lat < latmax & Lat > latmin;
            rows_to_keep = any(x,2) & any(y,2);
            Lon = Lon(rows_to_keep,:);
            Lat = Lat(rows_to_keep,:);
            
            
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
                latcorn = permute(latcorn,[3,1,2]);
                loncorn = corners(:,:,1,1:4); loncorn = squeeze(loncorn);
                loncorn = permute(loncorn,[3,1,2]);
                
                
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
                %ColumnAmountNO2Strat
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'ColumnAmountNO2Strat')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ColumnAmountNO2Strat = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ColumnAmountNO2Strat=double(ColumnAmountNO2Strat); ColumnAmountNO2Strat=ColumnAmountNO2Strat';
                %SlantColumnAmountNO2
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'SlantColumnAmountNO2')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); SlantColumnAmountNO2 = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); SlantColumnAmountNO2=double(SlantColumnAmountNO2); SlantColumnAmountNO2=SlantColumnAmountNO2';
                %ColumnAmountNO2TropStd
                datasetID = H5D.open(fileID, h5dsetname(hinfo,1,2,1,1,'ColumnAmountNO2TropStd')); dataspaceID = H5D.get_space(datasetID); H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); ColumnAmountNO2TropStd = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); ColumnAmountNO2TropStd=double(ColumnAmountNO2TropStd); ColumnAmountNO2TropStd=ColumnAmountNO2TropStd';
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
                
                % We already identified what rows to keep, so we'll reuse
                % that here 
                lon(~rows_to_keep,:)=[];                         
                lat(~rows_to_keep,:)=[];                         
                loncorn(:,~rows_to_keep,:)=[];                   
                latcorn(:,~rows_to_keep,:)=[];                   
                FoV75CornerLatitude(:,~rows_to_keep,:)=[];       
                FoV75CornerLongitude(:,~rows_to_keep,:)=[];      
                SolarAzimuthAngle(~rows_to_keep,:)=[];           
                SolarZenithAngle(~rows_to_keep,:)=[];            
                ViewingAzimuthAngle(~rows_to_keep,:)=[];         
                ViewingZenithAngle(~rows_to_keep,:)=[];          
                Time(~rows_to_keep,:)=[];                        
                AMFStrat(~rows_to_keep,:)=[];                    
                AMFTrop(~rows_to_keep,:)=[];                     
                CloudFraction(~rows_to_keep,:)=[];               
                CloudPressure(~rows_to_keep,:)=[];               
                CloudRadianceFraction(~rows_to_keep,:)=[];       
                ColumnAmountNO2(~rows_to_keep,:)=[];             
                SlantColumnAmountNO2(~rows_to_keep,:)=[];        
                TerrainHeight(~rows_to_keep,:)=[];               
                TerrainPressure(~rows_to_keep,:)=[];             
                TerrainReflectivity(~rows_to_keep,:)=[];         
                TropopausePressure(~rows_to_keep,:) = [];
                vcdQualityFlags(~rows_to_keep,:)=[];             
                XTrackQualityFlags(~rows_to_keep,:)=[];          
                RelativeAzimuthAngle(~rows_to_keep,:)=[];        
                ColumnAmountNO2Trop(~rows_to_keep,:)=[];         
                ColumnAmountNO2Strat(~rows_to_keep,:)=[];
                ColumnAmountNO2TropStd(~rows_to_keep,:) = [];
                Row(~rows_to_keep,:)=[];                         
                Pixel(~rows_to_keep,:)=[];
                Swath(~rows_to_keep,:)=[];                       
                
                if DEBUG_LEVEL > 0; disp(' Saving imported OMI fields to "Data"'); end
                %Save the imported items to the structure 'Data'.  Changed
                %on 21 May 2015 so that these matrices will retain their
                %shape
                Data(E).Latitude = lat;                                  Data(E).LatBdy = [latmin latmax];
                Data(E).Longitude = lon;                                 Data(E).LonBdy = [lonmin lonmax];
                Data(E).Loncorn = loncorn;                               Data(E).FoV75CornerLongitude = FoV75CornerLongitude;
                Data(E).Latcorn = latcorn;                               Data(E).FoV75CornerLatitude = FoV75CornerLatitude;
                Data(E).SolarAzimuthAngle = SolarAzimuthAngle;           Data(E).AMFTrop = AMFTrop;
                Data(E).SolarZenithAngle = SolarZenithAngle;             Data(E).AMFStrat = AMFStrat;
                Data(E).ViewingAzimuthAngle = ViewingAzimuthAngle;       Data(E).TerrainHeight = TerrainHeight;
                Data(E).ViewingZenithAngle = ViewingZenithAngle;         Data(E).TerrainPressure = TerrainPressure;
                Data(E).Time = Time;                                     Data(E).TerrainReflectivity = TerrainReflectivity;
                Data(E).ColumnAmountNO2 = ColumnAmountNO2;               Data(E).vcdQualityFlags = vcdQualityFlags;
                Data(E).ColumnAmountNO2Trop = ColumnAmountNO2Trop;       Data(E).SlantColumnAmountNO2 = SlantColumnAmountNO2;
                Data(E).ColumnAmountNO2TropStd = ColumnAmountNO2TropStd; Data(E).ColumnAmountNO2Strat = ColumnAmountNO2Strat;
                Data(E).CloudRadianceFraction = CloudRadianceFraction;   Data(E).CloudPressure = CloudPressure;
                Data(E).RelativeAzimuthAngle = RelativeAzimuthAngle;     Data(E).CloudFraction = CloudFraction;
                Data(E).Row = Row;                                       Data(E).XTrackQualityFlags = XTrackQualityFlags;
                Data(E).Swath = Swath;                                   Data(E).Date=date;
                Data(E).TropopausePressure = TropopausePressure;         Data(E).Pixel=Pixel;
                
                %Add MODIS cloud info to the files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if DEBUG_LEVEL > 0; fprintf('\n Adding MODIS cloud data \n'); end
                
                %Convert the OMI date to a Julian calendar day
                d2=1+datenum(str2double(year),str2double(month),str2double(day))-datenum(str2double(year),1,1);
                julian_day = sprintf('%03d',d2);
                
                
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
                if e < n % If there is at least one more swath, get its start time from the file name
                    next_omi_swath_time = str2double(sp_files(e+1).name(29:32));
                else % otherwise add 99 minutes to the start time for this swath
                    omi_hr = str2double(sp_files(e).name(29:30));
                    omi_min = str2double(sp_files(e).name(31:32));
                    next_omi_swath_time = 100*(floor((omi_min + 99)/60)+omi_hr) + mod(omi_min + 99,60);
                    if next_omi_swath_time > 2359; error(E.callError('modis_cloud','Next OMI swath time for MODIS cloud binning calculated to be in the next day. \nManual attention recommended')); end
                end
                
                % Initialize the mod_data structure (JLL 15 Jan 2015 -
                % Parallelization)
                mod_Data = struct('Longitude',[],'Latitude',[],'CloudFraction',[]);
                
                for ii=1:length(modis_files);
                    mod_filename=modis_files(ii).name;
                    % Skip any modis files that do not occur during the
                    % time period of the current swath
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
                            mod_Data.Longitude=[mod_Data.Longitude;Longitude];
                            mod_Data.Latitude=[mod_Data.Latitude;Latitude];
                            mod_Data.CloudFraction=[mod_Data.CloudFraction;CloudFraction];
                        end
                    end
                end
                
                %If there is no "mod_Data" available for this loop, fill
                %the regular Data field with -127. Otherwise, find all the
                %MODIS cloud pixels in each OMI pixel and average them
                %together.
                if isempty(mod_Data.Longitude)
                    Data(E).MODISCloud=behr_fill_val()*ones(size(Data(E).Latitude));
                else
                    Data(E).MODISCloud=nan(size(Data(E).Latitude));
                    for jj=1:numel(Data(E).Latitude);
                        x1 = Data(E).Loncorn(1,jj);   y1 = Data(E).Latcorn(1,jj);
                        x2 = Data(E).Loncorn(2,jj);   y2 = Data(E).Latcorn(2,jj);
                        x3 = Data(E).Loncorn(3,jj);   y3 = Data(E).Latcorn(3,jj);
                        x4 = Data(E).Loncorn(4,jj);   y4 = Data(E).Latcorn(4,jj);
                        
                        xall=[x1;x2;x3;x4;x1];
                        yall=[y1;y2;y3;y4;y1];
                        xx_cld = inpolygon(mod_Data.Latitude,mod_Data.Longitude,yall,xall);
                        
                        cld_vals=mod_Data.CloudFraction(xx_cld);
                        cld_vals(isnan(cld_vals))=[];
                        
                        Data(E).MODISCloud(jj)=mean(cld_vals);
                        Data(E).MODIS_Cloud_File=mod_filename;
                        
                    end
                end
                
                % "clear" cannot be used in parfor loops, so instead we'll
                % set all the arrays to empty to reduce the memory load.
                mod_Data(j).Longitude = [];
                mod_Data(j).Latitude = [];
                mod_Data(j).CloudFraction = [];
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
                        error(E.filenotfound('MCD43C3 (Albedo) file within 21 days'));
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
                    x1 = Data(E).Loncorn(1,k);   y1 = Data(E).Latcorn(1,k);
                    x2 = Data(E).Loncorn(2,k);   y2 = Data(E).Latcorn(2,k);
                    x3 = Data(E).Loncorn(3,k);   y3 = Data(E).Latcorn(3,k);
                    x4 = Data(E).Loncorn(4,k);   y4 = Data(E).Latcorn(4,k);
                    
                    
                    xall=[x1;x2;x3;x4;x1];
                    yall=[y1;y2;y3;y4;y1];
                    
                    xx_alb = inpolygon(band3_lats,band3_lons,yall,xall);
                    
                    band3_vals=band3(xx_alb);  band3_zeros=band3_vals==0;
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
                    % Since GLOBE data is on a grid where a row of
                    % latitudinal points all have the same longitude and
                    % vice versa, we can quickly reduce the number of
                    % points by comparing just one lat and lon vector to
                    % the extent of the pixel.
                    ai=find(globe_lat_matrix(:,1)>=min(yall) & globe_lat_matrix(:,1)<=max(yall));
                    bi=find(globe_lon_matrix(1,:)>=min(xall) & globe_lon_matrix(1,:)<=max(xall));
                    pressurex=terpres(ai,bi);
                    pressure_latx=globe_lat_matrix(ai,bi);
                    pressure_lonx=globe_lon_matrix(ai,bi);
                    %%%%%%%%%%%%%%%%%%%
                    
                    % inpolygon is slow compared to a simple logical test,
                    % so we only apply it to the subset of GLOBE heights
                    % immediately around our pixel.
                    xx_globe = inpolygon(pressure_latx,pressure_lonx,yall,xall);
                    
                    pres_vals=pressurex(xx_globe);  
                    GLOBETerpres(k)=1013.25 .* exp(-mean(pres_vals) / 7400 ); %Originally divided by 7640 m
                    
                    
                    if DEBUG_LEVEL > 2; telap = toc; fprintf('Time for GLOBE --> pixel %u/%u = %g sec \n',k,c,telap); end
                end
                
                Data(E).GLOBETerpres = GLOBETerpres;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            end %End the section carried out only if there are OMI pixels in the area of interest
        
        end %End the loop over all swaths in a day
        savename=[satellite,'_',retrieval,'_',year,month,day];
        saveData(fullfile(sp_mat_dir,savename), Data); % Saving must be handled as a separate function in a parfor loop because passing a variable name as a string upsets the parallelization monkey (it's not transparent).
 %       toc
 %       t=toc;
 %       if t>1200
 %           %error('Time exceeded 20 min. Stopping')
 %       end
    end %End the section checking if there are OMI files for the given time period
end %End the loop over all days
end 

function saveData(filename,Data)
    save(filename,'Data')
end
