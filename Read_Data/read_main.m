function read_main(varargin)
% READ_MAIN Reads in OMI, MODIS, and GLOBE data to .mat files
%
%   READ_MAIN is the first step in the BEHR workflow. It reads
%   in the satellite data from the various sources, include OMI NO2, MODIS
%   clouds, MODIS albedo, and GLOBE (a database, not a satellite) terrain
%   elevation. These data are cut down to the US domain and, in the case of
%   the MODIS and GLOBE data, averaged to the OMI pixels. The resulting
%   Data structures are saved as an OMI_SP .mat file.
%
%   This function is setup such that running it without arguments will
%   produce any new OMI_SP files required. This requires that the necessary
%   data be available either locally or via a mounted network drive. This
%   behavior can be changed with the following parameters:
%
%       'start' - set the first day to process using either a date string
%       that Matlab understands or a date number. Default is '2005-01-01'.
%
%       'end' - set the last day to process using either a date string that
%       Matlab understands or a date number. Default is today.
%
%       'sp_mat_dir' - the directory that the OMI_SP .mat files will be
%       saved to. Default is the path provided by the behr_paths class.
%
%       'omi_he5_dir' - the directory that contains the OMI NO2 HDF5 files,
%       sorted into subdirectories by year and month (i.e. this directory
%       itself should contain subdirectories 2005, 2006, etc., each of
%       which has subdirectories 01, 02, 03, etc. that contain the HDF5
%       files). Default is the path provided by the behr_paths class.
%
%       'modis_myd06_dir' - the directory that contains the MODIS MYD06
%       cloud HDF4 files, sorted by year. Default is the path provided by
%       the behr_paths class.
%
%       'modis_mcd43_dir' - the directory that contains the MODIS MCD43C1
%       BRDF parameters files, sorted into subdirectories by year. Default
%       is the path provided by the behr_paths class.
%
%       'globe_dir' - the directory that contains the GLOBE (Global Land
%       One-km Base Elevation) terrain elevation data. This will contain
%       files a10g through p10g and a10g.hdr through p10g.hdr. Default is
%       the path provided by the behr_paths class.
%
%       'region' - which region BEHR is running in. This controls both the
%       longitude and latitude limits and which orbits are skipped as
%       "nighttime" orbits. This must be a string. Default (and only option
%       at present) is 'US'.
%
%       'allow_no_myd' - boolean (default false) which allows the run to
%       process days for which no MODIS cloud fraction data is available.
%
%       'overwrite' - scalar logical which controls whether existing files
%       will be overwritten. If false, a day will be skipped if the
%       corresponding OMI_SP .mat file exists in the directory given as
%       'omi_he5_dir'. If true, no days will be skipped and the data in
%       omi_he5_dir will be overwritten.
%
%       'DEBUG_LEVEL' - verbosity. Default is 2; i.e. most progress
%       message, but no timing messages will be printed. 0 = no messages;
%       greater means more messages.

%****************************%
% CONSOLE OUTPUT LEVEL - 0 = none, 1 = minimal, 2 = all messages, 3 = times

%****************************%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION & INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

p = inputParser;
p.addParameter('start', '');
p.addParameter('end', '');
p.addParameter('sp_mat_dir', '');
p.addParameter('omi_he5_dir', '');
p.addParameter('omi_pixcor_dir', '');
p.addParameter('modis_myd06_dir', '');
p.addParameter('modis_mcd43_dir', '');
p.addParameter('globe_dir', '');
p.addParameter('region', 'US');
p.addParameter('allow_no_myd', false);
p.addParameter('overwrite', false)
p.addParameter('DEBUG_LEVEL', 2);

p.parse(varargin{:});
pout = p.Results;

date_start = pout.start;
date_end = pout.end;
sp_mat_dir = pout.sp_mat_dir;
omi_he5_dir = pout.omi_he5_dir;
omi_pixcor_dir = pout.omi_pixcor_dir;
modis_myd06_dir = pout.modis_myd06_dir;
modis_mcd43_dir = pout.modis_mcd43_dir;
globe_dir = pout.globe_dir;
allow_no_myd = pout.allow_no_myd;
region = pout.region; 
overwrite = pout.overwrite;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

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

% Cleanup object will safely exit if there's a problem
if onCluster
    cleanupobj = onCleanup(@() mycleanup());
end

% Defined the number of threads to run, this will be used to open a
% parallel pool. numThreads should be set in the calling run script,
% otherwise it will default to 1.
global numThreads;
if isempty(numThreads)
    numThreads = 1;
end

%%%%%%%%%%%%%%%%%%%%%%
%%%%% VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%

date_start = validate_date(date_start);
date_end = validate_date(date_end);

if ~ischar(sp_mat_dir)
    E.badinput('Paramter "sp_mat_dir" must be a string')
elseif ~ischar(omi_he5_dir)
    E.badinput('Paramter "omi_he5_dir" must be a string')
elseif ~ischar(omi_pixcor_dir)
    E.badinput('Paramter "omi_pixcor_dir" must be a string')
elseif ~ischar(modis_myd06_dir)
    E.badinput('Paramter "modis_myd06_dir" must be a string')
elseif ~ischar(modis_mcd43_dir)
    E.badinput('Paramter "modis_mcd43_dir" must be a string')
elseif ~ischar(globe_dir)
    E.badinput('Paramter "globe_dir" must be a string')
elseif ~ischar(region)
    E.badinput('Paramter "region" must be a string')
elseif (~islogical(overwrite) && ~isnumeric(overwrite)) || ~isscalar(overwrite)
    E.badinput('Parameter "overwrite" must be a scalar logical or number')
end


% Specify the longitude and latitude ranges of interest for this retrieval.
% Additionally, set the earliest and latest start time (in UTC) for the
% swaths that will be allowed. This will help 
switch lower(region)
    case 'us'
        lonmin = -125;    
        lonmax = -65;
        latmin = 25;    
        latmax = 50;
        earliest_omi_starttime = 1500;
        latest_omi_starttime = Inf;
    case 'hk'
        lonmin = 108;
        lonmax = 118;
        latmin = 19;
        latmax = 26;
        earliest_omi_starttime = -Inf;
        latest_omi_starttime = 1300;
    otherwise 
        E.badinput('Region "%s" not recognized', region)
end

if lonmin > lonmax %Just in case I enter something backwards...
    E.callError('bounds', 'Lonmin is greater than lonmax')
elseif latmin > latmax
    E.callError('bounds', 'Latmin is greater than latmax')
end

%Process all files between these dates, in yyyy/mm/dd format unless
%overriding dates are passed into the function.
%****************************%
if isempty(date_start) || isempty(date_end)
    date_start='2013/08/01';
    date_end='2013/08/06';
end
%****************************%

% This helps preallocate the Data structure, which is generally more
% efficient than expanding each time it's needed. Ideally this should be
% the maximum number of swaths expected for a given domain; extra swaths
% will be removed at the end.
estimated_num_swaths = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DATA DIRECTORIES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These are the directories to save or load data from. By default, they are
% taken from the behr_paths class, which can be created by running
% BEHR_initial_setup in the root folder of this repo. Alternately, paths
% can be specified as parameter inputs to this function (useful on the
% cluster where this should be called from a run script or in unit tests
% where the default save directories need to be overwritten).

%This is the directory where the final .mat file will be saved.
if isempty(sp_mat_dir)
    sp_mat_dir = behr_paths.SPMatSubdir(region);
end

%This is the directory where the OMI NASA SP (OMNO2) he5 files are
%saved. It should include subfolders organized by year, which in turn
%are organized by month.
if isempty(omi_he5_dir)
    omi_he5_dir = behr_paths.omno2_dir;
end

%This is the directory where the OMPIXCOR he5 files are saved. It
%should include subfolders organized by year, which in turn are
%organized by month.
if isempty(omi_pixcor_dir)
    omi_pixcor_dir = behr_paths.ompixcor_dir;
end

%This is the directory where the MODIS myd06_L2*.hdf files are saved.
%It should include subfolders organized by year.
if isempty(modis_myd06_dir)
    modis_myd06_dir = behr_paths.myd06_dir;
end

%This is the directory where the MODIS MCD43C3*.hdf files are saved. It
%should include subfolders organized by year.
if isempty(modis_mcd43_dir)
    modis_mcd43_dir = behr_paths.mcd43d_dir;
end

%This is the directory where the GLOBE data files and their headers
%(.hdr files) are saved.
if isempty(globe_dir)
    globe_dir = behr_paths.globe_dir;
end

% Verify the paths integrity.
nonexistant = {};

if ~exist(sp_mat_dir,'dir')
    nonexistant{end+1} = 'sp_mat_dir';
end
if ~exist(omi_he5_dir,'dir')
    nonexistant{end+1} = 'he5_dir';
end
if ~exist(omi_pixcor_dir, 'dir')
    nonexistant{end+1} = 'ompixcor_dir';
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
    E.callError('bad_paths',sprintf(msg,nonexistant{:}));
end




%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN BODY %%%%%
%%%%%%%%%%%%%%%%%%%%%

%Add a little buffer around the edges to make sure we have ancillary data
%everywhere that we have NO2 profiles.
ancillary_lonlim = [lonmin - 10, lonmax + 10];
ancillary_latlim = [latmin - 10, latmax + 10];

%Load the land classification map. We'll use this to decide when to use the
%ocean surface reflectance parameterization. 
if DEBUG_LEVEL > 1; fprintf('Loading land/ocean classification map\n'); end
if DEBUG_LEVEL > 2; t_load_land_ocean = tic; end

[ocean_mask.mask, ocean_mask.lon, ocean_mask.lat] = get_modis_ocean_mask(ancillary_lonlim, ancillary_latlim);

if DEBUG_LEVEL > 2; fprintf('    Time to load land/ocean classification map: %f\n', toc(t_load_land_ocean)); end

%Go ahead and load the terrain pressure data - only need to do this once
if DEBUG_LEVEL > 1; fprintf('Loading globe elevations\n'); end
if DEBUG_LEVEL > 2; t_load_globe = tic; end

[globe_elevations, globe_lon_matrix, globe_lat_matrix] = load_globe_alts(ancillary_lonlim, ancillary_latlim);
globe_elevations(isnan(globe_elevations)) = 0;
if DEBUG_LEVEL > 2; fprintf('    Time to load GLOBE elevations: %f\n', toc(t_load_globe)); end

if DEBUG_LEVEL > 1; fprintf('Loading COART sea reflectances\n'); end
if DEBUG_LEVEL > 2; t_load_coart = tic; end
[~, coart_lut] = coart_sea_reflectance(0);
if DEBUG_LEVEL > 2; fprintf('    Time to load COART look up table: %f\n', toc(t_load_coart)); end

%For loop over all days from the starting or last finished date to the end
%date. We will give the absolute paths to files rather than changing the
%active directory, as MATLAB seems to run slightly slower if the current
%working directory is on the server.
% total_days=datenum(date_end)-datenum(last_date)+1;
% for j=1:total_days;

if onCluster
    n_workers=numThreads;
else
    n_workers=0;
end

if onCluster && isempty(gcp('nocreate'))
    parpool(numThreads);
end

onCluster_local = onCluster;

% Setup some values that either need to be computed to determine the loop
% indices or which are better calculated outside the loop.

datenums = datenum(date_start):datenum(date_end);
core_githead = git_head_hash(behr_paths.behr_core);
behrutils_githead = git_head_hash(behr_paths.behr_utils);
genutils_githead = git_head_hash(behr_paths.utils);
behr_grid = GlobeGrid(0.05, 'domain', [lonmin, lonmax, latmin, latmax]);


if DEBUG_LEVEL > 1; fprintf('Staring main loop\n'); end

t_comm = tic;
parfor(j=1:length(datenums), n_workers)
%for j=1:length(datenums)
    this_task = getCurrentTask();
    if isempty(this_task)
        this_task.ID = -1;
    end

    if DEBUG_LEVEL > 2; fprintf('Worker %d: Time to enter parfor loop: %f s\n', this_task.ID, toc(t_comm)); end
    if DEBUG_LEVEL > 2; t_day = tic; end

    %Read the desired year, month, and day
    this_dnum = datenums(j);
    this_year = year(this_dnum);
    this_year_str = sprintf('%04d', this_year);
    this_month=month(this_dnum);
    this_month_str = sprintf('%02d', this_month);
    this_day=day(this_dnum);
    this_day_str = sprintf('%02d', this_day);
    
    % Check if the file already exists. If it does, and if we're set
    % to not overwrite, we don't need to process this day.
    savename = sp_savename(this_dnum, region);
    if exist(fullfile(sp_mat_dir, savename), 'file') && ~overwrite
        if DEBUG_LEVEL > 0; fprintf('File %s exists, skipping this day\n', savename); end
        continue
    end
    
    % List variables that should be read directly from the OMI OMNO2 files.
    % If you want an additional variable, adding it here should be
    % sufficient to make that happen as long as it is under the
    % /HDFEOS/SWATHS/ColumnAmountNO2 group and is spelled exactly how the
    % dataset is named.
    sp_variables = {'Longitude', 'Latitude', 'SpacecraftAltitude', 'SpacecraftLatitude',...
        'SpacecraftLongitude', 'Time', 'ViewingZenithAngle',...
        'SolarZenithAngle', 'ViewingAzimuthAngle', 'SolarAzimuthAngle',...
        'AmfStrat', 'AmfTrop', 'CloudFraction', 'CloudRadianceFraction',...
        'TerrainHeight', 'TerrainPressure', 'TerrainReflectivity',...
        'CloudPressure', 'ColumnAmountNO2', 'SlantColumnAmountNO2',...
        'ColumnAmountNO2Trop', 'ColumnAmountNO2TropStd', 'ColumnAmountNO2Strat',...
        'TropopausePressure', 'VcdQualityFlags', 'XTrackQualityFlags'};
    
    % Variables from the OMPIXCOR files. As with the SP variables, adding a
    % variable here that is under the '/HDFEOS/SWATHS/OMI Ground Pixel
    % Corners VIS' group should be sufficient to cause it to be read in.
    pixcor_variables = {'TiledArea', 'TiledCornerLongitude', 'TiledCornerLatitude',...
        'FoV75Area', 'FoV75CornerLongitude', 'FoV75CornerLatitude'};
    
    % Variables that will be added by BEHR. These will need manual
    % intervention if you choose to add more variables since they're not
    % being copied directly from existing files.
    behr_variables = {'Date', 'Grid', 'LatBdy', 'LonBdy', 'Row', 'Swath', 'RelativeAzimuthAngle',...
        'MODISCloud',  'MODISAlbedo', 'MODISAlbedoQuality','MODISAlbedoFillFlag', 'GLOBETerpres',...
        'IsZoomModeSwath', 'AlbedoOceanFlag','OMPIXCORFile', 'MODISCloudFiles', 'MODISAlbedoFile',...
        'GitHead_Core_Read', 'GitHead_BEHRUtils_Read', 'GitHead_GenUtils_Read', 'OMNO2File', 'BEHRRegion'};
    
    sub_data = make_empty_struct_from_cell([sp_variables, pixcor_variables, behr_variables],0);
    Data = repmat(make_empty_struct_from_cell([sp_variables, pixcor_variables, behr_variables],0), 1, estimated_num_swaths);
    
    %Set the file path and name, assuming that the file structure is
    %<he5_directory>/<year>/<month>/...files...  Then figure out how many
    %files there are
    short_filename = sprintf('OMI-Aura_L2-OMNO2_%04dm%02d%02d*.he5', this_year, this_month, this_day);
    file_dir = fullfile(omi_he5_dir, this_year_str, this_month_str); %Used both here to find all he5 files and in the swath for loop to identify each file.
    file_pattern=fullfile(file_dir,short_filename);
    sp_files = dir(file_pattern);
    sp_files = remove_duplicate_orbits(sp_files);
    n = length(sp_files);
    
    
    if isempty(sp_files);
        fprintf('No data available for %s\n', datestr(this_dnum));
        continue
    end
    
    % Read in the MODIS albedo data for this day. We do it outside the loop
    % over orbits to limit the number of reads of the (fairly large) MCD43D
    % files.
    if DEBUG_LEVEL > 2; t_alb_read = tic; end
    modis_brdf_data = read_modis_albedo(modis_mcd43_dir, this_dnum, ancillary_lonlim, ancillary_latlim);
    if DEBUG_LEVEL > 2; fprintf('Worker %d: Time to read MODIS BRDF = %f\n', this_task.ID, toc(t_alb_read)); end
    
    data_ind = 0;
    for a=1:n %For loop over all the swaths in a given day.
        if DEBUG_LEVEL > 2; t_orbit = tic; end

        if DEBUG_LEVEL > 0
            if a==1 || mod(a,10)==0; fprintf('Swath %u of %s \n', a, datestr(this_dnum)); end
        end
        %Read in each file, saving the hierarchy as 'hinfo'
        this_sp_filename = fullfile(omi_he5_dir, this_year_str, this_month_str, sp_files(a).name);
        [omi_starttime, omi_next_starttime] = get_omi_swath_times(sp_files, a);
        
        if omi_starttime < earliest_omi_starttime || omi_starttime > latest_omi_starttime
            %If start time is < 1500 and we want to look at the US, reject
            %the file, as it is probably descending nodes only.
            if DEBUG_LEVEL > 0; fprintf(' Swath %d: Nighttime granule skipped\n',a); end
            continue
        end
       
        if DEBUG_LEVEL > 2; t_sp = tic; end 
        [this_data, pixels_in_domain] = read_omi_sp(this_sp_filename, '/HDFEOS/SWATHS/ColumnAmountNO2', sp_variables, sub_data, [lonmin, lonmax], [latmin, latmax]);
        if DEBUG_LEVEL > 2; fprintf('      Time to read SP data on worker %d: %f\n', this_task.ID, toc(t_sp)); end
        
        if ~pixels_in_domain
            if DEBUG_LEVEL > 1; disp('No points within lat/lon boundaries'); end
            continue
        end
        
        this_data.OMNO2File = this_sp_filename;
        
        % If we've gotten here, then there are pixels in the swath that lie
        % within the domain of interest. Add the OMPIXCOR data
        if DEBUG_LEVEL > 2; t_pixcor = tic; end
        pixcor_name = make_pixcorn_name_from_sp(this_sp_filename, fullfile(omi_pixcor_dir, this_year_str, this_month_str));
        this_data = read_omi_sp(pixcor_name, '/HDFEOS/SWATHS/OMI Ground Pixel Corners VIS', pixcor_variables, this_data, [lonmin, lonmax], [latmin, latmax], 'match_data', true);
        if DEBUG_LEVEL > 2; fprintf('      Time to read OMPIXCOR data on worker %d: %f\n', this_task.ID, toc(t_pixcor)); end
        
        this_data.OMPIXCORFile = pixcor_name;
        
        if DEBUG_LEVEL > 2; t_pixclean = tic; end
        this_data = handle_corner_zeros(this_data, DEBUG_LEVEL);
        
        % The OMNO2 and OMPIXCOR products place the zoom-mode pixels
        % differently in the swath - fix that here.
        this_data = align_zoom_mode_pixels(this_data);
        if DEBUG_LEVEL > 2; fprintf('      Time to clean up pixel corners on worker %d: %f\n', this_task.ID, toc(t_pixclean)); end        


        % Add a few pieces of additional information
        % Swath is given in the file name as a five digits number following
        % the "-o" in the name (swath == orbit number)
        swath = str2double(regexp(this_sp_filename, '(?<=-o)\d\d\d\d\d', 'match', 'once'));
        this_data.Swath = swath;
        
        raa_tmp=abs(this_data.SolarAzimuthAngle + 180 - this_data.ViewingAzimuthAngle); % the extra factor of 180 corrects for the definition of RAA in the scattering weight lookup table
        raa_tmp(raa_tmp > 180)=360-raa_tmp(raa_tmp > 180);
        this_data.RelativeAzimuthAngle = raa_tmp;
        
        % Add our calculated lat and lon corners. This should only be
        % temporary until the unit testing to verify that the version 3
        % read function produces identical files to the version 2 read
        % function is complete, since the MODIS and GLOBE variables rely on
        % the lat/lon corners to be averaged to the pixel.
        %this_data = add_behr_corners(this_data, this_sp_filename);
        
        % Add MODIS cloud info to the files 
        if DEBUG_LEVEL > 0; fprintf('\n Adding MODIS cloud data \n'); end
        
        if DEBUG_LEVEL > 2; t_modis_cld = tic; end
        this_data = read_modis_cloud(modis_myd06_dir, this_dnum, this_data, omi_starttime, omi_next_starttime, [lonmin, lonmax], [latmin, latmax],...
            'AllowNoFile', allow_no_myd, 'DEBUG_LEVEL', DEBUG_LEVEL);
        if DEBUG_LEVEL > 2; fprintf('      Time to average MODIS clouds on worker %d: %f\n', this_task.ID, toc(t_modis_cld)); end
        
        
        % Add MODIS albedo info to the files
        if DEBUG_LEVEL > 0; fprintf('\n Adding MODIS albedo information \n'); end
        if DEBUG_LEVEL > 2; t_modis_alb = tic; end
        this_data = avg_modis_alb_to_pixels(modis_brdf_data, coart_lut, ocean_mask, this_data, 'QualityLimit', 3, 'DEBUG_LEVEL', DEBUG_LEVEL);

        if DEBUG_LEVEL > 2; fprintf('      Time to average MODIS albedo on worker %d: %f\n', this_task.ID, toc(t_modis_alb)); end
        
        % Add GLOBE terrain pressure to the files
        if DEBUG_LEVEL > 0; fprintf('\n Adding GLOBE terrain data \n'); end
        if DEBUG_LEVEL > 2; t_globe = tic; end
        this_data = avg_globe_data_to_pixels(this_data, globe_elevations, globe_lon_matrix, globe_lat_matrix,...
            'DEBUG_LEVEL', DEBUG_LEVEL);
        if DEBUG_LEVEL > 2; fprintf('      Time to average GLOBE data on worker %d: %f\n', this_task.ID, toc(t_globe)); end
        
        % Add the few attribute-like variables
        this_data.Date = datestr(this_dnum, 'yyyy/mm/dd');
        this_data.LonBdy = [lonmin, lonmax];
        this_data.LatBdy = [latmin, latmax];
        this_data.GitHead_Core_Read = core_githead;
        this_data.GitHead_BEHRUtils_Read = behrutils_githead;
        this_data.GitHead_GenUtils_Read = genutils_githead;
        this_data.Grid = behr_grid;
        this_data.BEHRRegion = lower(region);
        
        data_ind = data_ind + 1;
        Data(data_ind) = this_data;

        if DEBUG_LEVEL > 2; fprintf('    Time for one orbit on worker %d: %f\n', this_task.ID, toc(t_orbit)); end
        
    end %End the loop over all swaths in a day
    
    % Remove preallocated but unused swaths
    Data(data_ind+1:end) = [];
    if DEBUG_LEVEL > 2; t_save = tic; end
    saveData(fullfile(sp_mat_dir,savename), Data); % Saving must be handled as a separate function in a parfor loop because passing a variable name as a string upsets the parallelization monkey (it's not transparent).
    if DEBUG_LEVEL > 2; fprintf('       Time to save on worker %d: %f\n', this_task.ID, toc(t_save)); end
    
    if DEBUG_LEVEL > 2; fprintf('    Time for one day on worker %d: %f\n', this_task.ID, toc(t_day)); end
end %End the loop over all days
end

function saveData(filename,Data)
save(filename,'Data')
end

function mycleanup()
err=lasterror;
if ~isempty(err.message)
    fprintf('MATLAB exiting due to problem: %s\n', err.message);
    if ~isempty(gcp('nocreate'))
        delete(gcp)
    end
    
    exit(1)
end
end

function [this_swath_time, next_swath_time] = get_omi_swath_times(sp_files, current_file_index)
% Calculate the start time for the next OMI swath.
% Usually there is at least one swath starting after
% the one overflying the US west coast for that day,
% but occasionally that swath is not present in the
% OMNO2 data.  In those cases, we need to calculate the
% time it should have started, knowing that the Aura
% orbit period is ~99 min.  If for some reason the
% calculated start time should end up being in the next
% day, error out so the user is aware of that.

E = JLLErrors;
% Find the orbit time in the file name: look for four numbers (hhmm)
% preceeded by "t" and succeeded by "-o". The file names have the format,
% e.g. OMI-Aura_L2-OMNO2_2013m0808t1715-o48223_v003-2013m0809t125937.he5,
% where in this case 2013m0808t1715 gives the date and UTC time at the
% beginning of the swath as Aug 8th, 2013, 17:15 UTC. Checking that it's
% followed by -o is necessary to distinguish from the processing time of
% 2013m0809t125937,
time_regex = '(?<=t)\d\d\d\d(?=-o)'; 
this_swath_time = regexp(sp_files(current_file_index).name, time_regex, 'match', 'once');
this_swath_time = str2double(this_swath_time);

if current_file_index < numel(sp_files) % If there is at least one more swath, get its start time from the file name
    next_swath_time = regexp(sp_files(current_file_index+1).name, time_regex, 'match', 'once');
    next_swath_time = str2double(next_swath_time);
else % otherwise add 99 minutes to the start time for this swath
    omi_hr = floor(this_swath_time/100);
    omi_min = mod(this_swath_time, 100);
    next_swath_time = 100*(floor((omi_min + 99)/60)+omi_hr) + mod(omi_min + 99,60);
end


end

function sp_files = remove_duplicate_orbits(sp_files)
% From Oct 23, 2013 to Nov 4, 2013, several orbits have two files.  Since I
% could find nothing about this, I'm assuming that the one with the later
% processing date is the best one to use.  This subfunction will take a
% structure of OMNO2 files returned from dir() and check for duplicate
% orbits. If any are found, a warning is issued and only the one with the
% most recent processing date is kept.

orbits = zeros(size(sp_files));
for a=1:numel(sp_files)
    % The orbit is identified in the file name as "-o#####"
    orbits(a) = str2double(sp_files(a).name(35:39));
end

uorbits = unique(orbits);
if numel(uorbits) == numel(orbits)
    % All orbits only present once, return now.
    return
end

% Otherwise we need to identify the doubled orbits
rm_bool = false(size(sp_files)); % will be set to true for the orbits to be removed
for a=1:numel(uorbits)
    % For each orbit that has >1 file, find the one with the most recent
    % processing time; the rest for that orbit will be removed.
    xx = find(orbits == uorbits(a));
    if numel(xx) > 1
        names = {sp_files(xx).name};
        latest_proc_ind = 0;
        latest_proc_datenum = 0;
        for b = 1:numel(names)
            proc_time = names{b}([46:49,51:54]);
            proc_datenum = datenum(proc_time,'yyyymmdd');
            if proc_datenum > latest_proc_datenum
                latest_proc_ind = b;
                latest_proc_datenum = proc_datenum;
            end
        end
        yy = true(size(xx));
        yy(latest_proc_ind) = false;
        rm_bool(xx(yy)) = true;
    end
end

% Give a warning for the log
if any(rm_bool)
    rm_files = {sp_files(rm_bool).name};
    fspec = repmat('\t%s\n',1,sum(rm_bool));
    wmsg = sprintf('Duplicate orbits detected, the following files will not be used:\n%s',fspec);
    warning(wmsg, rm_files{:});
end

sp_files(rm_bool) = [];

end

function pixcor_name = make_pixcorn_name_from_sp(sp_name, pixcor_dir_with_year_month)
% OMPIXCOR file names are the same as OMNO2 names, except OMNO2 is replaced
% by OMPIXCOR and the processing time is likely different
E = JLLErrors;
[~,sp_name] = fileparts(sp_name);
stem_regex = 'OMI-Aura_L2-OMNO2_\d\d\d\dm\d\d\d\dt\d\d\d\d-o\d\d\d\d\d';
pixcor_pattern = regexp(sp_name, stem_regex, 'match', 'once');
pixcor_pattern = strrep(pixcor_pattern, 'OMNO2', 'OMPIXCOR');
pixcor_pattern = strcat(pixcor_pattern, '*.he5');
F = dir(fullfile(pixcor_dir_with_year_month, pixcor_pattern));
if numel(F) == 1
    pixcor_name = fullfile(pixcor_dir_with_year_month, F(1).name);
    return
elseif numel(F) < 1
    E.filenotfound('Could not find OMPIXCOR file matching %s in %s', pixcor_pattern, pixcor_dir_with_year_month);
elseif numel(F) > 1
    E.toomanyfiles('Multiple files matched pattern %s in %s', pixcor_pattern, pixcor_dir_with_year_month);
end
end

function data = align_zoom_mode_pixels(data)
% As of 21 Apr 2017, during zoom mode operation, OMNO2 places the 30
% available pixels in rows 0-29. OMPIXCOR places them in rows 15-44. Yay
% consistency. This subfunction fixes that so that the corner coordinates
% lie in rows 0-29. 
%
% I will assume that any day with entire rows of NaNs in Latitude/Longitude
% is a zoom mode day.

E = JLLErrors;

corner_fields = {'FoV75CornerLongitude', 'FoV75CornerLatitude';...
                 'TiledCornerLongitude', 'TiledCornerLatitude'};
area_fields = {'FoV75Area'; 'TiledArea'};

if size(corner_fields, 2) ~= 2
    E.callError('bad_field_def', 'corner_fields must be n-by-2, with the Longitude field in the field column, latitude field in the second');
elseif size(corner_fields, 1) ~= size(area_fields,1)
    E.notimplemented('corner fields without corresponding area fields or vice versa')
end

data.IsZoomModeSwath = false;
pix_nans = all(isnan(data.Longitude),1) & all(isnan(data.Latitude),1);
pix_nans_for_zoom = 31:60;
corn_nans_for_zoom = [1:15, 46:60];
if ~any(pix_nans)
    return
end

for a=1:size(corner_fields,1)
    loncorn = data.(corner_fields{a,1});
    latcorn = data.(corner_fields{a,2});
    corn_nans = squeeze(all(all(isnan(loncorn),1),2) & all(all(isnan(latcorn),1),2));
    pixarea = data.(area_fields{a});
    area_nans = isnan(pixarea);
    
    if ~isequal(size(loncorn), [4, size(data.Longitude)]) || ~isequal(size(latcorn), [4, size(data.Latitude)])
        E.callError('zoom_corners', 'The %s and %s fields are not 4 in the first dimension, and the same size as Longitude/Latitude in the second and third', corner_fields{a,1}, corner_fields{a,2});
    end
    
    if ~isvector(pixarea)
        E.callError('zoom_area', 'The %s field is expected to be a vector; it is not', area_fields{a});
    elseif numel(pixarea) ~= size(data.Longitude,2)
        E.callError('zoom_area', 'The %s field is not the same length as Longitude in the second dimension', area_fields{a})
    end
    
    if isvecequal(find(corn_nans), corn_nans_for_zoom) && isvecequal(find(area_nans), corn_nans_for_zoom) && isvecequal(find(pix_nans), pix_nans_for_zoom)
        % If all the NaNs are where we expect them to be for zoom mode,
        % then assum this is a zoom mode swath and we need to move the
        % corner coordinates to line up with the regular pixels.
        data.IsZoomModeSwath = true;
        
        loncorn(:,:,~pix_nans) = loncorn(:,:,~corn_nans);
        loncorn(:,:,pix_nans) = nan;
        latcorn(:,:,~pix_nans) = latcorn(:,:,~corn_nans);
        latcorn(:,:,pix_nans) = nan;
        
        pixarea(~pix_nans) = pixarea(~area_nans);
        pixarea(pix_nans) = nan;
    else
        % Otherwise, we may have pixels with valid center coordinates but
        % not corners, or pixels with valid corners but not center
        % coordinates. I decided to NaN the corners of pixels with invalid
        % center coordinates but not the other way around b/c a failure of
        % the corner algorithm does not mean the NO2 retrieval is bad
        % (although it does preclude BEHR retrieval b/c we won't be able to
        % average MODIS albedo and GLOBE terrain pressure to the pixel) but
        % if the center coordinate is NaNed for some reason, I don't want
        % to do anything to suggest that that pixel is good.
        loncorn(:,:,pix_nans) = nan;
        latcorn(:,:,pix_nans) = nan;
        pixarea(pix_nans) = nan;
    end
    
    
    
    data.(corner_fields{a,1}) = loncorn;
    data.(corner_fields{a,2}) = latcorn;
    data.(area_fields{a}) = pixarea;
end

end


function data = handle_corner_zeros(data, DEBUG_LEVEL)
fns = fieldnames(data);
ff = ~iscellcontents(regexpi(fns, 'corner', 'once'), 'isempty');
fns = fns(ff);
for a=1:numel(fns)
    xx = all(data.(fns{a}) == 0, 1);
    if any(xx(:)) && DEBUG_LEVEL > 0
        fprintf('    Pixels with all corners == 0 found for field %s, setting corners to NaN\n', fns{a});
    end
    data.(fns{a})(:,xx) = nan;
end
end


function data = add_behr_corners(data, sp_filename)
spacecraft_vars = {'Longitude', 'Latitude', 'SpacecraftAltitude', 'SpacecraftLatitude', 'SpacecraftLongitude'};
spacecraft = make_empty_struct_from_cell(spacecraft_vars);
spacecraft = read_omi_sp(sp_filename, '/HDFEOS/SWATHS/ColumnAmountNO2', spacecraft_vars, spacecraft, [-180 180], [-90 90]);

% If dealing with a zoom mode day, then rows 30-59 will be NaNs in the
% pixel coordinates. NaNs mess up find_submatrix2, so we need to identify
% the non-NaN piece and work with that, then put the corners in the right
% place.
datanans = all(isnan(data.Longitude),1) & all(isnan(data.Latitude),1);
scnans = all(isnan(spacecraft.Longitude),1) & all(isnan(spacecraft.Latitude),1);
subinds = find_submatrix2(data.Longitude(:,~datanans), data.Latitude(:,~datanans), spacecraft.Longitude(:,~scnans), spacecraft.Latitude(:,~scnans));

corners = fxn_corner_coordinates(spacecraft.Latitude(:,~scnans), spacecraft.Longitude(:,~scnans), spacecraft.SpacecraftLatitude, spacecraft.SpacecraftLongitude, spacecraft.SpacecraftAltitude);

xx = subinds(1,1):subinds(1,2);
yy = subinds(2,1):subinds(2,2);
loncorn = squeeze(corners(xx,yy,1,1:4));
loncorn = permute(loncorn, [3 1 2]);
latcorn = squeeze(corners(xx,yy,2,1:4));
latcorn = permute(latcorn, [3 1 2]);

data.Loncorn = nan([4, size(data.Longitude)]);
data.Latcorn = nan([4, size(data.Latitude)]);
data.Loncorn(:,:,~datanans) = loncorn;
data.Latcorn(:,:,~datanans) = latcorn;

end
