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

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DEPENDENCIES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
if DEBUG_LEVEL > 1; fprintf('Adding folders\n'); end
%Add the 'Utils' folder and all subfolders to MATLAB's search path. Within
%the Git repository for BEHR, this is the /Utils folder.
mpath = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(mpath,'..','Utils')));


% Add the paths needed to run on the cluster
if onCluster;
    addpath(genpath('~/MATLAB/Classes'));
    addpath(genpath('~/MATLAB/Utils'));
else
    addpath(genpath(BEHR_paths('classes')));
    addpath(genpath(BEHR_paths('utils')));
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
%****************************%

% Set to 1 to overwrite existing files in the time range given,
% set to 0 to only produce missing files.
overwrite = 1;

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
    if DEBUG_LEVEL > 1; fprintf('Setting data paths from global variables\n'); end
    
    global sp_mat_dir;
    global omi_he5_dir;
    global omi_pixcor_dir
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
        error(E.callError('bad_cluster_path',sprintf(msg,nonexistant{:})));
    end
    
    
else
    %This is the directory where the final .mat file will be saved.
    sp_mat_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Production tests/v2-1D';%BEHR_paths('sp_mat_dir');
    
    %This is the directory where the OMI NASA SP (OMNO2) he5 files are
    %saved. It should include subfolders organized by year, which in turn
    %are organized by month.
    omi_he5_dir = BEHR_paths('omno2_dir');
    
    %This is the directory where the OMPIXCOR he5 files are saved. It
    %should include subfolders organized by year, which in turn are
    %organized by month.
    omi_pixcor_dir = BEHR_paths('ompixcor_dir');
    
    %This is the directory where the MODIS myd06_L2*.hdf files are saved.
    %It should include subfolders organized by year.
    modis_myd06_dir = BEHR_paths('myd06_dir');
    
    %This is the directory where the MODIS MCD43C3*.hdf files are saved. It
    %should include subfolders organized by year.
    modis_mcd43_dir = BEHR_paths('mcd43c3_dir');
    
    %This is the directory where the GLOBE data files and their headers
    %(.hdr files) are saved.
    globe_dir = BEHR_paths('globe_dir');
end



%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN BODY %%%%%
%%%%%%%%%%%%%%%%%%%%%

%Go ahead and load the terrain pressure data - only need to do this once
%Add a little buffer around the edges to make sure we have terrain data
%everywhere that we have NO2 profiles.
glonlim = [lonmin - 10, lonmax + 10];
glatlim = [latmin - 10, latmax + 10];
[terpres, globe_lon_matrix, globe_lat_matrix] = load_globe_alts(glonlim, glatlim);
terpres(isnan(terpres)) = 0;

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

datenums = datenum(date_start):datenum(date_end);

%parfor(j=1:length(datenums), n_workers)
for j=1:length(datenums)
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
    savename = sp_savename(this_dnum);
    if exist(fullfile(sp_mat_dir, savename), 'file') && ~overwrite
        if DEBUG_LEVEL > 0; fprintf('File %s exists, skipping this day\n', savename); end
        continue
    end
    
    %Prepare a data structure to receive the final data.
    % As I understand this currently (15 Jan 2015) "Data" should be a local
    % variable on each worker for parallelization, and this line will reset
    % its value within each loop. Thus we shouldn't need to worry about
    % cross-communication between workers.
    sp_variables = {'Longitude', 'Latitude', 'Time', 'ViewingZenithAngle',...
        'SolarZenithAngle', 'AmfStrat', 'AmfTrop', 'CloudFraction',...
        'CloudRadianceFraction', 'TerrainHeight', 'TerrainPressure',...
        'TerrainReflectivity', 'CloudPressure', 'ColumnAmountNO2',...
        'SlantColumnAmountNO2', 'ColumnAmountNO2Trop', 'VcdQualityFlags',...
        'XTrackQualityFlags'};
    
    pixcor_variables = {'TiledArea', 'TiledCornerLongitude', 'TiledCornerLatitude',...
        'FoV75Area', 'FoV75CornerLongitude', 'FoV75CornerLatitude'};
    
    behr_variables = {'Date', 'LatBdy', 'LonBdy', 'Row', 'MODISCloud', 'MODISCloudFiles',...
        'MODISAlbedo', 'MODIS_Albedo_File', 'GLOBETerpres'};
    
    sub_data = make_empty_struct_from_cell([sp_variables, pixcor_variables, behr_variables],0);
    Data = make_empty_struct_from_cell([sp_variables, pixcor_variables, behr_variables],0);
    %Set the file path and name, assuming that the file structure is
    %<he5_directory>/<year>/<month>/...files...  Then figure out how many
    %files there are
    short_filename = sprintf('OMI-Aura_L2-OMNO2_%04dm%02d%02d*.he5', this_year, this_month, this_day);
    file_dir = fullfile(omi_he5_dir, this_year_str, this_month_str); %Used both here to find all he5 files and in the swath for loop to identify each file.
    file_pattern=fullfile(file_dir,short_filename);
    sp_files = dir(file_pattern);
    sp_files = remove_duplicate_orbits(sp_files);
    n = length(sp_files);
    E=0;
    if isempty(sp_files);
        fprintf('No data available for %s\n', datestr(this_dnum));
        continue
    end
    
    for e=1:n %For loop over all the swaths in a given day.
        if DEBUG_LEVEL > 0
            if e==1 || mod(e,10)==0; fprintf('Swath %u of %s \n', e, datestr(this_dnum)); end
        end
        %Read in each file, saving the hierarchy as 'hinfo'
        this_sp_filename = fullfile(omi_he5_dir, this_year_str, this_month_str, sp_files(e).name);
        [omi_starttime, omi_next_starttime] = get_omi_swath_times(sp_files, e);
        
        if omi_starttime < 1500 && strcmp(region,'US')
            %If start time is < 1500 and we want to look at the US, reject
            %the file, as it is probably descending nodes only.
            if DEBUG_LEVEL > 0; fprintf(' Swath %d: Nighttime granule skipped\n',e); end
            continue
        end
        
        [this_data, pixels_in_domain] = read_omi_sp(this_sp_filename, '/HDFEOS/SWATHS/ColumnAmountNO2', sp_variables, sub_data, [lonmin, lonmax], [latmin, latmax]);
        
        if ~pixels_in_domain
            if DEBUG_LEVEL > 1; disp('No points within lat/lon boundaries'); end
            continue
        end
        
        
        % If we've gotten here, then there are pixels in the swath that lie
        % within the domain of interest. Add the OMPIXCOR data
        pixcor_name = make_pixcorn_name_from_sp(this_sp_filename, fullfile(omi_pixcor_dir, this_year_str, this_month_str));
        this_data = read_omi_sp(pixcor_name, '/HDFEOS/SWATHS/OMI Ground Pixel Corners VIS', pixcor_variables, this_data, [lonmin, lonmax], [latmin, latmax]);
        
        
        % Add MODIS cloud info to the files 
        if DEBUG_LEVEL > 0; fprintf('\n Adding MODIS cloud data \n'); end
        
        this_data = read_modis_cloud(modis_myd06_dir, this_dnum, this_data, omi_starttime, omi_next_starttime, [lonmin, lonmax], [latmin, latmax], 'DEBUG_LEVEL', DEBUG_LEVEL);
        
        
        %Add MODIS albedo info to the files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if DEBUG_LEVEL>0; fprintf('\n Adding MODIS albedo information \n'); end
        alb_dir = fullfile(modis_mcd43_dir,this_year);
        
        %Find the closest MCD file
        in=[0 1 -1 2 -2 3 -3 4 -4 5 -5 6 -6 7 -7 8 -8 9 -9 10 -10 11 -11 12 -12 13 -13 14 -14 15 -15 16 -16 17 -17 18 -18 19 -19 20 -20 21 -21];
        for ii=1:length(in);
            mcd_date = num2str(str2double(julian_day) + in(ii),'%03g');
            alb_filename = fullfile(alb_dir,['MCD43C3.A',this_year,mcd_date,'*.hdf']);
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
        
        
    end %End the loop over all swaths in a day
    saveData(fullfile(sp_mat_dir,savename), Data); % Saving must be handled as a separate function in a parfor loop because passing a variable name as a string upsets the parallelization monkey (it's not transparent).
    
    
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
    if next_swath_time > 2359; 
        E.notimplemented('Next OMI swath time calculated to be in the next day for file %s.', sp_files(current_file_index).name); 
    end
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
