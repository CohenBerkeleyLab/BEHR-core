function read_tempo_simulator(date_start, date_end)
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
DEBUG_LEVEL = 2;
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
    date_start='2013/06/10';
    date_end='2013/06/10';
end
%****************************%

%These will be included in the file name
%****************************%
satellite='TEMPO';
retrieval='SIM_POLYALB_US';
%****************************%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DATA DIRECTORIES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The location of the directories to read or save data to.  If onCluster is
% true, these will need to be set in the runscript - I figured this would
% be easier than setting them as shell environmental variables and using
% getenv - JLL 15 Jan 2015

if onCluster
    global sp_mat_dir;
    global tempo_base_file;
    global modis_mcd43_dir;
    global globe_dir;
    
    % Verify the paths integrity.
    nonexistant = {};
    
    if ~exist(sp_mat_dir,'dir')
        nonexistant{end+1} = 'sp_mat_dir';
    end
    if ~exist(tempo_base_file,'file')
        nonexistant{end+1} = 'tempo_base_file';
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
    sp_mat_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed';
    
    %This is the file that contains the basic Data structure with the TEMPO
    %lon/lat coordinates
    tempo_base_file = '/Users/Josh/Documents/MATLAB/BEHR/Utils/TEMPO/tempo_pix_geocoords.mat';
    
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

% What hours (in UTC) that we'll assume TEMPO will retrieve
utc_hrs = 13:22;
% The position of the TEMPO satellite (assumed for now)
tempo_lon = -100;
tempo_lat = 0;
tempo_alt = 35786; % in km

% Type of albedo to use. Options are 'black-sky' which will use the normal
% MCD43C3 albedo (fixed SZA at noon, integrated over all viewing angles).
% 'poly' will adjust for the solar zenith angle using a polynomial and the
% kernels given in the MCD43C1 product. 'brdf' will do a full kernel implementation
% that accounts for all of the viewing geometry.
alb_type = 'poly';
if ~ismember(alb_type, {'black-sky','poly','brdf'})
    E.callError('bad_setting','The albedo type is not recognized')
end

%File names will be prefixed with "<satellite>_<retrieval>_", e.g. for OMI
%satellite SP retrieval, the prefix will be "OMI_SP_" and then the date in
%year, month, date order.  This section checks to see if the last file in
%the mat directory has the expected prefix.  If so, that date is taken as
%the last date completed, otherwise it is assumed that the retrieval will
%need to start from the specified start date. This allows he5 reading to be
%stopped and restarted with minimal intervention.
file_prefix = [satellite,'_',retrieval,'_']; l = length(file_prefix);


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

datenums = datenum(date_start):datenum(date_end);

local_tempo_file = tempo_base_file;
local_modis_mcd34_dir = modis_mcd43_dir;
local_sp_dir = sp_mat_dir;

for j=1:length(datenums)
    %Read the desired year, month, and day
    R=datenums(j);
    this_date=datestr(R,26);
    year=this_date(1:4);
    month=this_date(6:7);
    day=this_date(9:10);
    
    D = load(local_tempo_file);
    Data = D.Data;
    
    % Cut down to just the area we want to deal with b/c seriously this
    % takes way too long otherwise.
    xx = any(Data.Longitude > lonmin, 1) & any(Data.Longitude < lonmax, 1);
    yy = any(Data.Latitude > latmin, 2) & any(Data.Latitude < latmax, 2);
    % Sometimes there only a few points in a column or row that aren't nans
    % that get removed by one of the criteria
    xx2 = any(~isnan(Data.Longitude(yy,:)), 1);
    yy2 = any(~isnan(Data.Latitude(:,xx)), 2);
    
    Data.Longitude = Data.Longitude(yy&yy2, xx&xx2);
    Data.Latitude = Data.Latitude(yy&yy2, xx&xx2);
    Data.Loncorn = Data.Loncorn(:, yy&yy2, xx&xx2);
    Data.Latcorn = Data.Latcorn(:, yy&yy2, xx&xx2);
    
    sz = size(Data.Longitude);
    Data.SolarZenithAngle = nan(sz);
    Data.ViewingZenithAngle = nan(sz);
    Data.RelativeAzimuthAngle = nan(sz);
    Data.CloudFraction = zeros(sz);
    Data.CloudPressure = 1013*ones(sz);
    Data.CloudRadianceFraction = zeros(sz);
    Data.MODISAlbedo = nan(sz);
    Data.MODIS_Albedo_File = '';
    Data.GLOBETerpres = nan(sz);
    
    
    n = numel(utc_hrs);
    Data = repmat(Data,n,1);
    
    for e=1:n %For loop over all the scans in a given day.
        if DEBUG_LEVEL > 0; fprintf('Adding data for %d UTC\n',utc_hrs(e)); end
        % Need surface pressure in order to get the zenith angles, which we
        % need in turn to get the MODIS albedo (at least over the ocean)
        if e == 1
            % Terrain height better not change hour to hour!
            Data(e) = addGLOBETerpress(Data(e), globe_lon_matrix, globe_lat_matrix, terpres, DEBUG_LEVEL);
        else
            Data(e).GLOBETerpres = Data(1).GLOBETerpres;
        end
        
        time = sprintf('%s-%s-%s %02d:00:00',year,month,day,utc_hrs(e));
        [Data(e).SolarZenithAngle, Data(e).ViewingZenithAngle, Data(e).RelativeAzimuthAngle] = sat_angles(Data(e).Longitude, Data(e).Latitude, Data(e).GLOBETerpres, tempo_lon, tempo_lat, tempo_alt, time);
        
        if e == 1 || ~strcmpi(alb_type,'black-sky')
            Data(e) = addMODISAlbedo(Data(e), datenums(j), local_modis_mcd34_dir, alb_type, DEBUG_LEVEL);
        else
            % Similarly, if we don't use a BRDF product, albedo
            % shouldn't change with SZA. (That may be something to think
            % about too though.)
            Data(e).MODISAlbedo = Data(1).MODISAlbedo;
            Data(e).MODIS_Albedo_File = Data(1).MODIS_Albedo_File;
        end
        
        % For now, we're just going to leave the cloud fractions as 0. We
        % could pull the geometric ones from MODIS, but that would only be
        % at 10a and 2p at best, and we can't get the radiance fraction
        % without radiative transfer as best as I can tell.
    end
    
    savename=[satellite,'_',retrieval,'_',year,month,day];
    if DEBUG_LEVEL > 0; fprintf('Saving %s\n',fullfile(local_sp_dir,savename)); end
    saveData(fullfile(local_sp_dir,savename), Data); % Saving must be handled as a separate function in a parfor loop because passing a variable name as a string upsets the parallelization monkey (it's not transparent).

end


end



function saveData(filename,Data)
save(filename,'Data','-v7.3')
end

function Data = addMODISAlbedo(Data, this_date, modis_mcd43_dir, alb_type, DEBUG_LEVEL)
E = JLLErrors;

%Convert the OMI date to a Julian calendar day
julian_day = modis_date_to_day(this_date);
yr = num2str(year(this_date));

%Add MODIS albedo info to the files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DEBUG_LEVEL>0; fprintf('\n Adding MODIS albedo information \n'); end
alb_dir = fullfile(modis_mcd43_dir,yr);

%Find the closest MCD file
in=[0 1 -1 2 -2 3 -3 4 -4 5 -5 6 -6 7 -7 8 -8 9 -9 10 -10 11 -11 12 -12 13 -13 14 -14 15 -15 16 -16 17 -17 18 -18 19 -19 20 -20 21 -21];
for ii=1:length(in);
    mcd_date = num2str(julian_day + in(ii),'%03g');
    if strcmpi(alb_type,'black-sky')
        alb_prefix = 'MCD43C3.A';
    elseif ismember(alb_type,{'poly','brdf'})
        alb_prefix = 'MCD43C1.A';
    else
        error('read_tempo_simulator:add_modis:bad_alb_type','Albedo prefix not recognized')
    end
    alb_filename = fullfile(alb_dir,[alb_prefix,yr,mcd_date,'*.hdf']);
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

% I think it causes a problem for parallel loops if a variable that *may*
% be used in the loop never gets defined before the loop starts, causing
% MATLAB to crash because it tries to send the variable just in case and 
% can't. That happens when which variable is used it selected by an if 
% statement like the versions of band3 here. Get around that by initializing
% them as empty variables.

band3 = [];
band3_iso = [];
band3_vol = [];
band3_geo = [];

mcd43_info = hdfinfo(fullfile(alb_dir,alb_files(1).name));
if strcmpi(alb_type,'black-sky')
    band3 = hdfread(hdf_dsetID(mcd43_info,1,1,'Albedo_BSA_Band3'));
    band3 = double(band3);
    band3 = flipud(band3);
    band3(band3==32767)=NaN; %JLL 11 Apr 2014: 32767 is the fill value for this data set; we will remove the NaNs further down
    band3 = band3 * 1e-3; %JLL 11-Apr-2014 Albedo matrix needs to have the scale factor applied
elseif ismember(alb_type, {'poly','brdf'})
    band3_iso = hdfread(hdf_dsetID(mcd43_info,1,1,'BRDF_Albedo_Parameter1_Band3'));
    band3_vol = hdfread(hdf_dsetID(mcd43_info,1,1,'BRDF_Albedo_Parameter2_Band3'));
    band3_geo = hdfread(hdf_dsetID(mcd43_info,1,1,'BRDF_Albedo_Parameter3_Band3'));

    band3_iso = flipud(double(band3_iso));
    band3_vol = flipud(double(band3_vol));
    band3_geo = flipud(double(band3_geo)); 

    % Remove fill values (as specified in the HDF attributes)
    band3_iso(band3_iso==32767)=nan;
    band3_vol(band3_vol==32767)=nan;
    band3_geo(band3_geo==32767)=nan;

    % Applied the scale factor (as specified in the HDF
    % attributes)
    band3_iso = band3_iso * 1e-3;
    band3_vol = band3_vol * 1e-3;
    band3_geo = band3_geo * 1e-3;
end

%MODIS albedo is given in 0.05 degree cells and a single file covers the
%full globe, so figure out the lat/lon of the middle of the grid cells as:
band3_lat=-90+0.05/2:0.05:90-0.05/2; band3_lats=band3_lat'; band3_lats=repmat(band3_lats,1,7200);
band3_lon=-180+0.05/2:0.05:180-0.05/2; band3_lons=repmat(band3_lon,3600,1);

%To speed up processing, restrict the MODIS albedo data to
%only the area we need to worry about.  This will
%significantly speed up the search for albedo values within
%each pixel.
lat_min=Data.Latcorn(:); lat_min(lat_min==0)=[]; lat_min=floor(min(lat_min));
lat_max=Data.Latcorn(:); lat_max(lat_max==0)=[]; lat_max=ceil(max(lat_max));
lon_min=Data.Loncorn(:); lon_min(lon_min==0)=[]; lon_min=floor(min(lon_min));
lon_max=Data.Loncorn(:); lon_max(lon_max==0)=[]; lon_max=ceil(max(lon_max));

in_lats = find(band3_lat>=lat_min & band3_lat<=lat_max);
in_lons = find(band3_lon>=lon_min & band3_lon<=lon_max);
if strcmpi(alb_type,'black-sky')
    band3=band3(in_lats,in_lons);
elseif ismember(alb_type,{'poly','brdf'})
    band3_iso = band3_iso(in_lats,in_lons); 
    band3_vol = band3_vol(in_lats,in_lons);
    band3_geo = band3_geo(in_lats,in_lons);
end
band3_lats=band3_lats(in_lats,in_lons);
band3_lons=band3_lons(in_lats,in_lons);
s=size(Data.Latitude);
c=numel(Data.Latitude);
MODISAlbedo=zeros(s);

if DEBUG_LEVEL > 3; MODISAlb_Ocn = zeros(s); end %JLL

%Now actually average the MODIS albedo for each OMI pixel
if DEBUG_LEVEL > 0; disp(' Averaging MODIS albedo to OMI pixels'); end
parfor k=1:c;
    t = getCurrentTask();
%    t.ID = 0;
    if DEBUG_LEVEL > 1 && mod(k,10000)==1; fprintf('Adding albedo data to pixel %d of %d\n',k,c); end
    if DEBUG_LEVEL > 2; tic; end
    x1 = Data.Loncorn(1,k);   y1 = Data.Latcorn(1,k);
    x2 = Data.Loncorn(2,k);   y2 = Data.Latcorn(2,k);
    x3 = Data.Loncorn(3,k);   y3 = Data.Latcorn(3,k);
    x4 = Data.Loncorn(4,k);   y4 = Data.Latcorn(4,k);
    
    
    xall=[x1;x2;x3;x4;x1];
    yall=[y1;y2;y3;y4;y1];
    
    xx_alb = inpolygon(band3_lats,band3_lons,yall,xall);
    
    if strcmpi(alb_type,'black-sky')
        if DEBUG_LEVEL>1 && mod(k,10000)==1; fprintf('W%d: Calculating average black-sky albedo\n',t.ID); end
        band3_vals=band3(xx_alb);  band3_zeros=band3_vals==0;
        band3_vals(band3_zeros)=NaN; band3_vals(isnan(band3_vals))=[];
        band3_avg=mean(band3_vals);
    elseif strcmpi(alb_type,'brdf')
        if DEBUG_LEVEL>1 && mod(k,10000)==1; fprintf('W%d: Calculating average BRDF albedo\n',t.ID); end
        band3_vals_brdf = modis_brdf_alb(band3_iso(xx_alb), band3_vol(xx_alb), band3_geo(xx_alb), Data.SolarZenithAngle(k), Data.ViewingZenithAngle, Data.RelativeAzimuthAngle(k));
        band3_avg = nanmean(band3_vals_brdf(band3_vals_brdf>0));
    elseif strcmpi(alb_type,'poly')    
        if DEBUG_LEVEL>1 && mod(k,10000)==1; fprintf('W%d: Calculating average polynomial albedo\n',t.ID); end
        band3_vals_poly = modis_brdf_alb_poly(band3_iso(xx_alb), band3_vol(xx_alb), band3_geo(xx_alb), Data.SolarZenithAngle(k));
        band3_avg = nanmean(band3_vals_poly);
    end
    %put in ocean surface albedo from LUT
    if isnan(band3_avg)==1;
        sza_vec = [5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 89];
        alb_vec = [0.038 0.038 0.039 0.039 0.040 0.042 0.044 0.046 0.051 0.058 0.068 0.082 0.101 0.125 0.149 0.158 0.123 0.073];
        alb = interp1(sza_vec,alb_vec,Data.SolarZenithAngle(k));
        band3_avg = alb;
        if DEBUG_LEVEL > 3; MODISAlb_Ocn(k) = 1; end
    end
    
    MODISAlbedo(k) = band3_avg;
    if DEBUG_LEVEL > 2; telap = toc; fprintf(' Time for MODIS alb --> pixel %u/%u = %g sec \n',k,c,telap); end
end

Data.MODISAlbedo = MODISAlbedo;
Data.MODIS_Albedo_File = fullfile(alb_dir,alb_files(1).name);
end

function Data = addGLOBETerpress(Data, globe_lon_matrix, globe_lat_matrix, terpres, DEBUG_LEVEL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Add GLOBE terrain pressure to the files%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if DEBUG_LEVEL > 0; fprintf('\n Adding GLOBE terrain data \n'); end

GLOBETerpres = zeros(size(Data.Latitude));

%GLOBE matrices are arrange s.t. terpres(1,1) is in the SW
%corner and terpres(end, end) is in the NE corner.

parfor k=1:numel(Data.Longitude)
    
    if DEBUG_LEVEL > 1 && mod(k,10000)==1; fprintf('Averaging GLOBE data to pixel %u of %u \n',k,numel(Data.Longitude)); end
    if DEBUG_LEVEL > 2; tic; end
    x1 = Data.Loncorn(1,k);   y1 = Data.Latcorn(1,k);
    x2 = Data.Loncorn(2,k);   y2 = Data.Latcorn(2,k);
    x3 = Data.Loncorn(3,k);   y3 = Data.Latcorn(3,k);
    x4 = Data.Loncorn(4,k);   y4 = Data.Latcorn(4,k);
    
    
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

Data.GLOBETerpres = GLOBETerpres;
end
