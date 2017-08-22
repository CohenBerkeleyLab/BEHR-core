function BEHR_main(varargin)
%Josh Laughner <joshlaugh5@gmail.com>
%Based on BEHR_nwus by Ashley Russell (02/09/2012)
%Takes "OMI_SP_yyyymmdd.m" files produced by read_omno2_v_aug2012.m as it's
%main input.

%****************************%
% CONSOLE OUTPUT LEVEL - 0 = none, 1 = minimal, 2 = all messages, 3 = times %
% Allows for quick control over the amount of output to the console.
% Choose a higher level to keep track of what the script is doing.
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

% Cleanup object will safely exit if there's a problem
if onCluster
    cleanupobj = onCleanup(@() mycleanup());
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DEPENDENCIES %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

%Add the 'Utils' folder and all subfolders to MATLAB's search path. Within
%the Git repository for BEHR, this is the /Utils folder.
mpath = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(mpath,'..','Utils')));


% Add the paths needed to run on the cluster. Modify these manually if
% needed.
addpath(genpath(BEHR_paths_old('classes')));
addpath(genpath(BEHR_paths_old('utils')));

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION AND INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addParameter('start', '');
p.addParameter('end', '');
p.addParameter('behr_mat_dir', '');
p.addParameter('sp_mat_dir', '');
p.addParameter('amf_tools_path', '');
p.addParameter('no2_profile_path', '');
p.addParameter('overwrite', false);

p.parse(varargin{:});
pout = p.Results;

date_start = pout.start;
date_end = pout.end;
behr_mat_dir = pout.behr_mat_dir;
sp_mat_dir = pout.sp_mat_dir;
amf_tools_path = pout.amf_tools_path;
no2_profile_path = pout.no2_profile_path;
overwrite = pout.overwrite;

%%% Validation %%%
date_start = validate_date(date_start);
date_end = validate_date(date_end);

if ~ischar(behr_mat_dir)
    E.badinput('Parameter "behr_mat_dir" must be a string');
elseif ~ischar(sp_mat_dir)
    E.badinput('Parameter "behr_mat_dir" must be a string');
elseif ~ischar(amf_tools_path)
    E.badinput('Parameter "amf_tools_path" must be a string');
elseif ~ischar(no2_profile_path)
    E.badinput('Parameter "no2_profile_path" must be a string');
elseif (~islogical(overwrite) && ~isnumeric(overwrite)) || ~isscalar(overwrite)
    E.badinput('Parameter "overwrite" must be a scalar logical or number')
end

%This is the directory where the final .mat file will be saved. This will
%need to be changed to match your machine and the files' location.
if isempty(behr_mat_dir)
    behr_mat_dir = BEHR_paths_old('behr_mat_dir');
end

%This is the directory where the "OMI_SP_*.mat" files are saved. This will
%need to be changed to match your machine and the files' location.
if isempty(sp_mat_dir)
    sp_mat_dir = BEHR_paths_old('sp_mat_dir');
end

%Add the path to the AMF_tools folder which contains rNmcTmp2.m,
%omiAmfAK2.m, integPr2.m and others.  In the Git repository for BEHR, this
%is the 'AMF_tools' folder.
if isempty(amf_tools_path)
    amf_tools_path = BEHR_paths_old('amf_tools_dir');
end

%This is the directory where the NO2 profiles are stored. This will
%need to be changed to match your machine and the files' location.
%no2_profile_path = '/Volumes/share/GROUP/SAT/BEHR/Monthly_NO2_Profiles';
if isempty(no2_profile_path)
    no2_profile_path = BEHR_paths_old('no2_profile_path');
end

% Verify the paths integrity.
nonexistant = {};

if ~exist(behr_mat_dir,'dir')
    nonexistant{end+1} = 'behr_mat_dir';
end
if ~exist(sp_mat_dir,'dir')
    nonexistant{end+1} = 'sp_mat_dir';
end
if ~exist(amf_tools_path,'dir')
    nonexistant{end+1} = 'amf_tools_path';
end
if ~exist(no2_profile_path,'dir')
    nonexistant{end+1} = 'no2_profile_path';
end

if numel(nonexistant)>0
    string_spec = [repmat('\n\t%s',1,numel(nonexistant)),'\n\n'];
    msg = sprintf('The following paths are not valid: %s Please double check them in the run file',string_spec);
    error(E.callError('bad_cluster_path',sprintf(msg,nonexistant{:})));
end

%Store paths to relevant files
addpath(amf_tools_path)
fileTmp = fullfile(amf_tools_path,'nmcTmpYr.txt');
fileDamf = fullfile(amf_tools_path,'damf.txt');
%****************************%

%****************************%
%Process all files between these dates, in yyyy/mm/dd format
%****************************%
if isempty(date_start) || isempty(date_end)
    date_start='2013/08/01';
    date_end='2013/08/06';
    fprintf('BEHR_main: Used hard-coded start and end dates (%s, %s)\n', date_start, date_end);
end
%****************************%

%****************************%
% Which cloud product to use to calculate the AMF: OMI or MODIS
%****************************%
cloud_amf = 'omi';
%****************************%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a parallel pool if one doesn't exist and we are on a cluster
if onCluster && isempty(gcp('nocreate'))
    parpool(numThreads);
end

if onCluster
    n_workers=numThreads;
else
    % Running a parfor loop with 0 workers makes it run in serial mode,
    % which means it doesn't waste time sending data to and from the
    % workers
    n_workers=0;
end

datenums = datenum(date_start):datenum(date_end);
parfor(j=1:length(datenums), n_workers)
%for j=1:length(datenums)
    savename = behr_filename(datenums(j));
    
    if exist(fullfile(behr_mat_dir, savename),'file') && ~overwrite
        fprintf('%s already exists, skipping\n', savename);
        continue
    end
    
    
    if DEBUG_LEVEL > 0
        fprintf('Processing data for %s', datestr(datenums(j)));
    end
    sp_mat_name = sp_savename(datenums(j));
    
    if DEBUG_LEVEL > 1
        fprintf('Looking for SP file %s ...', fullfile(sp_mat_dir,sp_mat_name));
    end
    
    if ~exist(fullfile(sp_mat_dir,sp_mat_name),'file')
        if DEBUG_LEVEL > 0; disp('No SP file exists for given day'); end
        continue
    end
    if DEBUG_LEVEL > 1; fprintf('\t ...Found.\n'); end
    S=load(fullfile(sp_mat_dir,sp_mat_name)); %JLL 17 Mar 2014: Will load the variable 'Data' into the workspace
    Data=S.Data;
    
    if ~exist('profile_file','var') || ~strcmp(profile_file(2:3),month_str)
        month_str=datestr(datenums(j),'mm');
        profile_file=['m',month_str,'_NO2_profile'];
        
        if DEBUG_LEVEL > 1; disp(['Loading ',fullfile(no2_profile_path,profile_file)]); end
        S=load(fullfile(no2_profile_path,profile_file));
        PROFILE = S.PROFILE;
    end
    for d=1:length(Data);
        % Data is initialized in read_omno2_v_aug2012 with a single 0
        % in the Longitude field.  Since points outside the lat/lons of
        % interest are removed completely, we should also check if all
        % points are gone.
        if numel(Data(d).Longitude)==1 || isempty(Data(d).Longitude);
            if DEBUG_LEVEL > 1; fprintf('  Note: Data(%u) is empty\n',d); end
            continue %JLL 17 Mar 2014: Skip doing anything if there's really no information in this data
        end
        if DEBUG_LEVEL>0; fprintf('  Swath %u of %s \n',d,datestr(datenums(j))); end
        c=numel(Data(d).Longitude);
        
        %JLL 17 Mar 2014: Load some of the variables from 'Data' to
        %make referencing them less cumbersome. Also convert some
        %to column vectors to work with rNmcTmp2 and rDamf2
        lon = Data(d).Longitude(:);
        lat = Data(d).Latitude(:);
        loncorns = Data(d).FoV75CornerLongitude;
        latcorns = Data(d).FoV75CornerLatitude;
        
        sza = Data(d).SolarZenithAngle(:);
        vza = Data(d).ViewingZenithAngle(:);
        phi = Data(d).RelativeAzimuthAngle(:);
        
        mon = month(datenums(j)) * ones(size(Data(d).Latitude(:)));
        pressure = behr_pres_levels();
        if DEBUG_LEVEL > 1; fprintf('   Interpolating temperature data\n'); end
        temperature = rNmcTmp2(fileTmp, pressure, lon, lat, mon); %JLL 17 Mar 2014: Interpolates temperature values to the pressures and lat/lon coordinates desired
        
        surfPres = Data(d).GLOBETerpres(:);
        albedo = Data(d).MODISAlbedo(:);
        
        surfPres(surfPres>=1013)=1013; %JLL 17 Mar 2014: Clamp surface pressure to sea level or less.
        cldPres = Data(d).CloudPressure(:);
        cldPres(cldPres>=1013)=1013; % JLL 13 May 2016: Also clamp cloud pressure. Whenever this is >1013, the AMF becomes a NaN because the lookup table cannot handle "surface" pressure >1013
        
        if DEBUG_LEVEL > 1; fprintf('   Calculating clear and cloudy AMFs\n'); end
        dAmfClr = rDamf2(fileDamf, pressure, sza, vza, phi, albedo, surfPres); %JLL 18 Mar 2014: Interpolate the values in dAmf to the albedo and other conditions input
        cloudalbedo=0.8*ones(size(Data(d).CloudFraction(:))); %JLL 18 Mar 2014: Assume that any cloud has an albedo of 0.8
        dAmfCld = rDamf2(fileDamf, pressure, sza, vza, phi, cloudalbedo, cldPres); %JLL 18 Mar 2014: Interpolate dAmf again, this time taking the cloud top and albedo as the bottom pressure
        
        if DEBUG_LEVEL > 1; fprintf('   Reading NO2 profiles\n'); end
        [no2_bins] = rProfile_US(PROFILE, loncorns, latcorns, c); %JLL 18 Mar 2014: Bins the NO2 profiles to the OMI pixels; the profiles are averaged over the pixel
        no2_bins = reshape(no2_bins,length(pressure),size(vza,1),size(vza,2));
        no2Profile = no2_bins ./ (10^6); % NO2 from WRF is in ppm in these files
        prof_i = isnan(squeeze(no2Profile(1, :, :))); %JLL 18 Mar 2014: prof_i is a matrix of 0 or 1s that is 1 wherever the bottom of NO2 profile is NaN
        
        pTerr = surfPres;
        pCld = cldPres;
        if strcmpi(cloud_amf,'omi')
            cldFrac = Data(d).CloudFraction(:);
        else
            cldFrac = Data(d).MODISCloud(:);
        end
        
        cldRadFrac = Data(d).CloudRadianceFraction(:);
        
        
        if DEBUG_LEVEL > 1; disp('   Calculating BEHR AMF'); end
        [amf, amfVis, ~, ~, scattering_weights, avg_kernels, no2_prof_interp, sw_plevels] = omiAmfAK2(pTerr, pCld, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile); %JLl 18 Mar 2014: The meat and potatoes of BEHR, where the TOMRAD AMF is adjusted to use the GLOBE pressure and MODIS cloud fraction
        amf(prof_i==1)=NaN;
        amfVis(prof_i==1)=NaN;
        scattering_weights(:,prof_i==1)=NaN;
        avg_kernels(:,prof_i==1)=NaN;
        sw_plevels(:,prof_i==1)=NaN;
        no2_prof_interp(:,prof_i==1)=NaN;
        
        sz = size(Data(d).Longitude);
        len_vecs = size(scattering_weights,1);  % JLL 26 May 2015 - find out how many pressure levels there are. Will often be 30, but might change.
        % Need this to properly reshape the scattering weights, AKs, pressure levels, and (soon) profiles
        
        Data(d).BEHRAMFTrop = reshape(amf,sz); %JLL 18 Mar 2014: Save the resulting AMF of the pixel
        Data(d).BEHRAMFTropVisOnly = reshape(amfVis,sz);
        Data(d).BEHRScatteringWeights = reshape(scattering_weights, [len_vecs, sz]);
        Data(d).BEHRAvgKernels = reshape(avg_kernels, [len_vecs, sz]);
        Data(d).BEHRNO2apriori = reshape(no2_prof_interp, [len_vecs, sz]);
        Data(d).BEHRPressureLevels = reshape(sw_plevels, [len_vecs, sz]);
        
    end
    
    b=length(Data);
    for z=1:b;
        if ~isfield(Data,'BEHRAMFTrop') || isempty(Data(z).BEHRAMFTrop)
            continue
        end
        Data(z).BEHRColumnAmountNO2Trop=Data(z).ColumnAmountNO2Trop.*Data(z).AmfTrop./Data(z).BEHRAMFTrop;
        Data(z).BEHRColumnAmountNO2TropVisOnly=Data(z).ColumnAmountNO2Trop.*Data(z).AmfTrop./Data(z).BEHRAMFTropVisOnly;
        % make sure fill values in the original column or AMF are
        % fill values in BEHR.
        Data(z).BEHRColumnAmountNO2Trop(Data(z).ColumnAmountNO2Trop < -1e29 | Data(z).AmfTrop < -30000) = nan;
        Data(z).BEHRColumnAmountNO2TropVisOnly(Data(z).ColumnAmountNO2Trop < -1e29 | Data(z).AmfTrop < -30000) = nan;
        if DEBUG_LEVEL > 0; fprintf('   BEHR [NO2] stored for swath %u\n',z); end
        
    end
    
    
    addpath('/Users/Josh/Documents/MATLAB/BEHR/Utils/m_map'); %JLL 18 Mar 2014: Adds the path to the m_map toolkit, needed for hdf_quadrangle_5km_new
    
    %*********************************%
    %JLL 19 Mar 2014: These will be used to define the quadrangles -
    %the quads will be smaller than the OMI pixel, and multiple quads
    %will take on the value for the same (closest) OMI pixel.  By
    %keeping the quads' centers the same over all retrievals you wish
    %to average, this will allow easier averaging over multiple OMI
    %swaths. This is a form of oversampling.
    %*********************************%
    lonmin = -125;  lonmax = -65;
    latmin = 25;   latmax = 50;
    resolution = 0.05; resolution2 = 0.05;
    %*********************************%
    %
    if lonmin > lonmax %Just in case I enter something backwards...
        error(E.badinput('Lonmin is greater than lonmax'))
    elseif latmin > latmax
        error(E.badinput('Latmin is greater than latmax'))
    end
    
    %*********************************%
    %JLL 19 Mar 2014: Save all relevant values produced by add2grid to
    %a new structure called 'OMI'
    %*********************************%
    
    if DEBUG_LEVEL > 0; disp('  Preparing OMI structure'); end
    s=numel(Data);
    
    % Prepare the OMI data structure which will receive the gridded
    % data - this will be passed to the gridding functions to keep the
    % field names in the right order.
    OMI=struct('Date','','Longitude', [], 'Latitude', [], 'Time', [], 'ViewingZenithAngle', [], 'SolarZenithAngle', [], 'ViewingAzimuthAngle', [], 'SolarAzimuthAngle', [],...
        'RelativeAzimuthAngle', [], 'AmfStrat', [], 'AmfTrop',[], 'CloudFraction', [], 'CloudRadianceFraction', [], 'CloudPressure', [], 'ColumnAmountNO2', [],...
        'SlantColumnAmountNO2', [], 'ColumnAmountNO2Trop', [], 'ColumnAmountNO2TropStd',[],'ColumnAmountNO2Strat',[],'TerrainHeight', [], 'TerrainPressure', [], 'TerrainReflectivity', [], 'VcdQualityFlags',{{}},...
        'MODISCloud', [], 'MODISAlbedo', [], 'GLOBETerpres', [], 'XTrackQualityFlags', {{}}, 'Row', [], 'Swath', [], 'TropopausePressure', [], 'BEHRColumnAmountNO2Trop',[],...
        'BEHRAMFTrop', [], 'BEHRColumnAmountNO2TropVisOnly', [], 'BEHRAMFTropVisOnly', [], 'Count', [], 'Area', [], 'Areaweight', [], 'MapData', struct);
    % Matlab treats structures as matrices, so we can duplicate our
    % structure to have the required number of entries just like a
    % matrix.
    OMI = repmat(OMI,1,s);
    hh=0;
    for d=1:s;
        if Data(d).ViewingZenithAngle==0;
        elseif numel(Data(d).ViewingZenithAngle)==1;
            continue
        else
            if DEBUG_LEVEL > 1; fprintf('   Gridding data for swath %u\n',d); end
            hh=hh+1;
            OMI(hh) = add2grid_BEHR(Data(d),OMI(hh),resolution,resolution2,[lonmin, lonmax],[latmin, latmax]); %JLL 20 Mar 2014: Superimpose data to a grid determined by lat & lon min/max and resolution above. Default resolution is 0.05 degree
        end
    end
    
    % Clean up any unused elements in OMI
    OMI(hh+1:end) = [];
    
    savename = behr_filename(datenums(j));
    if DEBUG_LEVEL > 0; disp(['   Saving data as',fullfile(behr_mat_dir,savename)]); end
    saveData(fullfile(behr_mat_dir,savename),Data,OMI)
    
end
end

function saveData(filename,Data,OMI)
save(filename,'OMI','Data')
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
