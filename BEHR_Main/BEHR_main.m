function BEHR_main(varargin)
% BEHR_MAIN: primary BEHR algorithm
%
%   This function is the primary BEHR algorithm, it takes the OMI, MODIS,
%   and GLOBE data read in by read_main.m and uses it to
%   recalculated the BEHR AMFs and VCDs. There are a number of input
%   parameters that control it's operation; the defaults are set such that
%   it should run if you simply execute this script, but in most cases you
%   will want to change at least the start and end dates.
%
%   Parameters:
%       'start' - the first date to process as a date number or date string
%       that Matlab recognized implicitly. If not given, defaults to
%       2005-01-01.
%
%       'end' - the last date to process, same format requirements as
%       starting date. If not given, defaults to today.
%
%       'behr_mat_dir' - the directory that the final .mat file should be
%       saved in. If not given, it defaults to
%       fullfile(behr_paths.behr_mat_dir, lower(region), lower(prof_mode)).
%
%       'sp_mat_dir' - the directory that the .mat files resulting from
%       read_main.m are stored in. If not given, it defaults to
%       fullfile(behr_paths.sp_mat_dir, lower(region))
%
%       'amf_tools_path' - the directory that contains the files
%       nmcTmpYr.txt and damf.txt. If not given, defaults to the path
%       stored in behr_paths.m
%
%       'no2_profile_path' - the directory to look for WRF output files in.
%       If not given, or given as an empty string, this is determined
%       automatically.
%
%       'overwrite' - a boolean that controls whether existing files in the
%       behr_mat_dir should be overwritten or not. Defaults to false (if a
%       file exists for a given day, that day will not be reprocessed).
%
%       'profile_mode' - must be the string 'daily' or 'monthly' (defaults
%       to 'monthly'). Controls whether daily or monthly profiles will be
%       used, which also controls whether rProfile_WRF.m looks for files
%       named 'WRF_BEHR_monthly_yyyy-mm.nc' (monthly) or
%       'wrfout_*_yyyy-mm-dd_hh-00-00' (daily).
%
%       'use_psm_gridding' - if false (default), uses CVM gridding for all
%       fields. If true, then NO2 fields will be gridded using the PSM
%       method (specifically, fields specified as psm_gridded_vars in
%       BEHR_publishing_gridded_fields will be gridded by PSM).
%
%       'DEBUG_LEVEL' - level of progress messaged printed to the console.
%       0 = none, 1 = minimal, 2 = all, 3 = processing times are added.
%       Default is 2.
%
%Josh Laughner <joshlaugh5@gmail.com>

%Based on BEHR_nwus by Ashley Russell (02/09/2012)


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION AND INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% You may change the default values here if you just want to be able to
% click 'Run', but it's usually better to pass these as parameters.
p = inputParser;
p.addParameter('start', '2005-01-01');
p.addParameter('end', today);
p.addParameter('behr_mat_dir', '');
p.addParameter('sp_mat_dir', '');
p.addParameter('amf_tools_path', behr_paths.amf_tools_dir);
p.addParameter('no2_profile_path', '');
p.addParameter('region', 'us');
p.addParameter('overwrite', false);
p.addParameter('profile_mode', 'monthly');
p.addParameter('use_psm_gridding', false);
p.addParameter('DEBUG_LEVEL', 2);

p.parse(varargin{:});
pout = p.Results;

date_start = pout.start;
date_end = pout.end;
behr_mat_dir = pout.behr_mat_dir;
sp_mat_dir = pout.sp_mat_dir;
no2_profile_path = pout.no2_profile_path;
region = pout.region;
overwrite = pout.overwrite;
prof_mode = pout.profile_mode;
use_psm = pout.use_psm_gridding;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

%%% Validation %%%
allowed_prof_modes = {'daily','monthly'};

date_start = validate_date(date_start);
date_end = validate_date(date_end);

if ~ischar(behr_mat_dir)
    E.badinput('Parameter "behr_mat_dir" must be a string');
elseif ~ischar(sp_mat_dir)
    E.badinput('Parameter "sp_mat_dir" must be a string');
elseif ~ischar(no2_profile_path)
    E.badinput('Parameter "no2_profile_path" must be a string');
elseif (~islogical(overwrite) && ~isnumeric(overwrite)) || ~isscalar(overwrite)
    E.badinput('Parameter "overwrite" must be a scalar logical or number')
elseif ~ismember(prof_mode,allowed_prof_modes)
   	E.badinput('prof_mode (if given) must be one of %s', strjoin(allowed_prof_modes,', '));
elseif ~isscalar(use_psm) || (~islogical(use_psm) && ~isnumeric(use_psm))
    E.badinput('use_psm_gridding must be a scalar logical or number')
end

% If using the default SP file directory, look in the right region
% subfolder.
if isempty(sp_mat_dir)
    sp_mat_dir = behr_paths.SPMatSubdir(region);
end

% Set behr_mat_dir to the daily or monthly directory, with region
% subdirectory, if using the default path
if isempty(behr_mat_dir)
    behr_mat_dir = behr_paths.BEHRMatSubdir(region, prof_mode);
end

% Bring the AMF tools path (where damf.txt is) into this workspace
amf_tools_path = behr_paths.amf_tools_dir;

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
if ~isempty(no2_profile_path) && ~exist(no2_profile_path,'dir')
    nonexistant{end+1} = 'no2_profile_path';
end

if numel(nonexistant)>0
    string_spec = [repmat('\n\t%s',1,numel(nonexistant)),'\n\n'];
    msg = sprintf('The following paths are not valid: %s Please double check them in the run file',string_spec);
    error(E.callError('bad_cluster_path',sprintf(msg,nonexistant{:})));
end

%Store paths to relevant files
fileDamf = fullfile(amf_tools_path,'damf.txt');

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
    fprintf('Assuming onCluster is false\n');
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


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a parallel pool if one doesn't exist and we are on a cluster
if onCluster    
    if isempty(gcp('nocreate'))
        parpool(numThreads);
    end    
    n_workers = Inf;
else
    n_workers = 0;
end

if onCluster
    n_workers=numThreads;
else
    % Running a parfor loop with 0 workers makes it run in serial mode,
    % which means it doesn't waste time sending data to and from the
    % workers
    n_workers=0;
end

core_githead = git_head_hash(behr_paths.behr_core);
behrutils_githead = git_head_hash(behr_paths.behr_utils);
genutils_githead = git_head_hash(behr_paths.utils);
psm_githead = git_head_hash(behr_paths.psm_dir);
imatpy_githead = git_head_hash(behr_paths.python_interface);
wrfutils_githead = git_head_hash(behr_paths.wrf_utils);

datenums = datenum(date_start):datenum(date_end);
%parfor(j=1:length(datenums), n_workers)
for j=1:length(datenums)
    savename = behr_filename(datenums(j), prof_mode, region);
    
    if exist(fullfile(behr_mat_dir, savename),'file') && ~overwrite
        fprintf('%s already exists, skipping\n', savename);
        continue
    end
    
    
    if DEBUG_LEVEL > 0
        fprintf('Processing data for %s\n', datestr(datenums(j)));
    end
    sp_mat_name = sp_savename(datenums(j), region);
    
    if DEBUG_LEVEL > 1
        fprintf('Looking for SP file %s ...', fullfile(sp_mat_dir,sp_mat_name));
    end
    
    if ~exist(fullfile(sp_mat_dir,sp_mat_name),'file')
        if DEBUG_LEVEL > 0; disp('No SP file exists for given day'); end
        continue
    end
    
    if DEBUG_LEVEL > 1; fprintf('\t ...Found.\n'); end
    S=load(fullfile(sp_mat_dir,sp_mat_name));
    Data=S.Data;
    
    if isempty(Data)
        % If there is no data read in from read_main.m, then there's
        % nothing to do for this day.
        continue
    end
    
    % Double check that the loaded SP file is for the same region as we're
    % trying to process
    if ~strcmpi(Data(1).BEHRRegion, region)
        E.callError('behr_region', 'Somehow I loaded a file with a different region specified in Data.BEHRRegion (%s) than I am trying to process (%s)', Data(1).BEHRRegion, region)
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE OUR AMFS %
    %%%%%%%%%%%%%%%%%%%%%%
    
    for d=1:length(Data);
        % Data is initialized in read_main with a single 0 in the Longitude
        % field.  Since points outside the lat/lons of interest are removed
        % completely, we should also check if all points are gone.
        if numel(Data(d).Longitude)==1 || isempty(Data(d).Longitude);
            if DEBUG_LEVEL > 1; fprintf('  Note: Data(%u) is empty\n',d); end
            continue %JLL 17 Mar 2014: Skip doing anything if there's really no information in this data
        end
        if DEBUG_LEVEL>0; fprintf('  Swath %u of %s \n',d,datestr(datenums(j))); end
        c=numel(Data(d).Longitude);
        
        %JLL 17 Mar 2014: Load some of the variables from 'Data' to
        %make referencing them less cumbersome. Also convert some
        %to column vectors to work with rNmcTmp2 and rDamf2
        loncorns = Data(d).FoV75CornerLongitude;
        latcorns = Data(d).FoV75CornerLatitude;
        time = Data(d).Time;       
 
        sza = Data(d).SolarZenithAngle;
        vza = Data(d).ViewingZenithAngle;
        phi = Data(d).RelativeAzimuthAngle;
        surfPres = Data(d).GLOBETerpres;
        albedo = Data(d).MODISAlbedo;
        cldFrac = Data(d).CloudFraction;
        cldRadFrac = Data(d).CloudRadianceFraction;
        
        pressure = behr_pres_levels();
        
        surfPres(surfPres>=1013)=1013; %JLL 17 Mar 2014: Clamp surface pressure to sea level or less.
        cldPres = Data(d).CloudPressure;
        cldPres(cldPres>=1013)=1013; % JLL 13 May 2016: Also clamp cloud pressure. Whenever this is >1013, the AMF becomes a NaN because the lookup table cannot handle "surface" pressure >1013
        
        if DEBUG_LEVEL > 1; fprintf('   Reading NO2 and temperature profiles\n'); end
        [no2Profile, temperature, wrf_profile_file, wrf_pres_mode, wrf_temp_mode] = rProfile_WRF(datenums(j), prof_mode, region, loncorns, latcorns, time, surfPres, pressure, no2_profile_path); %JLL 18 Mar 2014: Bins the NO2 profiles to the OMI pixels; the profiles are averaged over the pixel
        bad_profs = squeeze(all(isnan(no2Profile),1));
        
        if DEBUG_LEVEL > 1; fprintf('   Calculating clear and cloudy AMFs\n'); end
        dAmfClr = rDamf2(fileDamf, pressure, sza, vza, phi, albedo, surfPres); %JLL 18 Mar 2014: Interpolate the values in dAmf to the albedo and other conditions input
        cloudalbedo=0.8*ones(size(Data(d).CloudFraction)); %JLL 18 Mar 2014: Assume that any cloud has an albedo of 0.8
        dAmfCld = rDamf2(fileDamf, pressure, sza, vza, phi, cloudalbedo, cldPres); %JLL 18 Mar 2014: Interpolate dAmf again, this time taking the cloud top and albedo as the bottom pressure

        if DEBUG_LEVEL > 1; disp('   Calculating BEHR AMF'); end
        [amf, amfVis, ~, ~, scattering_weights_clear, scattering_weights_cloudy, avg_kernels, no2_prof_interp, sw_plevels] = omiAmfAK2(surfPres, cldPres, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile); %JLl 18 Mar 2014: The meat and potatoes of BEHR, where the TOMRAD AMF is adjusted to use the GLOBE pressure and MODIS cloud fraction
        amf(bad_profs)=NaN;
        amfVis(bad_profs)=NaN;
        scattering_weights_clear(:,bad_profs)=NaN;
        scattering_weights_cloudy(:,bad_profs)=NaN;
        avg_kernels(:,bad_profs)=NaN;
        sw_plevels(:,bad_profs)=NaN;
        no2_prof_interp(:,bad_profs)=NaN;
        
        sz = size(Data(d).Longitude);
        len_vecs = size(scattering_weights_clear,1);  % JLL 26 May 2015 - find out how many pressure levels there are. Will often be 30, but might change.
        % Need this to properly reshape the scattering weights, AKs, pressure levels, and (soon) profiles
        
        Data(d).BEHRAMFTrop = reshape(amf,sz); %JLL 18 Mar 2014: Save the resulting AMF of the pixel
        Data(d).BEHRAMFTropVisOnly = reshape(amfVis,sz);
        Data(d).BEHRScatteringWeightsClear = reshape(scattering_weights_clear, [len_vecs, sz]);
        Data(d).BEHRScatteringWeightsCloudy = reshape(scattering_weights_cloudy, [len_vecs, sz]);
        Data(d).BEHRAvgKernels = reshape(avg_kernels, [len_vecs, sz]);
        Data(d).BEHRNO2apriori = reshape(no2_prof_interp, [len_vecs, sz]);
        Data(d).BEHRWRFFile = wrf_profile_file;
        Data(d).BEHRWRFPressureMode = wrf_pres_mode;
        Data(d).BEHRWRFTemperatureMode = wrf_temp_mode;
        Data(d).BEHRProfileMode = prof_mode;
        Data(d).BEHRPressureLevels = reshape(sw_plevels, [len_vecs, sz]);
        Data(d).BEHRQualityFlags = behr_quality_flags(Data(d)); 
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATE VCDS FROM NASA SCDS AND OUR AMFS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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
        
        Data(z).GitHead_Core_Main = core_githead;
        Data(z).GitHead_BEHRUtils_Main = behrutils_githead;
        Data(z).GitHead_GenUtils_Main = genutils_githead;
        Data(z).GitHead_PSM_Main = psm_githead;
        Data(z).GitHead_MatPyInt_Main = imatpy_githead;
        Data(z).GitHead_WRFUtils_Main = wrfutils_githead;
    end
    
    
    %%%%%%%%%%%%%%%%%
    % GRIDDING DATA %
    %%%%%%%%%%%%%%%%%
    
    OMI = psm_wrapper(Data, Data(1).Grid, 'only_cvm', ~use_psm, 'DEBUG_LEVEL', DEBUG_LEVEL);
    
    %%%%%%%%%%%%%
    % SAVE FILE %
    %%%%%%%%%%%%%
    
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
