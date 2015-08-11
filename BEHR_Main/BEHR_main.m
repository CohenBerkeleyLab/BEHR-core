function BEHR_main(date_start, date_end)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FILE LOCATIONS %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The location of the directories to read or save data to.  If onCluster is
% true, these will need to be set in the runscript - I figured this would
% be easier than setting them as shell environmental variables and using
% getenv - JLL 15 Jan 2015

if onCluster
    global behr_mat_dir;
    global sp_mat_dir;
    global amf_tools_path;
    global no2_profile_path;
    
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
    
    
else
    %This is the directory where the final .mat file will be saved. This will
    %need to be changed to match your machine and the files' location.
    behr_mat_dir = '/Volumes/share-sat/SAT/BEHR/BEHR_Files_2014';
    
    %This is the directory where the "OMI_SP_*.mat" files are saved. This will
    %need to be changed to match your machine and the files' location.
    sp_mat_dir = '/Volumes/share-sat/SAT/BEHR/SP_Files_2014';
    
    %Add the path to the AMF_tools folder which contains rNmcTmp2.m,
    %omiAmfAK2.m, integPr2.m and others.  In the Git repository for BEHR, this
    %is the 'AMF_tools' folder.
    amf_tools_path = '/Users/Josh/Documents/MATLAB/BEHR/AMF_tools';

    %This is the directory where the NO2 profiles are stored. This will
    %need to be changed to match your machine and the files' location.
    %no2_profile_path = '/Volumes/share/GROUP/SAT/BEHR/Monthly_NO2_Profiles';
    no2_profile_path = '/Volumes/share-sat/SAT/BEHR/Monthly_NO2_Profiles';
end

%Store paths to relevant files
addpath(amf_tools_path)
fileTmp = fullfile(amf_tools_path,'nmcTmpYr.txt');
fileDamf = fullfile(amf_tools_path,'damf.txt');
fileNO2 = fullfile(amf_tools_path,'PRFTAV.txt');
%****************************%

%****************************%
%Process all files between these dates, in yyyy/mm/dd format
%****************************%
if nargin < 2
    date_start='2013/08/01';
    date_end='2013/08/06';
    fprintf('BEHR_main: Used hard-coded start and end dates\n');
end
%****************************%

%These will be included in the file name
%****************************%
satellite='OMI';
retrieval='BEHR_omiCloudAMF';
%****************************%

%****************************%
% Which cloud product to use to calculate the AMF: OMI or MODIS
%****************************%
cloud_amf = 'omi';
cloud_rad_amf = 'omi';
%****************************%


% Check that all directories given ARE directories
if ~exist(behr_mat_dir,'dir')
    E.filenotfound(behr_mat_dir)
elseif ~exist(sp_mat_dir,'dir')
    E.filenotfound(sp_mat_dir)
elseif ~exist(amf_tools_path,'dir')
    E.filenotfound(amf_tools_path)
elseif ~exist(no2_profile_path,'dir')
    E.filenotfound(no2_profile_path)
end

%Find the last file completed and set the start date to the next day.  This
%will allow the process to be stopped and started with minimum
%intervention.
file_prefix = [satellite,'_',retrieval,'_']; l = length(file_prefix);
last_file=dir(fullfile(behr_mat_dir,sprintf('%s*.mat',file_prefix)));

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

tic
% Create a parallel pool if one doesn't exist and we are on a cluster
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
    if DEBUG_LEVEL > 0; disp(['Processing data for ', date]); end
    
    filename = ['OMI_SP_',year,month,day,'.mat'];

    if DEBUG_LEVEL > 1; disp(['Looking for SP file ',fullfile(sp_mat_dir,filename),'...']); end %#ok<PFGV> % The concern with using global variables in a parfor is that changes aren't synchronized.  Since I'm not changing them, it doesn't matter.
    if isequal(exist(fullfile(sp_mat_dir,filename),'file'),0)
        if DEBUG_LEVEL > 0; disp('No SP file exists for given day'); end
        continue
    else
        if DEBUG_LEVEL > 1; fprintf('\t ...Found.\n'); end
        S=load(fullfile(sp_mat_dir,filename)); %JLL 17 Mar 2014: Will load the variable 'Data' into the workspace
        Data=S.Data;
        
        if exist('profile_file','file')==1 && strcmp(profile_file(2:3),month)==1; %JLL 20 Mar 2014: 
        else
            profile_file=['m',month,'_NO2_profile'];
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
            else
                if DEBUG_LEVEL>0; fprintf('  Swath %u of %s \n',d,date); end
                c=numel(Data(d).Longitude);
                
                Data(d).MODISAlbedo(isnan(Data(d).MODISAlbedo)==1)=0; %JLL 17 Mar 2014: replace NaNs with fill values
                Data(d).GLOBETerpres(isnan(Data(d).GLOBETerpres)==1)=1013.0000;
                
                %JLL 17 Mar 2014: Load some of the variables from 'Data' to
                %make referencing them less cumbersome. Also convert some
                %to column vectors to work with rNmcTmp2 and rDamf2
                lon = Data(d).Longitude(:);
                lat = Data(d).Latitude(:);
                loncorns=Data(d).Loncorn;
                latcorns=Data(d).Latcorn;
                
                sza = Data(d).SolarZenithAngle(:);
                vza = Data(d).ViewingZenithAngle(:);
                phi = Data(d).RelativeAzimuthAngle(:);
                
                mon = str2double(month)*ones(size(Data(d).Latitude(:)));
                pressure = [1020 1015 1010 1005 1000 990 980 970 960 945 925 900 875 850 825 800 770 740 700 660 610 560 500 450 400 350 280 200];% 100 50 20 5];
                if DEBUG_LEVEL > 1; disp('   Interpolating temperature data'); end
                [temperature, tmpSAVE] = rNmcTmp2(fileTmp, pressure, lon, lat, mon); %JLL 17 Mar 2014: Interpolates temperature values to the pressures and lat/lon coordinates desired
                
                surfPres = Data(d).GLOBETerpres(:);
                albedo = Data(d).MODISAlbedo(:);
                
                surfPres(surfPres>=1013)=1013; %JLL 17 Mar 2014: Clamp surface pressure to sea level or less.
                cldPres = Data(d).CloudPressure(:);
                
                if DEBUG_LEVEL > 1; disp('   Calculating clear and cloudy AMFs'); end
                dAmfClr = rDamf2(fileDamf, pressure, sza, vza, phi, albedo, surfPres); %JLL 18 Mar 2014: Interpolate the values in dAmf to the albedo and other conditions input
                cloudalbedo=0.8*ones(size(Data(d).CloudFraction(:))); %JLL 18 Mar 2014: Assume that any cloud has an albedo of 0.8
                dAmfCld = rDamf2(fileDamf, pressure, sza, vza, phi, cloudalbedo, cldPres); %JLL 18 Mar 2014: Interpolate dAmf again, this time taking the cloud top and albedo as the bottom pressure
                
                if DEBUG_LEVEL > 1; disp('   Reading NO2 profiles'); end
                [no2_bins] = rProfile_US(PROFILE, loncorns, latcorns, c); %JLL 18 Mar 2014: Bins the NO2 profiles to the OMI pixels; the profiles are averaged over the pixel
                no2_bins = reshape(no2_bins,length(pressure),size(vza,1),size(vza,2));
                no2Profile1 = no2_bins ./ (10^6); % NO2 from WRF is in ppm in these files
                prof_i = zeros(size(Data(d).Latitude)); prof_i(isnan(squeeze(no2Profile1(1,:,:)))==1)=1; %JLL 18 Mar 2014: prof_i is a matrix of 0 or 1s that is 1 wherever the bottom of NO2 profile is NaN
                
                no2Profile2 = no2Profile1;
                pTerr = surfPres;
                pCld = cldPres;
                if strcmpi(cloud_amf,'omi')
                    cldFrac = Data(d).CloudFraction(:); 
                else
                    cldFrac = Data(d).MODISCloud(:);
                end
                
                cldRadFrac = Data(d).CloudRadianceFraction(:);
                
                
                if DEBUG_LEVEL > 1; disp('   Calculating BEHR AMF'); end
                noGhost=0; ak=1;
                [amf, ~, ~, scattering_weights, avg_kernels, no2_prof_interp, sw_plevels, ghost_fraction] = omiAmfAK2(pTerr, pCld, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile1, no2Profile2, noGhost, ak); %JLl 18 Mar 2014: The meat and potatoes of BEHR, where the TOMRAD AMF is adjusted to use the GLOBE pressure and MODIS cloud fraction
                amf(prof_i==1)=NaN;
                scattering_weights(:,prof_i==1)=NaN;
                avg_kernels(:,prof_i==1)=NaN;
                sw_plevels(:,prof_i==1)=NaN;
                no2_prof_interp(:,prof_i==1)=NaN;
                ghost_fraction(prof_i==1)=NaN;
                
                sz = size(Data(d).Longitude);
                len_vecs = size(scattering_weights,1);  % JLL 26 May 2015 - find out how many pressure levels there are. Will often be 30, but might change.
                                                        % Need this to properly reshape the scattering weights, AKs, pressure levels, and (soon) profiles
                
                Data(d).BEHRAMFTrop = reshape(amf,sz); %JLL 18 Mar 2014: Save the resulting AMF of the pixel
                Data(d).BEHRGhostFraction = reshape(ghost_fraction,sz);
                Data(d).BEHRScatteringWeights = reshape(scattering_weights, [len_vecs, sz]);
                Data(d).BEHRAvgKernels = reshape(avg_kernels, [len_vecs, sz]);
                Data(d).BEHRNO2apriori = reshape(no2_prof_interp, [len_vecs, sz]);
                Data(d).BEHRPressureLevels = reshape(sw_plevels, [len_vecs, sz]);
            end
        end
        
        b=length(Data);
        for z=1:b;
            if isfield(Data,'BEHRAMFTrop')==0 || isempty(Data(z).BEHRAMFTrop)==1;
                continue
            else
                Data(z).BEHRColumnAmountNO2Trop=Data(z).ColumnAmountNO2Trop.*Data(z).AMFTrop./Data(z).BEHRAMFTrop;
                if DEBUG_LEVEL > 0; fprintf('   BEHR [NO2] stored for swath %u\n',z); end
            end
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
            'RelativeAzimuthAngle', [], 'AMFStrat', [], 'AMFTrop',[], 'CloudFraction', [], 'CloudRadianceFraction', [], 'CloudPressure', [], 'ColumnAmountNO2', [],...
            'SlantColumnAmountNO2', [], 'ColumnAmountNO2Trop', [], 'ColumnAmountNO2TropStd',[],'ColumnAmountNO2Strat',[],'TerrainHeight', [], 'TerrainPressure', [], 'TerrainReflectivity', [], 'vcdQualityFlags',{{}},...
            'MODISCloud', [], 'MODISAlbedo', [], 'GLOBETerpres', [], 'XTrackQualityFlags', {{}}, 'Row', [], 'Swath', [], 'TropopausePressure', [], 'BEHRColumnAmountNO2Trop',[],...
            'BEHRAMFTrop', [], 'Count', [], 'Area', [], 'Areaweight', [], 'MapData', struct);
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

        savename=[file_prefix,year,month,day];  
        if DEBUG_LEVEL > 0; disp(['   Saving data as',fullfile(behr_mat_dir,savename)]); end
        saveData(fullfile(behr_mat_dir,savename),Data,OMI)
        toc
        t=toc;
        %if t>1200
            %error('Time exceeded 20 min. Stopping')
        %end
    end
end
end

function saveData(filename,Data,OMI)
    save(filename,'OMI','Data')
end
