function [ Data, OMI ] = BEHR_main_one_day( Data, varargin )
%BEHR_MAIN_ONE_DAY The BEHR algorithm for a single day's data.
%   [ DATA, OMI ] = BEHR_main_one_day( DATA ) Takes as input a DATA
%   structure created by READ_MAIN() and calculates BEHR AMFs for it as
%   well as grids the data using the BEHR-PSM repository. Returns the
%   structure DATA with the native pixels including BEHR NO2 VCDs and the
%   OMI structure with the gridded quantities.
%
%   Additional parameters:
%
%       'no2_profile_path' - the top directory to look for WRF output files
%       in. If not given, or given as an empty string, this is determined
%       automatically based on "profile_mode" and the region specified in
%       Data as the field "BEHRRegion".
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
%       'lookup_sweights' - scalar logical, determines whether the
%       algorithm should look up scattering weights based on SZA, VZA, etc.
%       or use the scattering weights already stored in Data. This is
%       intended to allow for uncertainty testing by ensuring the same
%       scattering weights are used while the NO2 profile is varied.
%       Default is true, i.e. the scattering weights are looked up from the
%       TOMRAD table and not read from Data.
%
%       'lookup_profile' - scalar logical, determined whether the WRF NO2
%       profile should be read from WRF netCDF files (true, default) or
%       read in from Data. Similar to "lookup_sweights", intended for
%       uncertainty analysis. If false, Data must be the result of
%       BEHR_main to have the profiles stored in it.
%
%       'DEBUG_LEVEL' - level of progress messaged printed to the console.
%       0 = none, 1 = minimal, 2 = all, 3 = processing times are added.
%       Default is 2.

p = inputParser;
% Parameters relevant to the normal retrieval
p.addParameter('no2_profile_path', '');
p.addParameter('profile_mode', 'monthly');
p.addParameter('use_psm_gridding', false);

% Parameters relevant to error analysis
p.addParameter('lookup_sweights', true);
p.addParameter('lookup_profile', true);

% Other parameters
p.addParameter('DEBUG_LEVEL', 2);

p.parse(varargin{:});
pout = p.Results;

no2_profile_path = pout.no2_profile_path;
prof_mode = pout.profile_mode;
use_psm = pout.use_psm_gridding;
lookup_sweights = pout.lookup_sweights;
lookup_profile = pout.lookup_profile;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

allowed_prof_modes = {'daily','monthly'};

if ~ischar(no2_profile_path)
    E.badinput('Parameter "no2_profile_path" must be a string');
elseif ~ismember(prof_mode,allowed_prof_modes)
    E.badinput('prof_mode (if given) must be one of %s', strjoin(allowed_prof_modes,', '));
elseif ~isscalar(use_psm) || (~islogical(use_psm) && ~isnumeric(use_psm))
    E.badinput('use_psm_gridding must be a scalar logical or number')
end

%Store paths to relevant files
fileDamf = fullfile(behr_paths.amf_tools_dir,'damf.txt');

% Get the Git head hashes
core_githead = git_head_hash(behr_paths.behr_core);
behrutils_githead = git_head_hash(behr_paths.behr_utils);
genutils_githead = git_head_hash(behr_paths.utils);
psm_githead = git_head_hash(behr_paths.psm_dir);
imatpy_githead = git_head_hash(behr_paths.python_interface);
wrfutils_githead = git_head_hash(behr_paths.wrf_utils);

this_date = Data(1).Date;
region = Data(1).BEHRRegion;

for d=1:length(Data);
    % Data is initialized in read_main with a single 0 in the Longitude
    % field.  Since points outside the lat/lons of interest are removed
    % completely, we should also check if all points are gone.
    if numel(Data(d).Longitude)==1 || isempty(Data(d).Longitude);
        if DEBUG_LEVEL > 1; fprintf('  Note: Data(%u) is empty\n',d); end
        continue %JLL 17 Mar 2014: Skip doing anything if there's really no information in this data
    end
    if DEBUG_LEVEL>0; fprintf('  Swath %u of %s \n',d,datestr(this_date)); end
    
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
    [no2Profile, temperature, wrf_profile_file, TropoPres,pindx, wrf_pres_mode, wrf_temp_mode] = rProfile_WRF(this_date, prof_mode, region, loncorns, latcorns, time, surfPres, pressure, no2_profile_path); %JLL 18 Mar 2014: Bins the NO2 profiles to the OMI pixels; the profiles are averaged over the pixel
    if ~lookup_profile
        no2Profile_check = no2Profile;
        no2Profile = remove_nonstandard_pressures(Data(d).BEHRNO2apriori, Data(d).BEHRPressureLevels, pressure);
        if ~lookup_sweights && any(abs(no2Profile(:) - no2Profile_check(:)) > 1e-12)
            % If not using the scattering weights from the Data structure,
            % then we need to verify that we loaded the right temperature
            % profiles for the scattering weights.
            E.callError('profile_lookup', 'Looked up different temperature profiles than the NO2 profiles given in the file (difference exceeds 1 pptv)')
        end
    end
    bad_profs = squeeze(all(isnan(no2Profile),1));
    
    if lookup_sweights
        if DEBUG_LEVEL > 1; fprintf('   Calculating clear and cloudy AMFs\n'); end
        dAmfClr = rDamf2(fileDamf, pressure, sza, vza, phi, albedo, surfPres); %JLL 18 Mar 2014: Interpolate the values in dAmf to the albedo and other conditions input
        cloudalbedo=0.8*ones(size(Data(d).CloudFraction)); %JLL 18 Mar 2014: Assume that any cloud has an albedo of 0.8
        dAmfCld = rDamf2(fileDamf, pressure, sza, vza, phi, cloudalbedo, cldPres); %JLL 18 Mar 2014: Interpolate dAmf again, this time taking the cloud top and albedo as the bottom pressure
    else
        dAmfClr = remove_nonstandard_pressures(Data(d).BEHRScatteringWeightsClear, Data(d).BEHRPressureLevels);
        dAmfCld = remove_nonstandard_pressures(Data(d).BEHRScatteringWeightsCloudy, Data(d).BEHRPressureLevels);
        % The scattering weights in the Data structures already include the
        % temperature correction, so we need to set the temperature
        % profiles to something that makes that correction 1, i.e. no
        % effect.
        temperature = 220 * ones(size(dAmfClr));
    end
    
    if DEBUG_LEVEL > 1; disp('   Calculating BEHR AMF'); end

    [amf, amfVis, ~, ~, scattering_weights_clear, scattering_weights_cloudy, avg_kernels, no2_prof_interp, sw_plevels] = omiAmfAK2(surfPres, TropoPres, cldPres, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile); %JLl 18 Mar 2014: The meat and potatoes of BEHR, where the TOMRAD AMF is adjusted to use the GLOBE pressure and MODIS cloud fraction
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
    Data(d).TropoPresVSCldPres = (TropoPres-cldPres) > 0;
    Data(d).Interp_TropopausePressure = pindx;
    Data(d).BEHRTropopausePressure = TropoPres;
    Data(d).BEHRQualityFlags = behr_quality_flags(Data(d));   
end
% remove the field 'TropoPresVSCldPres' as it's only used in behr_quality_flags
Data = rmfield(Data,{'TropoPresVSCldPres','Interp_TropopausePressure'});
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

end


function val_out = remove_nonstandard_pressures(val, pres, std_pres)
E = JLLErrors;
if size(val,1) ~= 30 || ndims(val) ~= 3
    E.badinput('VAL should be 3D with the first dimension having length 30')
end
if ~isequal(size(pres), size(val))
    E.badinput('PRES must be the same size as VAL')
end

sz = size(val);
val_out = nan([28, sz(2:end)]);
for a=1:prod(sz(2:end))
    if all(isnan(pres(:,a)))
        if ~all(isnan(val(:,a)))
            E.callError('find_std_pres', 'Trying to remove nonstandard pressures, but PRES(:,a) is all NaNs and VAL(:,a) is not');
        end
        xx = false(size(pres(:,a)));
        xx(1:28) = true;
    else
        xx = ismember(pres(:,a), std_pres);
    end
    if sum(xx) ~= 28
        E.callError('find_std_pres','Did not find the 28 standard pressures')
    end
    val_out(:,a) = val(xx,a);
end
end
