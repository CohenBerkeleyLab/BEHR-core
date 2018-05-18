function [ Data, OMI ] = BEHR_main_one_day( Data, varargin )
%BEHR_MAIN_ONE_DAY The BEHR algorithm for a single day's data.
%   [ DATA, OMI ] = BEHR_main_one_day( DATA ) Takes as input a DATA
%   structure created by READ_MAIN() and calculates BEHR AMFs for it as
%   well as grids the data using the BEHR-PSM repository. Returns the
%   structure DATA with the native pixels including BEHR NO2 VCDs and the
%   OMI structure with the gridded quantities.
%
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
%       'err_wrf_missing_attr' - if true (default), then if WRF files are
%       missing attributes that are read in (usually units), an error is
%       thrown. However, if false, then default units are assumed. Use
%       "false" with caution, as if the necessary variables are given in
%       the wrong units, there will be no way to catch that if "false" is
%       given for this parameter.
%
%       'extra_gridding_fields' - a cell array of strings that lists extra
%       fields that you wish to have gridded, beyond the standard fields
%       listed in BEHR_publishing_gridded_fields.
%
%       'DEBUG_LEVEL' - level of progress messaged printed to the console.
%       0 = none, 1 = minimal, 2 = all, 3 = processing times are added.
%       Default is 2.
%
%
%   Parameters specific to error analysis:
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
%       'randomize_profile_time' - scalar logical (default false), if true,
%       the NO2 profile will be chosen from a different day in the same
%       month. The hour of day is still chosen properly, so the dominant
%       error that this is testing is if the wind direction is very wrong.
%
%       'randomize_profile_loc' - scalar logical (default false), if true,
%       then the pixel corners will be "jittered" by some amount (currently
%       0.2 degrees) so that a WRF profile from a nearby, but different,
%       location is used for that pixel. All corners of one pixel are
%       jittered in the same direction. This simulates if the emissions,
%       transport speed, or chemistry are in the incorrect place in WRF.
%       Currently, the randomization allows the pixels to not move, in
%       order to include the possibility that the
%       emissions/transport/chemistry are correct in the WRF simulation.


p = inputParser;
% Parameters relevant to the normal retrieval
p.addParameter('no2_profile_path', '');
p.addParameter('profile_mode', 'monthly');
p.addParameter('use_psm_gridding', false);
p.addParameter('err_wrf_missing_attr', true);
p.addParameter('extra_gridding_fields', {});

% Parameters relevant to error analysis
p.addParameter('lookup_sweights', true);
p.addParameter('lookup_profile', true);
p.addParameter('randomize_profile_time', false);
p.addParameter('randomize_profile_loc', false);

% Other parameters
p.addParameter('DEBUG_LEVEL', 2);

p.KeepUnmatched = true;
p.parse(varargin{:});
pout = p.Results;

no2_profile_path = pout.no2_profile_path;
prof_mode = pout.profile_mode;
use_psm = pout.use_psm_gridding;
err_wrf_missing_attr = pout.err_wrf_missing_attr;
extra_gridding_fields = pout.extra_gridding_fields;
lookup_sweights = pout.lookup_sweights;
lookup_profile = pout.lookup_profile;
randomize_profile_time = pout.randomize_profile_time;
randomize_profile_loc = pout.randomize_profile_loc;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

allowed_prof_modes = {'daily','monthly'};

if ~ischar(no2_profile_path)
    E.badinput('Parameter "no2_profile_path" must be a string');
elseif ~ismember(prof_mode,allowed_prof_modes)
    E.badinput('prof_mode (if given) must be one of %s', strjoin(allowed_prof_modes,', '));
elseif ~isscalar(use_psm) || (~islogical(use_psm) && ~isnumeric(use_psm))
    E.badinput('use_psm_gridding must be a scalar logical or number')
end

if ~iscellstr(extra_gridding_fields)
    E.badinput('extra_gridding_fields must be a cell array of char arrays')
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

for d=1:length(Data)
    % Data is initialized in read_main with a single 0 in the Longitude
    % field.  Since points outside the lat/lons of interest are removed
    % completely, we should also check if all points are gone.
    if numel(Data(d).Longitude)==1 || isempty(Data(d).Longitude)
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
    globe_terheight = Data(d).GLOBETerrainHeight;
    albedo = Data(d).MODISAlbedo;
    cldFrac = Data(d).CloudFraction;
    cldRadFrac = Data(d).CloudRadianceFraction;
    
    pressure = behr_pres_levels();
    
    if DEBUG_LEVEL > 1; fprintf('   Reading NO2 and temperature profiles\n'); end
    if randomize_profile_time
        % Used for uncertainty analysis to vary the day of month that the
        % profile is chosen from. 
        if strcmpi(prof_mode, 'daily')
            prof_date = random_day_in_month(this_date, true);
        else
            E.callError('incompatible_prof_mode', '"randomize_profile_time" cannot be used unless "prof_mode" is "daily"');
        end
    else
        prof_date = this_date;
    end
    
    if randomize_profile_loc
        % Used for uncertainty analysis. Move the pixel corners by roughly
        % one OMI pixel in each dimension to "jitter" which profiles are
        % used. This should emulate if e.g. emissions are allocated in the
        % wrong place, or transport is wrong, or chemistry is wrong,
        % because by moving the pixel we can alter how long the air parcel
        % has aged chemically by moving it nearer or farther from the
        % source.
        [prof_loncorns, prof_latcorns] = jitter_corners(loncorns, latcorns, 0.2);
    else
        prof_loncorns = loncorns;
        prof_latcorns = latcorns;
    end
        
    % For normal runs, we want the NO2 profiles and temperature profiles
    % set to NaN outside the integration limits (actually before the bin
    % just below the surface and after the bin just above the tropopause).
    % However, when doing error analysis, we sometimes run into issues
    % where the tropopause is different whether using precalculated or
    % online computed temperature profiles. (The difference in the
    % temperature profile itself is very small, but it is occasionally
    % enough to move the lapse rate to the other side of the 2 K/km
    % threshold.) Therefore in that case we need to keep the NO2 and
    % temperature profiles over all bins.
    keep_all_bins = ~lookup_profile;
    [no2Profile, temperature, wrf_profile_file, surfPres, surfPres_WRF, tropoPres, tropopause_interp_flag, wrf_pres_mode, wrf_temp_mode] = ...
        rProfile_WRF(prof_date, prof_mode, region, prof_loncorns, prof_latcorns, time, globe_terheight, pressure, no2_profile_path,...
        'err_missing_att', err_wrf_missing_attr, 'clip_at_int_limits', ~keep_all_bins); %JLL 18 Mar 2014: Bins the NO2 profiles to the OMI pixels; the profiles are averaged over the pixel
    
    surfPres(surfPres > 1013) = 1013;
    cldPres = Data(d).CloudPressure;
    cldPres = min(cldPres, surfPres); % Clamp cldPres to be <= surface pressure, should not have below surface clouds.
    
    if ~lookup_profile
        % If we want to use the exact same NO2 profiles as in the original
        % run, we can't use the ones in the Data file directly because we 
        % might need the profile extrapolated over all the standard
        % pressure levels. Since rProfile_WRF has an option to leave those
        % in, instead of trying to extrapolate the stored profiles, we just
        % check that there are no differences greater than 1 pptv. This will
        % not fail if there are NaNs in the stored profiles but not the
        % extrapolated ones because NaN > x will always be false.
        no2Profile_check = remove_nonstandard_pressures(Data(d).BEHRNO2apriori, Data(d).BEHRPressureLevels, pressure);
        if ~lookup_sweights && any(abs(no2Profile(:) - no2Profile_check(:)) > 1e-12)
            % If not using the scattering weights from the Data structure,
            % then we need to verify that we loaded the right temperature
            % profiles for the scattering weights.
            E.callError('profile_lookup', 'Looked up different temperature profiles than the NO2 profiles given in the file (difference exceeds 1 pptv)')
        end
        tropoPres = Data(d).BEHRTropopausePressure;
    end
    bad_profs = squeeze(all(isnan(no2Profile),1));
    
    if lookup_sweights || randomize_profile_loc || randomize_profile_time
        % If we change the profiles, we need to allow for the tropopause to
        % have changed and the effect of the different temperature profile.
        % The latter should be minimal, but I've run into cases where the
        % old temperature wasn't defined at enough pressure levels for the
        % new calculation to work (i.e. there were NaNs left during the
        % integration).
        if DEBUG_LEVEL > 1; fprintf('   Calculating clear and cloudy AMFs\n'); end
        dAmfClr = rDamf2(fileDamf, pressure, sza, vza, phi, albedo, surfPres); %JLL 18 Mar 2014: Interpolate the values in dAmf to the albedo and other conditions input
        cloudalbedo=0.8*ones(size(Data(d).CloudFraction)); %JLL 18 Mar 2014: Assume that any cloud has an albedo of 0.8
        dAmfCld = rDamf2(fileDamf, pressure, sza, vza, phi, cloudalbedo, cldPres); %JLL 18 Mar 2014: Interpolate dAmf again, this time taking the cloud top and albedo as the bottom pressure
    else
        % Really we will almost never get here; if we're running BEHR
        % normally, then we obviously need to look up the scattering
        % weights; if we're running an error analysis and perturb one of
        % the scattering weight input parameters, we need to look up the
        % scattering weights; if we're running an error analysis and
        % perturb the profiles, we need to look up the scattering weights
        % again because the tropopause may have changed. I'm leaving this
        % here in case it becomes useful in the future.
        dAmfClr = remove_nonstandard_pressures(Data(d).BEHRScatteringWeightsClear, Data(d).BEHRPressureLevels, pressure);
        dAmfCld = remove_nonstandard_pressures(Data(d).BEHRScatteringWeightsCloudy, Data(d).BEHRPressureLevels, pressure);
        % The scattering weights in the Data structures already include the
        % temperature correction, so we need to set the temperature
        % profiles to something that makes that correction 1, i.e. no
        % effect.
        temperature = 220 * ones(size(dAmfClr));
    end
    
    if DEBUG_LEVEL > 1; disp('   Calculating BEHR AMF'); end

    [amf, amfVis, ~, ~, scattering_weights_clear, scattering_weights_cloudy, avg_kernels, no2_prof_interp, sw_plevels] = omiAmfAK2(surfPres, tropoPres, cldPres, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile); %JLl 18 Mar 2014: The meat and potatoes of BEHR, where the TOMRAD AMF is adjusted to use the GLOBE pressure and MODIS cloud fraction
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
    % temporary fields, will be removed after the warning flag is set
    Data(d).TropoPresVSCldPres = (tropoPres-cldPres) > 0;
    Data(d).Interp_TropopausePressure = tropopause_interp_flag;
    %
    Data(d).BEHRSurfacePressure = surfPres;
    Data(d).WRFSurfacePressure = surfPres_WRF; % mainly for testing, I'm curious how much WRF's surface pressure differs when adjusted with GLOBE
    Data(d).BEHRTropopausePressure = tropoPres;
    Data(d).BEHRQualityFlags = behr_quality_flags(Data(d));   
end
% remove the field 'TropoPresVSCldPres' as it's only used in behr_quality_flags
Data = rmfield(Data,{'TropoPresVSCldPres','Interp_TropopausePressure'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATE VCDS FROM NASA SCDS AND OUR AMFS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b=length(Data);
for z=1:b
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

OMI = psm_wrapper(Data, Data(1).Grid, 'only_cvm', ~use_psm, 'extra_cvm_fields', extra_gridding_fields, 'DEBUG_LEVEL', DEBUG_LEVEL);

end


function val_out = remove_nonstandard_pressures(val, pres, std_pres)
n_std_pres = numel(std_pres);

E = JLLErrors;
if size(val,1) ~= n_std_pres + 3 || ndims(val) ~= 3
    E.badinput('VAL should be 3D with the first dimension having length %d', n_std_pres + 3)
end
if ~isequal(size(pres), size(val))
    E.badinput('PRES must be the same size as VAL')
end

sz = size(val);
val_out = nan([n_std_pres, sz(2:end)]);
for a=1:prod(sz(2:end))
    if all(isnan(pres(:,a)))
        if ~all(isnan(val(:,a)))
            E.callError('find_std_pres', 'Trying to remove nonstandard pressures, but PRES(:,a) is all NaNs and VAL(:,a) is not');
        end
        xx = false(size(pres(:,a)));
        % If all the pressures are NaNs for this profile, then there's no
        % way that we can find the standard pressures. This usually happens
        % when the pixel is outside of the WRF domain, so there's no
        % profile data of any sort. As long as that is true (i.e. VAL(:,a)
        % is all NaNs, checked above), just take the first n NaNs so that
        % the value out is the right size.
        xx(1:n_std_pres) = true;
    else
        xx = ismember(pres(:,a), std_pres);
    end
    if sum(xx) ~= n_std_pres
        E.callError('find_std_pres','Did not find the %d standard pressures', n_std_pres)
    end
    val_out(:,a) = val(xx,a);
end
end

function date_out = random_day_in_month(date_in, forbid_current_day)
if ~exist('forbid_current_day', 'var')
    forbid_current_day = false;
end
y = year(date_in);
m = month(date_in);
d = day(date_in);
while d == day(date_in)
    d = randi(eomday(y,m));
    if ~forbid_current_day
        % If we allow the randomization to pick the same day as the input
        % date, then we can always exit the first time.
        break
    end
end
date_out = datenum(y,m,d);
end

function [loncorn_out, latcorn_out] = jitter_corners(loncorn_in, latcorn_in, jitter_amt, force_jitter)
if ~exist('force_jitter', 'var')
    force_jitter = false;
end

E = JLLErrors;
if size(loncorn_in, 1) ~= 4 || size(latcorn_in, 1) ~= 4
    E.badinput('LONCORN_IN and LATCORN_IN must have length 4 in the first dimension')
elseif ~isequal(size(loncorn_in), size(latcorn_in))
    E.badinput('LONCORN_IN and LATCORN_IN must be the same size')
end
if ~isnumeric(jitter_amt) || ~isscalar(jitter_amt) || jitter_amt <= 0
    E.badinput('JITTER_AMT must be a scalar, positive number')
end

sz = size(loncorn_in);
x_sign = randi([-1 1], sz(2:3));
y_sign = randi([-1 1], sz(2:3));

% Okay, here's the slightly tricky part. If we want to force a pixel to
% jitter, then we need to rerandomize any locations where the x and y shift
% are both zero.
if force_jitter
    both_zero = x_sign(:) == 0 & y_sign(:) == 0;
    while any(both_zero(:))
        x_sign(both_zero) = randi([-1 1], sum(both_zero), 1);
        y_sign(both_zero) = randi([-1 1], sum(both_zero), 1);
        both_zero = x_sign(:) == 0 & y_sign(:) == 0;
    end
end

% Now we just need to create the actual jitter matrices so that all the
% corners for a given pixel move by the same amount.
x_shift = repmat(permute(x_sign, [3 1 2]) * jitter_amt, 4, 1, 1);
y_shift = repmat(permute(y_sign, [3 1 2]) * jitter_amt, 4, 1, 1);

loncorn_out = loncorn_in + x_shift;
latcorn_out = latcorn_in + y_shift;

end
