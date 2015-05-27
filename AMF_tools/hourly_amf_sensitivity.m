function [ percent_diff, amf, amf0, Vectors ] = hourly_amf_sensitivity( BinStruct, A0_method, month_in, varargin )
%HOURLY_AMF_SENSITIVITY Test changes in AMF with the hourly change in profile
%   To try to really characterize how sensitive BEHR is to changes in the a
%   priori, this function takes hourly binned profiles (output from
%   bin_profile_by_start_times in the NO2 Profile repo/utility_scripts
%   folder) and uses them to test the change in AMF relative to some base
%   a priori, which is specified by the second argument. Options are:
%       'flat' - a profile with no vertical variation
%       'wrf' - use the nearest WRF-Chem profile used in BEHR
%       'geos' - use a GEOS-Chem profile (currently unimplemented)
%
%   month_in is a numeric month number.
%
%   The parameters n_sza, n_vza, and n_alb let you set how many SZAs, VZAs,
%   and albedos to test. If set to 1, SZA and VZA will be set to 30 and
%   albedo to 0.05. If >1, linspace is used to generate the inputs (between
%   0 and 75 for SZA, 0 and 60 for VZA, and 0.02 and 0.2 for albedo).
%   Default for all of these n's is 5.
%
%   Josh Laughner <joshlaugh5@gmail.com> 21 May 2015

E=JLLErrors;
DEBUG_LEVEL = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3 
    month_in = 0;
end

if ~isstruct(BinStruct)
    E.badinput('BinStruct must be a structure, output from bin_profile_by_start_time');
elseif ~isfield(BinStruct, 'final_no2_bins')
    E.badinput('Only pass one site (or the AllSites field) from the BinStruct structure')
end

A0_method = lower(A0_method);
if ~ischar(A0_method)
    E.badinput('A0_method must be a string');
elseif ~ismember(A0_method,{'flat','wrf','geos'})
    E.badinput('A0_method must be either ''flat'', ''wrf'', or ''geos''');
end

if ~isnumeric(month_in) || ~isscalar(month_in) || month_in < 1 || month_in > 12
    E.badinput('month_in must be a scalar, valid, numerical month');
end

% Parameters dealt with using an inputParser.
p = inputParser;
p.addParameter('n_sza', 5, @(x) (isnumeric(x) && isscalar(x) && x > 0));
p.addParameter('n_vza', 5, @(x) (isnumeric(x) && isscalar(x) && x > 0));
p.addParameter('n_alb', 5, @(x) (isnumeric(x) && isscalar(x) && x > 0));
p.parse(varargin{:});
pout = p.Results;

n_sza = pout.n_sza;
n_vza = pout.n_vza;
n_albs = pout.n_alb;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% conversion to get the aircraft data to parts-per-part, i.e. if [NO2] is
% given in ppt, this should be 1e-12.
insitu_unit_conv = 1e-12;

% file paths needed to get the partial AMFs
[dAmf_file, tmp_corr_file] = amf_filepaths;

no2_bins_mat = BinStruct.final_no2_bins;
pres_bins_mat = BinStruct.final_pres_bins;
start_hours = BinStruct.bin_start_hours;
lon = BinStruct.Lon;
lat = BinStruct.Lat;

if isnan(lon)
    E.badvar('A longitude must be defined in BinStruct to get the proper temperature correction')
elseif isnan(lat)
    E.badvar('A latitude must be defined in BinStruct to get the proper temperature correction')
end


% We're going to just use clear sky assumptions; since we nearly always
% filter for cloud fraction < 20-30%, this is reasonable. Use 30 degrees
% for SZA and VZA and 0.05 for albedo because those are relatively common
% values in one OMI swath over the US.

if n_sza > 1
    szas = linspace(0,75,n_sza);
else
    szas = 30;
end

if n_vza > 1
    vzas = linspace(0,60,n_vza);
else
    vzas = 30;
end

if n_albs > 1
    albs = linspace(0.02,0.2,n_albs);
else
    albs = 0.05;
end

raa = 100; % 100 degrees is a fairly common rel. azimuth angle in an OMI swath

n_profs = numel(start_hours);
n_bins = size(no2_bins_mat,1);

amf = nan(n_profs, n_sza, n_vza, n_albs);
amf0 = nan(n_profs, n_sza, n_vza, n_albs);

% Load the profile that will be used as the A0 case
switch A0_method
    case 'flat'
        prof_A0 = 100e-12*ones(n_bins,1); % set to ~100 ppt (shouldn't matter - normalize by VCD, but why not)
    case 'wrf'
        no2_profile_path = '/Volumes/share-sat/SAT/BEHR/Monthly_NO2_Profiles';
        profile_file = sprintf('m%02d_NO2_profile',month_in);
        S=load(fullfile(no2_profile_path,profile_file));
        PROFILE = S.PROFILE;
        [model_prof, model_pres] = interp_wrf_prof(PROFILE,lat,lon);
        model_prof = model_prof .* 1e-6; % WRF profiles are in ppm in these files
    case 'geos'
        E.notimplemented('geos')
end

l=0; % variable for printing debugging messages
for s=1:n_sza
    for v=1:n_vza
        for a=1:n_albs
            for p=1:n_profs
                if DEBUG_LEVEL > 0
                    str = sprintf('SZA %d of %d; VZA %d of %d; ALB %d of %d; profile %d of %d',s,n_sza,v,n_vza,a,n_albs,p,n_profs);
                    if l > 0;
                        erasestr = repmat('\b',1,l);
                        fprintf(erasestr);
                    end
                    fprintf(str);
                    l = length(str);
                end
                not_nans = ~isnan(no2_bins_mat(:,p));
                if sum(not_nans) > 0
                    surfP = max(pres_bins_mat(not_nans,p));
                    surfP = clip(surfP,200,1013); % the scattering weight lookup cannot handle surface pressures outside this range
                    sweights = rDamf2(dAmf_file, pres_bins_mat(:,p), szas(s), vzas(v), raa, albs(a), surfP);
                    temperature = rNmcTmp2(tmp_corr_file, pres_bins_mat(:,p), lon, lat, month_in);
                    if ismember(A0_method,{'wrf','geos'})
                        prof_A0 = interp_model_prof(model_prof, model_pres, pres_bins_mat(:,p));
                    end
                    amf(p,s,v,a) = omiAmfAK2(surfP, surfP, 0, 0, pres_bins_mat(:,p), sweights, sweights, temperature, no2_bins_mat(:,p).*insitu_unit_conv, no2_bins_mat(:,p).*insitu_unit_conv, 0, 0);
                    amf0(p,s,v,a) = omiAmfAK2(surfP, surfP, 0, 0, pres_bins_mat(:,p), sweights, sweights, temperature, prof_A0, prof_A0, 0, 0);
                else
                    amf(p,s,v,a) = nan;
                    amf0(p,s,v,a) = nan;
                end
            end
        end
    end
end
if DEBUG_LEVEL > 0;
    fprintf('\n');
end

percent_diff = (amf ./ amf0 - 1)*100;
Vectors.SZA = szas;
Vectors.VZA = vzas;
Vectors.ALB = albs;
Vectors.Hours = start_hours;

end

function prof_out = interp_model_prof(model_prof, model_pres, pres_in)
if ~iscolumn(model_prof); model_prof = model_prof'; end
if ~iscolumn(model_pres); model_pres = model_pres'; end
if ~iscolumn(pres_in); pres_in = pres_in'; end
% Do the interpolation in log-log space since we generally expect NO2 to
% vary exponentially with pressure
model_prof = log(model_prof);
model_pres = log(model_pres);
pres_in = log(pres_in);
prof_out = interp1(model_pres, model_prof, pres_in);
prof_out = exp(prof_out);
end

