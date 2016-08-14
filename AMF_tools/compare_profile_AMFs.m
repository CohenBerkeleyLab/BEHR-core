function [ amf1, amf2, ancilliary ] = compare_profile_AMFs( profile1, pressures1, profile2, pressures2, ancilliary )
%COMPARE_PROFILE_AMFS Calculates AMFs for two profiles
%   [ AMF1, AMF2 ] = COMPARE_PROFILES_AMFS( PROFILE1, PRESSURES1, PROFILE2, PRESSURES2 )
%       Computes the AMFs for PROFILE1 and PROFILE2 that are defined at
%       PRESSURES1 and PRESSURES2 as the vertical coordinates. The profiles
%       must be given as mixing ratios (not number density). They are
%       typically given as unscaled mixing ratio (i.e. 1 ppb = 1e-9) but
%       this shouldn't matter as long as they are all in the same unit.
%
%   [ ___ ] = COMPARE_PROFILES_AMFS( PROFILE1, PRESSURES1, PROFILE2, PRESSURES2, ANCILLIARY )
%       Allows you to define some or all of the ancilliary parameters using
%       the structure ANCILLIARY with fields:
%           sza = solar zenith angle (0 to 90)
%           vza = viewing zenith angle (0 to 90)
%           raa = relative azimuth angle (0 to 180)
%           alb = albedo (0 to 1)
%           surfp = surface pressure (<= 1013)
%           cldfrac = geometric cloud fraction (0 to 1)
%           cldradfrac = radiance cloud fraction (0 to 1, usually about 50%
%               greater than the geometric fraction)
%           lon = longitude (-180 to 180)
%           lat = latitude (-90 to 90)
%           month = the numeric month (1 to 12), used for the temperature
%               correction
%       Any value not given as part of the structure will be assigned the
%       default value defined in this function.
%
%   [ AMF1, AMF2, ANCILLIARY ] = COMPARE_PROFILE_AMFS( ___ ) returns the
%   additional structure ANCILLIARY which contains the values of the
%   ancilliary parameters used. This lets you verify that your parameters
%   were used or see what the default values are.

    % Set the scattering weight parameters
    SZA = 55;
    VZA = 55;
    RAA = 100;
    ALB = 0.02;
    SurfP = 1013;
    cldFrac = 0.3;
    cldRadFrac = 0.5;
    
    % Set geographic parameters (used for temperature correction)
    lon = 0;
    lat = 0;
    this_month = 6;
    
if exist('ancilliary','var')
    if isfield(ancilliary,'sza')
        SZA = ancilliary.sza;
    end
    if isfield(ancilliary,'vza')
        VZA = ancilliary.vza;
    end
    if isfield(ancilliary,'raa')
        RAA = ancilliary.raa;
    end
    if isfield(ancilliary,'alb')
        ALB = ancilliary.alb;
    end
    if isfield(ancilliary,'surfp')
        SurfP = ancilliary.surfp;
    end
    if isfield(ancilliary,'cldfrac')
        cldFrac = ancilliary.cldfrac;
    end
    if isfield(ancilliary, 'cldradfrac');
        cldRadFrac = ancilliary.cldradfrac;
    end
    if isfield(ancilliary, 'lon')
        lon = ancilliary.lon;
    end
    if isfield(ancilliary, 'lat')
        lat = ancilliary.lat;
    end
    if isfield(ancilliary,'month')
        this_month = ancilliary.month;
    end
end

ancilliary.sza = SZA;
ancilliary.vza = VZA;
ancilliary.raa = RAA;
ancilliary.alb = ALB;
ancilliary.surfp = SurfP;
ancilliary.cldfrac = cldFrac;
ancilliary.cldradfrac = cldRadFrac;
ancilliary.lon = lon;
ancilliary.lat = lat;
ancilliary.month = this_month;

% File locations
[fileDamf,fileTmp] = amf_filepaths;


% Check the inputs; everything must be a column and the number of elements
% in the profile and its corresponding pressure must be the same
E = JLLErrors;
if ~isvector(profile1) || ~isvector(profile2) || ~isvector(pressures1) || ~isvector(pressures2)
    E.badinput('All inputs to this function must be column vectors (row vectors will be converted into column vectors)');
end
if isrow(profile1); profile1 = profile1'; end
if isrow(profile2); profile2 = profile2'; end
if isrow(pressures1); pressures1 = pressures1'; end
if isrow(pressures2); pressures2 = pressures2'; end
if ~all(size(profile1)==size(pressures1))
    E.dimMismatch('profile1','pressures1');
elseif ~all(size(profile2)==size(pressures2))
    E.dimMismatch('profile2','pressures2');
end
    

dAmfClr1 = rDamf2(fileDamf, pressures1, SZA, VZA, RAA, ALB, SurfP);
dAmfClr2 = rDamf2(fileDamf, pressures2, SZA, VZA, RAA, ALB, SurfP);

cloudalbedo=0.8;
cloudPres = 500;
dAmfCld1 = rDamf2(fileDamf, pressures1, SZA, VZA, RAA, cloudalbedo, cloudPres);
dAmfCld2 = rDamf2(fileDamf, pressures2, SZA, VZA, RAA, cloudalbedo, cloudPres);

temperature1 = rNmcTmp2(fileTmp, pressures1, lon, lat, this_month);
temperature2 = rNmcTmp2(fileTmp, pressures2, lon, lat, this_month);

noGhost = 1; ak = 0;
amf1 = omiAmfAK2(SurfP, cloudPres, cldFrac, cldRadFrac, pressures1, dAmfClr1, dAmfCld1, temperature1, profile1, profile1, noGhost, ak);
amf2 = omiAmfAK2(SurfP, cloudPres, cldFrac, cldRadFrac, pressures2, dAmfClr2, dAmfCld2, temperature2, profile2, profile2, noGhost, ak);


end

