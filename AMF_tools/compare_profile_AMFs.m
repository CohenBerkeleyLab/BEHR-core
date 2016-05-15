function [ amf1, amf2 ] = compare_profile_AMFs( profile1, pressures1, profile2, pressures2 )
%compare_profile_AMFs Calculate AMFs for two different profiles


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

