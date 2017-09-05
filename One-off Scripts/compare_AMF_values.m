% compare_AMF_values
%
%   Allows the user to get two AMF values using different profiles,
%   geometry, cloud fraction, etc.
%
%   Josh Laughner <joshlaugh5@gmail.com> 6 Aug 2014

lat = 39.4435;          lon = -76.3124;

pressures1 = wrfpres;   pressures2 = mpresbins;
profile1 = wrfprof;     profile2 = mno2bins;
sza1 = 29.2022;         sza2 = sza1;
vza1 = 57.0189;         vza2 = vza1;
phi1 = 155.13;        phi2 = phi1;
albedo1 = 0.0248;       albedo2 = albedo1;
surfPres1 = 999.17;     surfPres2 = surfPres1;
cloudPres1 = 1000;      cloudPres2 = cloudPres1;
cldFrac1 = 0;           cldFrac2 = cldFrac1;
cldRadFrac1 = 0;        cldRadFrac2 = cldFrac2;
month1 = 7;             month2 = 7;

noGhost = 1; ak = 0;

%%%%%%%%%%%%%%%%%%%%%%
%%%%% FILE PATHS %%%%%
%%%%%%%%%%%%%%%%%%%%%%

amf_tools_path = '/Users/Josh/Documents/MATLAB/BEHR/AMF_tools';
fileTmp = fullfile(amf_tools_path,'nmcTmpYr.txt');
fileDamf = fullfile(amf_tools_path,'damf.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The profiles and pressure must be columns for the AMF calculation to
% function correctly
if ~iscolumn(pressures1); pressures1 = pressures1'; end
if ~iscolumn(pressures2); pressures2 = pressures2'; end
if ~iscolumn(profile1); profile1 = profile1'; end
if ~iscolumn(profile2); profile2 = profile2'; end

%%%%%%%%%%%%%%%%%%%%%
%%%%%   AMF 1   %%%%%
%%%%%%%%%%%%%%%%%%%%%

[temperature, ~] = rNmcTmp2(fileTmp, pressures1, lon, lat, month1);
dAmfClr1 = rDamf2(fileDamf, pressures1, sza1, vza1, phi1, albedo1, surfPres1);
cloudalbedo=0.8;
dAmfCld1 = rDamf2(fileDamf, pressures1, sza1, vza1, phi1, cloudalbedo, cloudPres1);

amf1 = omiAmfAK2(surfPres1, cloudPres1, cldFrac1, cldRadFrac1, pressures1, dAmfClr1, dAmfCld1, temperature, profile1, profile1, noGhost, ak);

%%%%%%%%%%%%%%%%%%%%%
%%%%%   AMF 2   %%%%%
%%%%%%%%%%%%%%%%%%%%%

[temperature, ~] = rNmcTmp2(fileTmp, pressures2, lon, lat, month2);
dAmfClr2 = rDamf2(fileDamf, pressures2, sza2, vza2, phi2, albedo2, surfPres2);
cloudalbedo=0.8;
dAmfCld2 = rDamf2(fileDamf, pressures2, sza2, vza2, phi2, cloudalbedo, cloudPres2);

amf2 = omiAmfAK2(surfPres2, cloudPres2, cldFrac2, cldRadFrac2, pressures2, dAmfClr2, dAmfCld2, temperature, profile2, profile2, noGhost, ak);