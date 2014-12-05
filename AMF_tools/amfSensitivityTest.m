function [ AMFs, SZAvec, VZAvec, RAAvec, ALBvec, SURFPRESSvec ] = amfSensitivityTest(profile, profile_pressures, lon, lat, month  )
%amfSensitivityTest Given a latitude, longitude, and profile, constructs a
%five-dimension matrix of AMFs
%   Scattering weights in the OMNO2 algorithm depend on SZA, VZA, RAA
%   (phi), surface albedo, and surface pressure.  They also depend on
%   profile shape, but that cannot be varied systematically as easily as
%   the other five parameters.  Given a profile, profile pressures,
%   latitude, and longitude (the last two are for the temperature profile,
%   and can be varied as well) this function will create a 5-7 dimensional
%   matrix of AMF values for different values in SZA, VZA, etc.  The last
%   two dimensions will be singletons if only a single lat/lon is
%   specified.

%%%%%%%%%%%%%%%%%%%%%%
%%%%% USER INPUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%

% Set to 0 if running this script without display or JVM (i.e. from command
% line), that will prevent the progress from from trying to be displayed.
waitbar_bool = 1;

% # SZAs between 0 and 88 degrees to test
nSZA = 9;

% # VZAs between 0 and 70 degrees to test
nVZA = 6;

% # RAA (rel. azimuth angles) between 0 and 180 degrees to test
nRAA = 11;

% Surface albedo can be between 0 and 1; but you can reset the min and max
% to restrict to i.e. 0 and 0.1 - this allows you to examine effects of
% albedo in a range appropriate to a ground type or clouds.
minAlb = 0;
maxAlb = 0.1;
nAlb = 6;

% Surface pressure can vary between 1013 and 0.003 hPa; but like albedo,
% you may reset the min and max.
minSurfPres = 880;
maxSurfPres = 1013;
nSurfPres = 5;

% set to 'OMI' (eventually perhaps 'GOME' will be an option) to determine
% what pressure levels to use
satellite = 'OMI';

% Location of the OMI temperature profiles and scattering weights
amf_tools_path = '/Users/Josh/Documents/MATLAB/BEHR/AMF_tools';
fileTmp = fullfile(amf_tools_path,'nmcTmpYr.txt');
fileDamf = fullfile(amf_tools_path,'damf.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

E=JLLErrors;

if ~isscalar(month)
    error(E.notScalar('month'));
end

if numel(profile) ~= numel(profile_pressures);
    error(E.numelMismatch('profile','profile_pressures'));
end
if numel(lon) ~= numel(lat)
    error(E.numelMismatch('lon','lat'));
end

if minAlb > maxAlb;
    warning('Min and max albedo backwards. Fixing.');
    tmp = maxAlb;
    maxAlb = minAlb;
    minAlb = tmp;
    clear('tmp');
end
if minAlb < 0;
    warning('Minimum albedo cannot be less than 0.  Clamping to 0.');
    minAlb = 0;
end
if maxAlb > 1;
	warning('Maximum albedo cannot be greater than 1. Clamping to 1.');
    maxAlb = 1;
end
if minSurfPres > maxSurfPres;
    warning('Min and max surface pressure backwards. Fixing.');
    tmp = maxSurfPres;
    maxSurfPres = minSurfPres;
    minSurfPres = tmp;
    clear('tmp');
end
if minSurfPres < 0.003;
    warning('Minimum surface pressure cannot be less than 0.003 hPa. Clamping to 0.003.');
    minSurfPres = 0.003;
end
if maxSurfPres > 1013;
	warning('Maximum surface pressure cannot be greater than 1013 hPa. Clamping to 1013.');
    maxSurfPres = 1013;
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN PROGRAM %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% Generate the vectors of SZA, VZA, RAA, albedo, and surface pressure that
% will generate 

SZAvec = linspace(0,88,nSZA);
VZAvec = linspace(0,70,nVZA);
RAAvec = linspace(0,180,nRAA);
ALBvec = linspace(minAlb, maxAlb, nAlb);
SURFPRESSvec = linspace(maxSurfPres, minSurfPres, nSurfPres);

fill_val = -127;

% Create the matrix 
AMFs = fill_val .* ones(numel(SZAvec), numel(VZAvec), numel(RAAvec), numel(ALBvec), numel(SURFPRESSvec), numel(lon), numel(lat));

% Prepare the pressure levels to interpolate to
switch satellite
    case 'OMI'
        presLevels = OMNO2PresLev; presLevels = presLevels';
        interp_profile = interp1(profile_pressures, profile, presLevels); % making presLevels a column vector ensures that the interpolated profile is a column
    otherwise
        error(E.callError('satellite','Satellite not recognized; could not determine pressure levels to use'))
end

if waitbar_bool
    wb = waitbar(0,'AMF Sensitivity Progress');
    wbcnt = 0; n=numel(AMFs);
end

% Many many loops to iterate over each variable
for x6 = 1:numel(lon);
    for x7 = 1:numel(lat);
        % Read the temperature profile for this lat/lon
        temperature = rNmcTmp2(fileTmp, presLevels, lon(x6), lat(x7), month);
        
        for x1=1:numel(SZAvec)
            if x1 == 1; fprintf('Now on SZA 1\n');
            else fprintf('\b\b2\n'); end
            for x2=1:numel(VZAvec)
                for x3=1:numel(RAAvec)
                    for x4=1:numel(ALBvec)
                        for x5=1:numel(SURFPRESSvec)
                            if waitbar_bool
                                wbcnt = wbcnt + 1;
                                waitbar(wbcnt/n,wb,'AMF Sensitivity Progress');
                            end
                            sza_i = SZAvec(x1);
                            vza_i = VZAvec(x2);
                            phi_i = RAAvec(x3);
                            albedo_i = ALBvec(x4);
                            surfPres_i = SURFPRESSvec(x5);
                            
                            dAmfClr = rDamf2(fileDamf, presLevels, sza_i, vza_i, phi_i, albedo_i, surfPres_i);
                            
                            cloudalbedo=0.8;
                            cloudPres_i = 500;
                            cldFrac_i = 0; 
                            cldRadFrac_i = 0;
                            dAmfCld = rDamf2(fileDamf, presLevels, sza_i, vza_i, phi_i, cloudalbedo, cloudPres_i);
                            
                            noGhost = 1; ak = 0;
                            amf_i = omiAmfAK2(surfPres_i, cloudPres_i, cldFrac_i, cldRadFrac_i, presLevels, dAmfClr, dAmfCld, temperature, interp_profile, interp_profile, noGhost, ak);
                            
                            AMFs(x1,x2,x3,x4,x5,x6,x7) = amf_i;
                        end
                    end
                end
            end
        end
    end
end

if waitbar_bool; close(wb); end
end
