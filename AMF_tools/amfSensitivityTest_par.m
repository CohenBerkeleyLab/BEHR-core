function [ Out ] = amfSensitivityTest_par(profile, profile_pressures, lon, lat, month, skytype  )
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

DEBUG_LEVEL = 2;

%%%%%%%%%%%%%%%%%%%
%%%%% STARTUP %%%%%
%%%%%%%%%%%%%%%%%%%

global onCluster;
if isempty(onCluster);
    onCluster = false;
end

% Paths to the dependencies within this file - necessary to run on the
% cluster. Set the global variable 'onCluster' to true as part of a Matlab
% run script that calls this function.

if onCluster
    clusterMatlabDir = '~/MATLAB';
    addpath(genpath(fullfile(clusterMatlabDir,'BEHR')));
    addpath(genpath(fullfile(clusterMatlabDir,'Utils')));
    addpath(genpath(fullfile(clusterMatlabDir,'Classes')));
end

% Instantiate the custom error class
E=JLLErrors;

% The number of threads to allow to run. Likewise, set this in a MATLAB
% run script calling this function, otherwise it will default to 1.
global numThreads;
if isempty(numThreads); numThreads = 1; end
if ~isscalar(numThreads); E.badvartype(numThreads, 'scalar'); end


%%%%%%%%%%%%%%%%%%%%%%
%%%%% USER INPUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%

% # SZAs between 0 and 88 degrees to test
nSZA = 9;

% # VZAs between 0 and 70 degrees to test
nVZA = 6;

% # RAA (rel. azimuth angles) between 0 and 180 degrees to test
nRAA = 5;

% Surface albedo can be between 0 and 1; but you can reset the min and max
% to restrict to i.e. 0 and 0.1 - this allows you to examine effects of
% albedo in a range appropriate to a ground type or clouds.
if strcmpi(skytype,'clear')
    minAlb = 0;
    maxAlb = 0.08;
else
    minAlb = 0.7;
    maxAlb = 0.9;
end
nAlb = 10;

% Surface pressure can vary between 1013 and 0.003 hPa; but like albedo,
% you may reset the min and max.
if strcmpi(skytype,'clear')
    minSurfPres = 795;
    maxSurfPres = 1013;
else
    minSurfPres = 344;
    maxSurfPres = 1003;
end
nSurfPres = 10;

% set to 'OMI' (eventually perhaps 'GOME' will be an option) to determine
% what pressure levels to use
satellite = 'OMI';

% Location of the OMI temperature profiles and scattering weights
amf_tools_path = fileparts(mfilename('fullpath'));
fileTmp = fullfile(amf_tools_path,'nmcTmpYr.txt');
fileDamf = fullfile(amf_tools_path,'damf.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

SZAvec = reshape(linspace(0,88,nSZA),[],1,1,1,1);
VZAvec = reshape(linspace(0,70,nVZA),1,[],1,1,1);
RAAvec = reshape(linspace(0,180,nRAA),1,1,[],1,1);
ALBvec = reshape(linspace(minAlb, maxAlb, nAlb),1,1,1,[],1);
SURFPRESSvec = reshape(linspace(maxSurfPres, minSurfPres, nSurfPres),1,1,1,1,[]);

sz = [numel(SZAvec), numel(VZAvec), numel(RAAvec), numel(ALBvec), numel(SURFPRESSvec)];
SZAmat = repmat(SZAvec, 1, sz(2), sz(3), sz(4), sz(5));
VZAmat = repmat(VZAvec, sz(1), 1, sz(3), sz(4), sz(5));
RAAmat = repmat(RAAvec, sz(1), sz(2), 1, sz(4), sz(5));
ALBmat = repmat(ALBvec, sz(1), sz(2), sz(3), 1, sz(5));
SURFPRESSmat = repmat(SURFPRESSvec, sz(1), sz(2), sz(3), sz(4), 1);

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

% If no parallel pool is active, start one. If there wasn't one, go ahead
% and close it at the end.
if isempty(gcp('nocreate')) && onCluster
    parpool(numThreads);
    n_workers = numThreads;
    closeparpool = true;
else
    closeparpool = false;
    n_workers = 0;
end

nlon = numel(lon);
nlat = numel(lat);



% By constructing the inputs as a matrix we can fully parallelize it rather
% than only be able to parallelize say the SZA loop.
for xlon = 1:nlon;
    for xlat = 1:nlat;
        % Read the temperature profile for this lat/lon
        temperature = rNmcTmp2(fileTmp, presLevels, lon(xlon), lat(xlat), month);
        
        parfor(x=1:numel(AMFs), n_workers)
            sza_i = SZAmat(x);
            vza_i = VZAmat(x);
            phi_i = RAAmat(x);
            albedo_i = ALBmat(x);
            surfPres_i = SURFPRESSmat(x);
            
            if DEBUG_LEVEL > 1
                t = getCurrentTask;
                if isempty(t)
                    t.ID = 0;
                end
                fprintf('Worker %d: SZA = %f, VZA = %f, RAA = %f, ALB = %f, SURFP = %f\n', t.ID, sza_i, vza_i, phi_i, albedo_i, surfPres_i);
            end
            
            dAmfClr = rDamf2(fileDamf, presLevels, sza_i, vza_i, phi_i, albedo_i, surfPres_i);
            
            cloudalbedo=0.8;
            cloudPres_i = 500;
            cldFrac_i = 0;
            cldRadFrac_i = 0;
            dAmfCld = rDamf2(fileDamf, presLevels, sza_i, vza_i, phi_i, cloudalbedo, cloudPres_i);
            
            amf_i = omiAmfAK2(surfPres_i, cloudPres_i, cldFrac_i, cldRadFrac_i, presLevels, dAmfClr, dAmfCld, temperature, interp_profile);
            
            AMFs(x) = amf_i;
        end
    end
end

Out.AMFs = AMFs;
Out.SZAs = SZAmat;
Out.VZAs = VZAmat;
Out.RAAs = RAAmat;
Out.ALBs = ALBmat;
Out.SurfPs = SURFPRESSmat;

if closeparpool
    p = gcp('nocreate');
    delete(p);
end
end
