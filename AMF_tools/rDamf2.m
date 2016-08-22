%%rDamf
%%arr 07/23/2008

% Reads the dAMF file used in the FORTRAN OMI NO2 algorithm

%..........................................................................
% Returns vector of dAMFs on input pressure grid, interpolated
%  at input solar zenith angle, viewing zenith angle, relative azimuth angle
%  surface albedo (cloud or terrain) and surface pressure (cloud or terrain)
%
% Inputs:
%  fileDamf = complete directory/filename of input file
%  pressure = vector of pressures (hPa) monotonically decreasing
%  sza, vza, phi, (all in deg), albedo, surfPres (hPa) are scalars
%
% Outputs (optional: Return these to avoid having to re-read fileDamf on each call):
%  presSave = pressure grid (vector) from fileDamf
%  szaSave  = solar zenith angle vector from fileDamf
%  vzaSave  = viewing zenith angle vector from fileDamf 
%  phiSave  = relative azimuth angle vector from fileDamf
%  albedoSave = albedo vector from fileDamf
%  surfPresSave = surface pressure vector from fileDamf
%  dAmfSave = large array of dAMF values from fileDamf
%..........................................................................

%function [presSave, szaSave, vzaSave, phiSave, albedoSave, surfPresSave, dAmfSave, dAmf] = rDamf2(fileDamf, pressure, sza, vza, phi, albedo, surfPres)
function dAmf = rDamf2(fileDamf, pressure, sza, vza, phi, albedo, surfPres)

E=JLLErrors;

if any(diff(pressure) > 0);
    E.badinput('pressure must be a monotonically decreasing vector of pressures');
else
    if exist('dAmfSave','var')==0;

        fid = fopen(fileDamf,'r');
        
        header = fgets(fid);

        %nPresSave = 0;  nSzaSave = 0;  nVzaSave = 0;  nPhiSave = 0;  nAlbedoSave = 0;  nSurfPresSave = 0;
        
        n = fscanf(fid,'%g',[1 6]);
        
        nPresSave = n(1); nSzaSave = n(2); nVzaSave = n(3); nPhiSave = n(4); nAlbedoSave = n(5); nSurfPresSave = n(6);
        
        %presSave     = zeros(nPresSave,1);
        %szaSave      = zeros(nSzaSave,1);
        %vzaSave      = zeros(nVzaSave,1);
        %phiSave      = zeros(nPhiSave,1);
        %albedoSave   = zeros(nAlbedoSave,1);
        %surfPresSave = zeros(nSurfPresSave,1);
        %temperatureCoef = 0;  
        %dAmfSave     = zeros(nPresSave, nSzaSave, nVzaSave, nPhiSave, nAlbedoSave, nSurfPresSave);
        
        presSave = fscanf(fid, '%g', [nPresSave 1]);
        szaSave = fscanf(fid, '%g', [nSzaSave 1]);
        vzaSave = fscanf(fid, '%g', [nVzaSave 1]);
        phiSave = fscanf(fid, '%g', [nPhiSave 1]);
        albedoSave = fscanf(fid, '%g', [nAlbedoSave 1]);
        surfPresSave = fscanf(fid, '%g', [nSurfPresSave 1]);
        
        temperatureCoef = fscanf(fid, '%g', [1 1]);

        [dAmfSave_vec,count] = fscanf(fid, '%g %g', inf);
    
        dAmfSave = reshape(dAmfSave_vec, [nPresSave nSzaSave nVzaSave nPhiSave nAlbedoSave nSurfPresSave]);
    end

    % Now interpolate onto the input values

    nPresSave = numel(presSave);  nSzaSave = numel(szaSave);        nVzaSave = numel(vzaSave);
    nPhiSave = numel(phiSave);    nAlbedoSave = numel(albedoSave);  nSurfPresSave = numel(surfPresSave);
    dAmf1  = zeros(size(phi,1),size(phi,2),nPresSave);
    for i=1:nPresSave;
        dum = reshape(dAmfSave(i,:,:,:,:,:), [nSzaSave, nVzaSave,nPhiSave, nAlbedoSave, nSurfPresSave]);
        %dAmf1(:,i) = interpn(szaSave, vzaSave, phiSave, albedoSave, surfPresSave, dum, sza, vza, phi, albedo, surfPres,'linear');
        dAmf1(:,:,i) = interpn(szaSave, vzaSave, phiSave, albedoSave, surfPresSave, dum, sza, vza, phi, albedo, surfPres,'linear');
    end
    %dAmf = exp(interp1(log(presSave), log(dAmf1'), log(pressure),'linear','extrap'));
    dAmfx=shiftdim(dAmf1,2);
    dAmf = exp(interp1(log(presSave), log(dAmfx), log(pressure),'linear','extrap'));
end

status = fclose(fid);
    