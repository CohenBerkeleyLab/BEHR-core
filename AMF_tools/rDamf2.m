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
end

mat_test = [ismatrix(sza), ismatrix(vza), ismatrix(phi), ismatrix(albedo), ismatrix(surfPres)];
mat_error = {'SZA', 'VZA', 'PHI', 'ALBEDO', 'SURFPRES'};
if any(~mat_test)
    E.badinput('The inputs %s must be 2D (matrices or vectors). The following are not: %s', strjoin(mat_error, ', '), strjoin(mat_error(~mat_test), ', '));
end
sz = size(sza);
size_test = [isequal(sz, size(vza)), isequal(sz, size(phi)), isequal(sz, size(albedo)), isequal(sz, size(albedo))];
if any(~size_test)
    E.badinput('The inputs %s must be all the same size. The following are not the same size as SZA: %s', strjoin(mat_error, ', '), strjoin(mat_error(find(~size_test)+1), ', '));
end


if exist('dAmfSave','var')==0;
    
    fid = fopen(fileDamf,'r');
    
    header = fgets(fid);
    
    n = fscanf(fid,'%g',[1 6]);
    
    nPresSave = n(1); nSzaSave = n(2); nVzaSave = n(3); nPhiSave = n(4); nAlbedoSave = n(5); nSurfPresSave = n(6);
    
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
    dAmf1(:,:,i) = interpn(szaSave, vzaSave, phiSave, albedoSave, surfPresSave, dum, sza, vza, phi, albedo, surfPres,'linear');
end

dAmfx=shiftdim(dAmf1,2);
dAmf = exp(interp1(log(presSave), log(dAmfx), log(pressure),'linear','extrap'));

status = fclose(fid);
end
