%%rDamf
%%arr 07/23/2008

% Reads the dAMF file used in the FORTRAN OMI NO2 algorithm
% JLL 11 Jul 2014: This version outputs the 6-D matrix contained in the
% damf file, for exploring the AMF under different conditions.

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
function Damf_table = rDamf2_struct(fileDamf)


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

% Save the inputs as a 6-dimensional matrix
Damf_table.matrix = dAmfSave;
Damf_table.dimensions(1).variable = 'Pressure'; Damf_table.dimensions(1).vector = presSave;
Damf_table.dimensions(2).variable = 'SZA'; Damf_table.dimensions(2).vector = szaSave;
Damf_table.dimensions(3).variable = 'VZA'; Damf_table.dimensions(3).vector = vzaSave;
Damf_table.dimensions(4).variable = 'Phi (Relative Azimuth Angle)'; Damf_table.dimensions(4).vector = phiSave;
Damf_table.dimensions(5).variable = 'Albedo'; Damf_table.dimensions(5).vector = albedoSave;
Damf_table.dimensions(6).variable = 'Surface Pressure'; Damf_table.dimensions(6).vector = surfPresSave;




status = fclose(fid);
end
