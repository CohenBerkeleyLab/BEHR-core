function [ alb_mat, alb_lon_mat, alb_lat_mat ] = read_modis_mcd43c3_2014( filename )
%read_modis_mcd43c3_2014 Version of read_modis_mcd43c3 updated for the 2014 format of mcd43 files (minimal change) and altered to be a function.
%   This function requires the filename of the mcd43 file to be read as a
%   string; if the file is not on your MATLAB path, this must be the full
%   absolute path to the file.  It will return three matrices: the matrix
%   of albedo values in band 3 (scaled to their proper values), 
%   and the longitude and latitudes of the grid cell centers.
%
%   Josh Laughner <joshlaugh5@gmail.com> 30 Apr 2014

mcd43_info = hdfinfo(filename);
alb_mat = hdfread(hdf_dsetID(mcd43_info,1,1,'Albedo_BSA_Band3'));
alb_mat = double(alb_mat);
alb_mat = flipud(alb_mat);
alb_mat = alb_mat .* 1e-3;

%MODIS albedo is given in 0.05 degree cells and a single file covers the
%full globe, so figure out the lat/lon of the middle of the grid cells as:
alb_lat_vec=-90+0.05/2:0.05:90-0.05/2; 
alb_lon_vec=-180+0.05/2:0.05:180-0.05/2;

[alb_lon_mat, alb_lat_mat] = meshgrid(alb_lon_vec, alb_lat_vec);

end

