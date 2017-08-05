function [ terpres, globe_lon_matrix, globe_lat_matrix ] = load_globe_alts( lonlim, latlim, scalefactor )
%load_globe_alts Loads GLOBE terrain heights and generates the lat/lon coordinates
%   The GLOBE database is a 30-arcsec resolution terrain height database.
%   The Matlab function "globedem" will load the data, but does not
%   automatically generate the latitude/longitude coordinates for each data
%   point.  This function does so, returning the terrain height in meters,
%   longitude, and latitude as the first, second, and third outputs
%   respectively.
%
%   The third argument (scalefactor) is optional, if omitted will be set to
%   1; i.e. every data point will be returned.  If you need smaller
%   matrices, set this to a higher number, e.g. 4 will only return every
%   4th point.

narginchk(2,3);
if nargin < 3;
    scalefactor = 1;
end

globe_dir = behr_paths.globe_dir;
% Load the terrain data
[terpres, refvec] = globedem(globe_dir,scalefactor,[min(latlim), max(latlim)],[min(lonlim), max(lonlim)]);
    %refvec will contain (1) number of cells per degree, (2)
    %northwest corner latitude, (3) NW corner longitude.
    %(2) & (3) might differ from the input latmin & lonmin
    %because of where the globe cell edges fall
cell_count = refvec(1);
globe_latmax = refvec(2); globe_latmin = globe_latmax - size(terpres,1)*(1/cell_count);
globe_lat_matrix = (globe_latmin + 1/(2*cell_count)):(1/cell_count):globe_latmax;
globe_lat_matrix = globe_lat_matrix';
globe_lat_matrix = repmat(globe_lat_matrix,1,size(terpres,2));

globe_lonmin = refvec(3); globe_lonmax = globe_lonmin + size(terpres,2)*(1/cell_count);
globe_lon_matrix = globe_lonmin + 1/(2*cell_count):(1/cell_count):globe_lonmax;
globe_lon_matrix = repmat(globe_lon_matrix,size(terpres,1),1); 

end

