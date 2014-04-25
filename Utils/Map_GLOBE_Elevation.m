%Map_GLOBE_Elevation
%   Maps the elevation within the lat/lon boundaries specified using the
%   GLOBE 1-km terrain database
%
%   Josh Laughner <joshlaugh5@gmail.com> 24 Apr 2014

%Specify the longitude and latitude ranges of interest.
%****************************%
lonmin = -120;    lonmax = -115;
latmin = 32;    latmax = 36;
%****************************%

%Specify the projection type - 'conic' or 'mercator' and whether to include
%the grid lines or not (1 or 0)
%****************************%
proj_type = 'conic';
grid_lines = 0;
%****************************%

globe_dir = '/Volumes/share/GROUP/SAT/BEHR/GLOBE_files';

if ~exist(globe_dir,'dir'); error('globe:not_found','GLOBE Directory not found.  Is the server mounted?'); end

%Go ahead and load the terrain pressure data - only need to do this once
[terpres, refvec] = globedem(globe_dir,1,[latmin, latmax],[lonmin, lonmax]);
    %refvec will contain (1) number of cells per degree, (2)
    %northwest corner latitude, (3) NW corner longitude.
    %(2) & (3) might differ from the input latmin & lonmin
    %because of where the globe cell edges fall

terpres(isnan(terpres)) = -500;
    
cell_count = refvec(1);
globe_latmax = refvec(2); globe_latmin = globe_latmax - size(terpres,1)*(1/cell_count);
globe_lat_matrix = (globe_latmin + 1/(2*cell_count)):(1/cell_count):globe_latmax;
globe_lat_matrix = globe_lat_matrix';
globe_lat_matrix = repmat(globe_lat_matrix,1,size(terpres,2));

globe_lonmin = refvec(3); globe_lonmax = globe_lonmin + size(terpres,2)*(1/cell_count);
globe_lon_matrix = globe_lonmin + 1/(2*cell_count):(1/cell_count):globe_lonmax;
globe_lon_matrix = repmat(globe_lon_matrix,size(terpres,1),1); 

if strcmpi(proj_type,'conic'); m_proj('Albers Equal-Area Conic','lon',[lonmin, lonmax], 'lat', [latmin latmax]);
elseif strcmpi(proj_type, 'mercator'); m_proj('Mercator','long',[lonmin lonmax], 'lat', [latmin latmax]); 
end

m_pcolor(globe_lon_matrix, globe_lat_matrix, terpres);
shading flat;
m_coast('color','w');
m_states('w');

if grid_lines; m_grid;
else m_grid('linestyle','none');
end

cb = colorbar;
ylabel(cb,'Elevation (m)  ','fontsize',12);
title(sprintf('GLOBE elevation for %.2f ^oN %.2f ^oW to %.2f ^oN %.2f ^oW   ', latmin, lonmin, latmax, lonmax),'fontsize',16,'fontweight','bold');