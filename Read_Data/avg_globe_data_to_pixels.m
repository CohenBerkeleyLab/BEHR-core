function [ data ] = avg_globe_data_to_pixels( data, globe_elevations, globe_lon_matrix, globe_lat_matrix, varargin )
%AVG_GLOBE_DATA_TO_PIXELS Averages GLOBE elevation data to OMI pixels
%   DATA = AVG_GLOBE_DATA_TO_PIXELS( DATA, GLOBE_ELEV, GLOBE_LON,
%   GLOBE_LAT) Averages Global Land One-km Elevation (GLOBE) data given in
%   GLOBE_ELEV at coordinates GLOBE_LON and GLOBE_LAT to the pixels with
%   corners given by the fields FoV75CornerLongitude, FoV75CornerLatitude
%   in the structure DATA. GLOBE_LON and GLOBE_LAT must be matrices the
%   same shape as GLOBE_ELEV. The elevations are converted to pressure by
%   the equation: 
%
%       p = (1013 hPa) * exp( -z / 7400 )
%
%   where z is the GLOBE elevation in meters. This is stored as the field
%   GLOBETerpres in DATA.
%
%   Parameters:
%
%       'DEBUG_LEVEL' - increase the verbosity. Default is 0, higher
%       numbers print more information.
%
%       'LoncornField', 'LatcornField' - change which fields in DATA are
%       used as the definition of the pixel corners


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
p = inputParser;
p.addParameter('DEBUG_LEVEL', 0);
p.addParameter('LoncornField', 'FoV75CornerLongitude');
p.addParameter('LatcornField', 'FoV75CornerLatitude');
p.parse(varargin{:});
pout = p.Results;

DEBUG_LEVEL = pout.DEBUG_LEVEL;
loncorn_field = pout.LoncornField;
latcorn_field = pout.LatcornField;

GLOBETerrainHeight = nan(size(data.Latitude));

%GLOBE matrices are arrange s.t. terpres(1,1) is in the SW
%corner and terpres(end, end) is in the NE corner.

for a=1:numel(GLOBETerrainHeight)
    
    if DEBUG_LEVEL > 3; fprintf('Averaging GLOBE data to pixel %u of %u \n',a,numel(GLOBETerrainHeight)); end
    if DEBUG_LEVEL > 3; tic; end
    
    xall=[data.(loncorn_field)(:,a); data.(loncorn_field)(1,a)];
    yall=[data.(latcorn_field)(:,a); data.(latcorn_field)(1,a)];
    
    % If there is an invalid corner coordinate, skip because we cannot
    % be sure the correct polygon will be used.
    if any(isnan(xall)) || any(isnan(yall))
        continue
    end
    
    %%%%SPEED IT UP%%%%
    % Since GLOBE data is on a grid where a row of
    % latitudinal points all have the same longitude and
    % vice versa, we can quickly reduce the number of
    % points by comparing just one lat and lon vector to
    % the extent of the pixel.
    if DEBUG_LEVEL > 4; t_cut = tic; end
    ai=find(globe_lat_matrix(:,1)>=min(yall) & globe_lat_matrix(:,1)<=max(yall));
    bi=find(globe_lon_matrix(1,:)>=min(xall) & globe_lon_matrix(1,:)<=max(xall));
    elevation_x=globe_elevations(ai,bi);
    pressure_latx=globe_lat_matrix(ai,bi);
    pressure_lonx=globe_lon_matrix(ai,bi);
    if DEBUG_LEVEL > 4; fprintf('    Time to cut down GLOBE data = %f\n', toc(t_cut)); end
    %%%%%%%%%%%%%%%%%%%
    
    % inpolygon is slow compared to a simple logical test,
    % so we only apply it to the subset of GLOBE heights
    % immediately around our pixel.
    if DEBUG_LEVEL > 4; t_poly = tic; end
    xx_globe = inpolygon(pressure_latx,pressure_lonx,yall,xall);
    if DEBUG_LEVEL > 4; fprintf('    Time to apply inpolygon to GLOBE data = %f\n', toc(t_poly)); end
    
    elev_vals=elevation_x(xx_globe);
    GLOBETerrainHeight(a)=mean(elev_vals);
    
    
    if DEBUG_LEVEL > 3; telap = toc; fprintf('Time for GLOBE --> pixel %u/%u = %g sec \n',a,numel(GLOBETerrainHeight),telap); end
end

data.GLOBETerrainHeight = GLOBETerrainHeight;
end

