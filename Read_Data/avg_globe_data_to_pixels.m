function [ data ] = avg_globe_data_to_pixels( data, globe_elevations, globe_lon_matrix, globe_lat_matrix, varargin )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here


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

GLOBETerpres = zeros(size(data.Latitude));

%GLOBE matrices are arrange s.t. terpres(1,1) is in the SW
%corner and terpres(end, end) is in the NE corner.

for a=1:numel(GLOBETerpres)
    
    if DEBUG_LEVEL > 1; fprintf('Averaging GLOBE data to pixel %u of %u \n',a,c); end
    if DEBUG_LEVEL > 2; tic; end
    
    xall=[data.(loncorn_field)(:,a); data.(loncorn_field)(1,a)];
    yall=[data.(latcorn_field)(:,a); data.(latcorn_field)(1,a)];
    
    %%%%SPEED IT UP%%%%
    % Since GLOBE data is on a grid where a row of
    % latitudinal points all have the same longitude and
    % vice versa, we can quickly reduce the number of
    % points by comparing just one lat and lon vector to
    % the extent of the pixel.
    ai=find(globe_lat_matrix(:,1)>=min(yall) & globe_lat_matrix(:,1)<=max(yall));
    bi=find(globe_lon_matrix(1,:)>=min(xall) & globe_lon_matrix(1,:)<=max(xall));
    pressurex=globe_elevations(ai,bi);
    pressure_latx=globe_lat_matrix(ai,bi);
    pressure_lonx=globe_lon_matrix(ai,bi);
    %%%%%%%%%%%%%%%%%%%
    
    % inpolygon is slow compared to a simple logical test,
    % so we only apply it to the subset of GLOBE heights
    % immediately around our pixel.
    xx_globe = inpolygon(pressure_latx,pressure_lonx,yall,xall);
    
    pres_vals=pressurex(xx_globe);
    GLOBETerpres(a)=1013.25 .* exp(-mean(pres_vals) / 7400 ); %Originally divided by 7640 m
    
    
    if DEBUG_LEVEL > 2; telap = toc; fprintf('Time for GLOBE --> pixel %u/%u = %g sec \n',a,c,telap); end
end

data.GLOBETerpres = GLOBETerpres;
end

