function [ corners ] = geometric_corner_calc( lat_in, lon_in )
% Straightforward calculation of pixel corners using their spatial arrangement
% Does not account for satellite IFOV or pixel response; mainly intended for
% datasets that do not include satellite position with their datasets.

earthradius = 6378.5; %km

if ~all(size(lat_in)==size(lon_in)); error('geometric_corner_calc:mismatch','Input lat and lon matrices are not the same size'); end

latcorn = zeros(4,size(lat_in,1),size(lat_in,2));
loncorn = zeros(4,size(lat_in,1),size(lat_in,2));

dims = size(lat_in);
corners = zeros(dims(1),dims(2),2,5);

for y = 1:dims(1)-1;
    lat = lat_in(y,:);
    lon = lon_in(y,:);
    
    lat1 = lat_in(y+1,:);
    lon1 = lon_in(y+1,:);
    
    [latedge, lonedge] = compute_edge_pts(lat, lon);  
    [latedge1, lonedge1] = compute_edge_pts(lat1, lon1);
    
    flightvector = compute_flight_vector(latedge, lonedge, latedge1, lonedge1);
    
    n = length(latedge);
    latcorner = zeros(n,2);
	loncorner = zeros(n,2);
    
    for i = 1:numel(lat)
	    [x1,y1,z1] = sph2cart(lon(i)*pi/180,lat(i)*pi/180,earthradius); %JLL 17 Mar 2014: This time we need the earth's radius because we are interested in the actual physical distance between the points, not in units of earth_radius
    	FV = flightvector(i,:); %JLL 17 Mar 2014: The 3D vector between adjacent rows
    	
    	posx1 = x1 - 0.5 * FV(1); posy1 = y1 - 0.5 * FV(2); posz1 = z1 - 0.5 * FV(3); %JLL 17 Mar 2014: Find the corners for this side of the pixel
   		posx2 = x1 + 0.5 * FV(1); posy2 = y1 + 0.5 * FV(2); posz2 = z1 + 0.5 * FV(3);
    
	    [poslatlonx1, poslatlony1, poslatlonz1] = cart2sph(posx1,posy1,posz1); %JLL 17 Mar 2014: Return to spherical coordinates
    	[poslatlonx2, poslatlony2, poslatlonz2] = cart2sph(posx2,posy2,posz2);
    
	    latcorner(i,1) = poslatlony1 / pi * 180; %JLL 17 Mar 2014: Convert back to degrees
		latcorner(i,2) = poslatlony2 / pi * 180;
	    loncorner(i,1) = poslatlonx1 / pi * 180;
    	loncorner(i,2) = poslatlonx2 / pi * 180;
    end
    
    
    for x = 1:dims(2);
        corners(y,x,1,1) = loncorner(x,1);
        corners(y,x,2,1) = latcorner(x,1);
        corners(y,x,1,2) = loncorner(x+1,1);
        corners(y,x,2,2) = latcorner(x+1,1);
        corners(y,x,1,3) = loncorner(x+1,2);
        corners(y,x,2,3) = latcorner(x+1,2);
        corners(y,x,1,4) = loncorner(x,2);
        corners(y,x,2,4) = latcorner(x,2);
        corners(y,x,1,5) = lon(x);
        corners(y,x,2,5) = lat(x);        
    end   
end

end

