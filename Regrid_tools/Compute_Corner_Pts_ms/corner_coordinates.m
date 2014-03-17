%%corner_coordinates
%%arr 12/10/2007

%JLL 17 Mar 2014: This script is now deprecated; use
%"fxn_corner_coordinates" instead

%{
Latitude=zeros(60,y); Longitude=zeros(60,y); SpacecraftLatitude=zeros(60,y);
SpacecraftLongitude=zeros(60,y); SpacecraftAltitude=zeros(60,y)
where y = number of scanlines in the orbit, with 60 scenes per scanline
returns array (60,y,2,5), 2=lon,lat, 5=corners(5th=center)
%}

dims = size(Latitude); %JLL 2-18-2014: Latitude should be defined in the main file that calls this script.
corners = zeros(dims(1),dims(2),2,5);

for y = 1:dims(1)-1;
    lat = Latitude(y,:);
    lon = Longitude(y,:);
    
    lat1 = Latitude(y+1,:);
    lon1 = Longitude(y+1,:);
    
    satlat = SpacecraftLatitude(y);
    satlon = SpacecraftLongitude(y);
    satalt = SpacecraftAltitude(y);
    
    satalt = satalt * 1E-3; %convert to km 
    
    %compute the ground pixel edge lat, lon in the across track direction
    %JLL 2-18-2014: These points lie on the edge between pixels in their
    %own rows.
    [latedge, lonedge] = compute_edge_pts(lat, lon);  
    [latedge1, lonedge1] = compute_edge_pts(lat1, lon1);
    
    %compute flight vector
    %JLL 2-19-2014: The ith row is a 3D-vector connecting the ith edge
    %points of two adjacent rows
    flightvector = compute_flight_vector(latedge, lonedge, latedge1, lonedge1);
    
    %compute fwhm in the flight direction if not already calculated
    fwhm = compute_fwhm(latedge, lonedge, satlat, satlon, satalt);
    %compute corner points
    [latcorner, loncorner] = compute_corner_pts(latedge, lonedge, flightvector, fwhm);
    
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

  