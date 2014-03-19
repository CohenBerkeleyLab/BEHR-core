%%add2grid_5km_new
%%arr 11/19/2009

Dimensions=size(Data(d).ColumnAmountNO2);
swath=d;

x=1:1:Dimensions(1)*Dimensions(2); 
y=1;

Lon1=Data(d).Loncorn(1,x);          Lat1=Data(d).Latcorn(1,x);
if isrow(Lon1); Lon1 = Lon1'; end;  if isrow(Lat1); Lat1 = Lat1'; end;

Lon2=Data(d).Loncorn(2,x);          Lat2=Data(d).Latcorn(2,x);
if isrow(Lon2); Lon2 = Lon2'; end;  if isrow(Lat2); Lat2 = Lat2'; end;

Lon3=Data(d).Loncorn(3,x);          Lat3=Data(d).Latcorn(3,x);
if isrow(Lon3); Lon3 = Lon3'; end;  if isrow(Lat3); Lat3 = Lat3'; end;

Lon4=Data(d).Loncorn(4,x);          Lat4=Data(d).Latcorn(4,x);
if isrow(Lon4); Lon4 = Lon4'; end;  if isrow(Lat4); Lat4 = Lat4'; end;

Lon5=Data(d).Longitude(x);          Lat5=Data(d).Latitude(x);
if isrow(Lon5); Lon5 = Lon5'; end;  if isrow(Lat5); Lat5 = Lat5'; end;

CoordLon=cat(3,Lon1,Lon2,Lon3,Lon4,Lon5); %JLL 18 Mar 2014: cat concatenates along the dimension specified as the first argument
CoordLat=cat(3,Lat1,Lat2,Lat3,Lat4,Lat5);

reslat=resolution; %JLL 18 Mar 2014: resolution, resolution2, lonmin, lonmax, latmin, latmax are all specified in BEHR_main
reslon=resolution2;
lon1=lonmin; lon2=lonmax;  maxy=(abs(lonmin-lonmax))/reslon; miny=1; maxy=single(maxy);
lat1=latmin; lat2=latmax;  maxx=(abs(latmax-latmin))/reslat; minx=1; maxx=single(maxx);

lCoordLon=zeros(Dimensions(1)*Dimensions(2),y,5);
lCoordLat=zeros(Dimensions(1)*Dimensions(2),y,5); %JLL 2-14-2014: This line changed from "lCoordLat = zeros(Dimensions(1)*Dimensions,y,5);" as *Dimensions is not a scalar

for x=1:1:Dimensions(1)*Dimensions(2);
    for c=1:5;
    lCoordLon(x,y,c)=(CoordLon(x,y,c)-lon1)/reslon;
    lCoordLat(x,y,c)=(CoordLat(x,y,c)-lat1)/reslat;
    end
end

hdf_quadrangle_5km_2014
