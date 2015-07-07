%%rProfile_US
%%arr 04/22/2011
%Reads monthly WRF profiles from Lukas (for new retrieval, trends)
%NO2 is in ppm

function [no2_bins] = rProfile_US(PROFILE, loncorns, latcorns, c)

pressure = PROFILE.Pressure;
lat_prs  = PROFILE.Latitude;
lon_prs  = PROFILE.Longitude;

sz = size(PROFILE.NO2_profile);
%no2=reshape(PROFILE.NO2_profile,length(pressure),194300);
no2=reshape(PROFILE.NO2_profile,length(pressure),sz(2)*sz(3)); % reshape to a 2D matrix
no2_bins=zeros(length(pressure),c);
for j=1:c;
                    
    x = [];
    x1 = loncorns(1,j);   y1 = latcorns(1,j);
    x2 = loncorns(2,j);   y2 = latcorns(2,j);
    x3 = loncorns(3,j);   y3 = latcorns(3,j);
    x4 = loncorns(4,j);   y4 = latcorns(4,j);
    
    xall=[x1;x2;x3;x4;x1];
    yall=[y1;y2;y3;y4;y1];
    
    xx = inpolygon(lat_prs,lon_prs,yall,xall); 
    
    no2_bin=no2(:,xx);
    no2_bins(:,j)=mean(no2_bin,2);
end


%for pres_i=1:length(pressure)
%    no2_bins(:,pres_i)=(interp2(lon_prs, lat_prs, squeeze(PROFILE.NO2new(pres_i,:,:)), lon, lat,'linear')); 
%end
%{
old version of this file

function [no2_bins] = rProfile(loncorns, latcorns, c)

if exist('PROFILE','var')==0
    load PROFILEalt
end

pressure = [1020 1015 1010 1005 1000 990 980 970 960 945 925 900 875 850 825 800 770 740 700 660 610 560 500 450 400 350 280 200 100 50 20 5];
lat_prs  = 30.025:0.05:49.975;       lat_prs=lat_prs';  lat_prs=repmat(lat_prs,1,480);  nlat=numel(lat_prs);
lon_prs  = -123.975:0.05:-100.025;                      lon_prs=repmat(lon_prs,400,1);  nlon=numel(lon_prs);

for j=1:c;
                    
    %if mod(j,50)==0
    %    disp([num2str(j),' lines read so far out of ', num2str(c)])
    %end
    x = [];
    x1 = loncorns(1,j);   y1 = latcorns(1,j);
    x2 = loncorns(2,j);   y2 = latcorns(2,j);
    x3 = loncorns(3,j);   y3 = latcorns(3,j);
    x4 = loncorns(4,j);   y4 = latcorns(4,j);
    m1 = (y2-y1)/(x2-x1);   m2 = (y3-y2)/(x3-x2);
    m3 = (y4-y3)/(x4-x3);   m4 = (y1-y4)/(x1-x4);


    lat1 = (m1*(lon_prs-x1)+y1);     lat2 = (m2*(lon_prs-x2)+y2);
    lat3 = (m3*(lon_prs-x3)+y3);     lat4 = (m4*(lon_prs-x4)+y4);

    xx = find(lat_prs>lat1 & lat_prs<lat2 & lat_prs<lat3 & lat_prs>lat4); 
    
    no2=PROFILE.NO2new(:,:);
    no2_bin=no2(:,xx);
    no2_bins(:,j)=mean(no2_bin,2);
end

%}