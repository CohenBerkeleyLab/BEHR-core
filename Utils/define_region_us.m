%%define_region_us
%%arr 09/13/2011

addpath C:\Ashley\Mapping\tools
addpath C:\Ashley\Mapping\mapping\m_map1.4
addpath C:\Ashley\Mapping\mapping\m_map1.4\m_map
addpath C:\Ashley\Papers\wdwe_2009\wdwe_8.0

tic  

Latitude=25:0.05:50; Latitudes=Latitude'; Latitudes=repmat(Latitudes,1,1201);
Longitude=-125:0.05:-65; Longitudes=repmat(Longitude,501,1);


%us_region (united states)
us_region_us=ones(size(Latitudes)); 
us_region_us(Latitudes<27)=0;
us_region_us(Latitudes>49)=0;
us_region_us(Longitudes<-68)=0;
cd C:\Ashley\Trends
save us_region_us us_region_us

%sf_region (san francisco)
Distance=zeros(size(Latitudes));
for i=1:numel(Latitudes);
    Distance(i)=m_lldist([-122 Longitudes(i)],[37.6 Latitudes(i)])/1000;
end
sf_region_us=ones(size(Latitudes)); 
sf_region_us(Distance>40)=0;
cd C:\Ashley\Trends
save sf_region_us sf_region_us

%ho_region (houston)
Distance=zeros(size(Latitudes));
for i=1:numel(Latitudes);
    Distance(i)=m_lldist([-95.25 Longitudes(i)],[29.8 Latitudes(i)])/1000;
end
ho_region_us=ones(size(Latitudes)); 
ho_region_us(Distance>30)=0;
cd C:\Ashley\Trends
save ho_region_us ho_region_us


%bkgd_region (united states)
load us_summer05_wd
Data(isnan(Data)==1)=0;
Data(Data>4e15)=0;
bkgd_region_us=Data;
cd C:\Ashley\Trends
save bkgd_region_us bkgd_region_us
