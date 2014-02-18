%%regrid_sp_ca
%%arr 5/13/2010

tic

warning off all
addpath C:\Ashley\Mapping\tools
addpath C:\Ashley\Mapping\mapping\m_map
addpath C:\Ashley\Mapping\mapping\m_map1.4\m_map
addpath C:\Ashley\Tools\Regrid_tools

dates = ['060701';'060702';'060703';'060704';'060705';'060706';...
    '060707';'060708';'060709';'060710';'060711';'060712';'060713';...
    '060714';'060715';'060716';'060717';'060718';'060719';'060720';...
    '060721';'060722';'060723';'060724';'060725';'060726';'060727';...
    '060728';'060729';'060730';'060731'];


    %************************************************%
    lonmin = -126;  lonmax = -113;
    latmin = 31;    latmax = 44;
    resolution = 0.05;
    %************************************************%
 
retrieval='SP';
satellite='OMI';

for i=1:length(dates)
    
    year=dates(i,1:2);
    month=dates(i,3:4);
    day=dates(i,5:6);
    
    filename=(['OMI_NEW_20',year,month,day]);
    directory='C:\Ashley\';
    cd(directory)
    a=exist([filename,'.mat'],'file');
    if a==0;
        disp(['No Data Available For ',satellite,' ',month,' ',day,' ',year])
    else
        load(filename)
        s=size(Data);
        hh=0;
        for d=1:s(2);
            if Data(d).ViewingZenithAngle==0;
            continue
            end
            add2grid_5km
           
            hh=hh+1;
            OMI(hh).Time=Time;
            OMI(hh).ViewingZenithAngle=ViewingZenithAngle;
            OMI(hh).SolarZenithAngle=SolarZenithAngle;
            OMI(hh).ViewingAzimuthAngle=ViewingAzimuthAngle;
            OMI(hh).SolarAzimuthAngle=SolarAzimuthAngle;
            OMI(hh).AMFInitial=AMFInitial;
            OMI(hh).AMFPolluted=AMFPolluted;
            OMI(hh).AMFUnpolluted=AMFUnpolluted;
            OMI(hh).CloudFraction=CloudFraction;
            OMI(hh).CloudRadianceFraction=CloudRadianceFraction;
            OMI(hh).ColumnAmountNO2=ColumnAmountNO2;
            OMI(hh).ColumnAmountNO2Initial=ColumnAmountNO2Initial;
            OMI(hh).ColumnAmountNO2Polluted=ColumnAmountNO2Polluted;
            OMI(hh).SlantColumnAmountNO2=SlantColumnAmountNO2;
            OMI(hh).TerrainHeight=TerrainHeight;
            OMI(hh).TerrainPressure=TerrainPressure;
            OMI(hh).TerrainReflectivity=TerrainReflectivity;
            OMI(hh).vcdQualityFlags=vcdQualityFlags;
            OMI(hh).Areaweight=Areaweight;
            OMI(hh).CloudPressure=CloudPressure;
            OMI(hh).RelativeAzimuthAngle=RelativeAzimuthAngle;
                Latitude=(latmin+0.025):0.05:(latmax-0.025); Latitude=Latitude'; Latitude=repmat(Latitude,1,(abs(lonmax-lonmin)/0.05));
            OMI(hh).Latitude=Latitude;
                Longitude=(lonmin+0.025):0.05:(lonmax-0.025); Longitude=repmat(Longitude,((latmax-latmin)/0.05),1);
            OMI(hh).Longitude=Longitude;
            OMI(hh).Pixel=Pixel;
            OMI(hh).ColumnAmountNO2Trop=ColumnAmountNO2Trop;
            OMI(hh).NewAlbedo=NewAlbedo;
            OMI(hh).NewTerpres=NewTerpres;
            OMI(hh).NewAmf=NewAmf;
            OMI(hh).NewColumnAmountNO2Polluted=NewColumnAmountNO2Polluted;
            OMI(hh).NewColumnAmountNO2Trop=NewColumnAmountNO2Trop;
        end
    end
    if exist('OMI','var')==0;
        continue
    else
    cd('C:\Ashley\')
    savename=(['SP_NEW_20',year,month,day]);
    save(savename, 'OMI')
    clear OMI
    toc
    end
end
