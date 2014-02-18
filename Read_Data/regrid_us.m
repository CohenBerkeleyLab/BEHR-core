%regrid_us
%arr 11/15/2011

warning off all
addpath C:\Ashley
addpath C:\Ashley\NewCaRetrieval2
addpath C:\Ashley\NewCaRetrieval2\Regrid_tools\
addpath C:\Ashley\NewCaRetrieval2\Regrid_tools\Compute_Corner_Pts_ms
addpath C:\Ashley\Mapping\tools
addpath C:\Ashley\Mapping\mapping\m_map
addpath C:\Ashley\Mapping\mapping\m_map1.4\m_map
tic

%
%*********************************%
lonmin = -125;  lonmax = -65;  
latmin = 25;   latmax = 50;
resolution = 0.05; resolution2 = 0.05;
%*********************************%
%

date_start='2004/10/01';
date_end='2011/09/29';

last_file=dir(fullfile('J:','New_Retrieval_sp_mats','US_wCld'));
if length(last_file(end).name)<5;
    last_date='2004/09/30';
else
    last_date=last_file(end).name(9:16); last_date=([last_date(1:4),'/',last_date(5:6),'/',last_date(7:8)]);
end

total_days=datenum(date_end)-datenum(last_date)+1;
last_date='2011/06/26';
for j=1:total_days;
    R=addtodate(datenum(last_date), j, 'day');
    date=datestr(R,26);
    year=date(1:4);
    month=date(6:7);
    day=date(9:10);
     
    filename = ['OMI_NEW_',year,month,day,'.mat'];

    if exist(['J:\New_Retrieval_sp_mats\NEUS_wCld\BEHR_new\',filename],'file')
        load(['J:\New_Retrieval_sp_mats\NEUS_wCld\BEHR_new\',filename])
        DATA=Data; 
        if exist(['J:\New_Retrieval_sp_mats\SEUS_wCld\BEHR_new\',filename],'file')
            load(['J:\New_Retrieval_sp_mats\SEUS_wCld\BEHR_new\',filename])
            DATA=[DATA,Data];
            if exist(['J:\New_Retrieval_sp_mats\SWUS_wCld\BEHR_new\',filename],'file')
                load(['J:\New_Retrieval_sp_mats\SWUS_wCld\BEHR_new\',filename])
                DATA=[DATA,Data];
                if exist(['J:\New_Retrieval_sp_mats\NEUS_wCld\BEHR_new\',filename],'file')
                    load(['J:\New_Retrieval_sp_mats\NWUS_wCld\BEHR_new\',filename])
                    DATA=[DATA,Data];
        
        %Determine unique swaths 
        u=[]; DATA2=[];
        for s=1:size(DATA,2)
        u=unique([u;DATA(s).Swath]);
        DATA(s).Loncorn=DATA(s).Loncorn(:); DATA(s).Latcorn=DATA(s).Latcorn(:);
        end

        %Set up variable for each swath
        fld={'Longitude';'Latitude';'Loncorn';'Latcorn';'Time';'ViewingZenithAngle';'SolarZenithAngle';'ViewingAzimuthAngle';'SolarAzimuthAngle';'AMFPolluted';'AMFUnpolluted';'CloudFraction';'CloudRadianceFraction';'ColumnAmountNO2Polluted';'TerrainHeight';'TerrainPressure';'TerrainReflectivity';'vcdQualityFlags';'CloudPressure';'SlantColumnAmountNO2';'ColumnAmountNO2';'AMFInitial';'ColumnAmountNO2Initial';'RelativeAzimuthAngle';'ColumnAmountNO2Trop';'TropFractionUnpolluted';'Row';'Swath';'NewCloud';'NewAlbedo';'NewTerpres';'NewAmf';'NewColumnAmountNO2Polluted';'NewColumnAmountNO2Trop'};
        DATA2=struct('Longitude',0,'Latitude',0,'Loncorn',0,'Latcorn',0,'Time',0,'ViewingZenithAngle',0,'SolarZenithAngle',0,'ViewingAzimuthAngle',0,'SolarAzimuthAngle',0,'AMFPolluted',0,'AMFUnpolluted',0,'CloudFraction',0,'CloudRadianceFraction',0,'ColumnAmountNO2Polluted',0,'TerrainHeight',0,'TerrainPressure',0,'TerrainReflectivity',0,'vcdQualityFlags',0,'CloudPressure',0,'SlantColumnAmountNO2',0,'ColumnAmountNO2',0,'AMFInitial',0,'ColumnAmountNO2Initial',0,'RelativeAzimuthAngle',0,'ColumnAmountNO2Trop',0,'TropFractionUnpolluted',0,'Row',0,'Swath',0,'NewCloud',0,'NewAlbedo',0,'NewTerpres',0,'NewAmf',0,'NewColumnAmountNO2Polluted',0,'NewColumnAmountNO2Trop',0);
        for su=1:size(u)-1
            DATA2(su+1)=struct('Longitude',0,'Latitude',0,'Loncorn',0,'Latcorn',0,'Time',0,'ViewingZenithAngle',0,'SolarZenithAngle',0,'ViewingAzimuthAngle',0,'SolarAzimuthAngle',0,'AMFPolluted',0,'AMFUnpolluted',0,'CloudFraction',0,'CloudRadianceFraction',0,'ColumnAmountNO2Polluted',0,'TerrainHeight',0,'TerrainPressure',0,'TerrainReflectivity',0,'vcdQualityFlags',0,'CloudPressure',0,'SlantColumnAmountNO2',0,'ColumnAmountNO2',0,'AMFInitial',0,'ColumnAmountNO2Initial',0,'RelativeAzimuthAngle',0,'ColumnAmountNO2Trop',0,'TropFractionUnpolluted',0,'Row',0,'Swath',0,'NewCloud',0,'NewAlbedo',0,'NewTerpres',0,'NewAmf',0,'NewColumnAmountNO2Polluted',0,'NewColumnAmountNO2Trop',0);
        end
        
        %Combine data into new data variable according to swath number
        for su=1:size(u)
            for s=1:size(DATA,2)
                if DATA(s).Swath(1)==u(su);
                    for f=1:length(fld)
                        eval(strcat(['DATA2(su).',fld{f},'=[DATA2(su).', fld{f},' ;DATA(s).', fld{f},'];']));
                    end
                else
                continue 
                end
            end
        end
        
        %Remove zeros from setting up variable and reorder Latcorn and Loncorn
        for su=1:size(u)
            for f=1:length(fld)
                eval(strcat(['DATA2(su).',fld{f},'(1)=[];']))
            end
            DATA2(su).Loncorn=reshape(DATA2(su).Loncorn,4,numel(DATA2(su).Loncorn)/4);
            DATA2(su).Latcorn=reshape(DATA2(su).Latcorn,4,numel(DATA2(su).Latcorn)/4);
        end

        Data=DATA2;
        clear DATA DATA2 OMI
        
        OMI=[];
        s=size(Data);
        hh=0;
        for d=1:s(2);
            if Data(d).ViewingZenithAngle==0;
            elseif numel(Data(d).ViewingZenithAngle)==1;
            continue
            end
            add2grid_5km_new
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
            Latitude=(25+0.025):0.05:(50-0.025); Latitude=Latitude'; Latitude=repmat(Latitude,1,1200);
            OMI(hh).Latitude=Latitude;
            Longitude=(-125+0.025):0.05:(-65-0.025); Longitude=repmat(Longitude,500,1);
            OMI(hh).Longitude=Longitude;
            OMI(hh).Pixel=Pixel;
            OMI(hh).ColumnAmountNO2Trop=ColumnAmountNO2Trop;
            OMI(hh).TropFractionUnpolluted=TropFractionUnpolluted;
            OMI(hh).NewTerpres=NewTerpres;
            OMI(hh).NewAlbedo=NewAlbedo;
            OMI(hh).NewAmf=NewAmf;
            OMI(hh).NewColumnAmountNO2Trop=NewColumnAmountNO2Trop;
            OMI(hh).NewColumnAmountNO2Polluted=NewColumnAmountNO2Polluted;
            OMI(hh).NewCloud=NewCloud;
            OMI(hh).Row=Row;
            OMI(hh).Swath=Swath;
        end

        savename=['OMI_NEW_',year,month,day];  
        cd('J:\New_Retrieval_sp_mats\US_wCld')
        save(savename,'Data','OMI')
        toc
        t=toc;
                end
            end
        end
    end
end