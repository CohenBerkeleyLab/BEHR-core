%BEHR_hdf
%arr 7/6/2012

warning off all
addpath C:\Ashley
addpath C:\Ashley\Tools
addpath C:\Ashley\Tools\Regrid_tools\
addpath C:\Ashley\Tools\Regrid_tools\Compute_Corner_Pts_ms
addpath C:\Ashley\Mapping\tools
addpath C:\Ashley\Mapping\mapping\m_map
addpath C:\Ashley\Mapping\mapping\m_map1.4\m_map


date_start='2004/12/31';
date_end='2012/04/30';

%last_file=dir(fullfile('Z:','GROUP','SAT','BEHR','BEHR_HDF'));
%if length(last_file(end).name)<10;
%    last_date='2004/12/31';
%else
%    last_date=last_file(end).name(10:17); last_date=([last_date(1:4),'/',last_date(5:6),'/',last_date(7:8)]);
%end

total_days=datenum(date_end)-datenum(date_start)+1;
%last_date='2004/12/31';

for j=1:total_days;
    R=addtodate(datenum(date_start), j, 'day');
    date=datestr(R,26);
    year=date(1:4);
    month=date(6:7);
    day=date(9:10);
    
    mat_filename=['OMI_BEHR_',year,month,day];
    load(['Z:\GROUP\SAT\BEHR\BEHR_files\NEWER\',mat_filename],'Data','OMI')
    
    %{
    %
    %*********************************%
    lonmin = -125.05;  lonmax = -65;  
    latmin = 25;   latmax = 50.05;
    resolution = 0.05; resolution2 = 0.05;
    %*********************************%
    %

    OMI=struct('Time',0,'ViewingZenithAngle',0,'SolarZenithAngle',0,'ViewingAzimuthAngle',0,'SolarAzimuthAngle',0,'CloudFraction',0,'CloudRadianceFraction',0,'ColumnAmountNO2',0,'SlantColumnAmountNO2',0,'TerrainHeight',0,'TerrainPressure',0,'TerrainReflectivity',0,'vcdQualityFlags',0,'Areaweight',0,'CloudPressure',0,'RelativeAzimuthAngle',0,'Latitude',0,'Longitude',0,'Pixel',0,'ColumnAmountNO2Trop',0,'GLOBETerpres',0,'MODISAlbedo',0,'BEHRAMFTrop',0,'BEHRColumnAmountNO2Trop',0,'MODISCloud',0,'Row',0,'Swath',0,'AMFTrop',0,'AMFStrat',0,'ColumnAmountNO2Strat',0,'ColumnAmountNO2Initial',0,'XTrackQualityFlags',0);
    s=size(Data);
    hh=0;
    for d=1:s(2);
        if Data(d).ViewingZenithAngle==0;
        elseif numel(Data(d).ViewingZenithAngle)==1;
        continue
        else
        add2grid_5km_new
        hh=hh+1;
        OMI(hh).Time=Time;
        OMI(hh).ViewingZenithAngle=ViewingZenithAngle;
        OMI(hh).SolarZenithAngle=SolarZenithAngle;
        OMI(hh).ViewingAzimuthAngle=ViewingAzimuthAngle;
        OMI(hh).SolarAzimuthAngle=SolarAzimuthAngle;
        OMI(hh).CloudFraction=CloudFraction;
        OMI(hh).CloudRadianceFraction=CloudRadianceFraction;
        OMI(hh).ColumnAmountNO2=ColumnAmountNO2;
        OMI(hh).SlantColumnAmountNO2=SlantColumnAmountNO2;
        OMI(hh).TerrainHeight=TerrainHeight;
        OMI(hh).TerrainPressure=TerrainPressure;
        OMI(hh).TerrainReflectivity=TerrainReflectivity;
        OMI(hh).vcdQualityFlags=vcdQualityFlags;
        OMI(hh).Areaweight=Areaweight;
        OMI(hh).CloudPressure=CloudPressure;
        OMI(hh).RelativeAzimuthAngle=RelativeAzimuthAngle;
        Latitude=(25+0.025):0.05:(50+0.025); Latitude=Latitude'; Latitude=repmat(Latitude,1,1201);
        OMI(hh).Latitude=Latitude;
        Longitude=(-125-0.025):0.05:(-65-0.025); Longitude=repmat(Longitude,501,1);
        OMI(hh).Longitude=Longitude;
        OMI(hh).Pixel=Pixel;
        OMI(hh).ColumnAmountNO2Trop=ColumnAmountNO2Trop;
        OMI(hh).GLOBETerpres=GLOBETerpres;
        OMI(hh).MODISAlbedo=MODISAlbedo;
        OMI(hh).BEHRAMFTrop=BEHRAMFTrop;
        OMI(hh).BEHRColumnAmountNO2Trop=BEHRColumnAmountNO2Trop;
        OMI(hh).MODISCloud=MODISCloud;
        OMI(hh).Row=Row;
        OMI(hh).Swath=Swath;
        OMI(hh).AMFTrop=AMFTrop;
        OMI(hh).AMFStrat=AMFStrat;
        OMI(hh).ColumnAmountNO2Strat=ColumnAmountNO2Strat;
        OMI(hh).ColumnAmountNO2Initial=ColumnAmountNO2Initial;
        OMI(hh).XTrackQualityFlags=XTrackQualityFlags;
        end
    end
    %}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    
    hdf_filename=['OMI_BEHR_',year,month,day,'.hdf'];
    cd Z:\GROUP\SAT\BEHR\BEHR_REGRIDDED_HDF2 
    hdf5write(hdf_filename,'Data/',[])
    
    sn=length(OMI);
    
    for i=1:sn;
        if max(OMI(i).Swath)==0;
            continue
        else
            %varnms=fieldnames(OMI(i));
            varnms={'CloudFraction';'CloudRadianceFraction';'vcdQualityFlags';'Areaweight';'Latitude';'Longitude';'ColumnAmountNO2Trop';'GLOBETerpres';'MODISAlbedo';'BEHRAMFTrop';'BEHRColumnAmountNO2Trop';'MODISCloud';'Row';'AMFTrop';'XTrackQualityFlags';'ColumnAmountNO2TropStd'};
            for k=1:length(varnms);
                local=['/Data/Swath',num2str(max(OMI(i).Swath(:))),'/',varnms{k}];
                %var=eval(strcat(['OMI(i).',varnms{k},'(2:end,2:end)']));
                var=single(eval(strcat(['OMI(i).',varnms{k}])));
                hdf5write(hdf_filename,local,var,'WriteMode','append')
            end
        end
    end
end
    