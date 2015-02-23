%BEHR_hdf2
%arr 7/12/2012

warning off all
tic

date_start='2011/04/01';
date_end='2012/04/29';

%last_file=dir(fullfile('Z:','GROUP','SAT','BEHR','SP_files'));
%if length(last_file(end).name)<5;
%    last_date='2006/04/29';
%else
%    last_date=last_file(end).name(8:15); last_date=([last_date(1:4),'/',last_date(5:6),'/',last_date(7:8)]);
%end

total_days=datenum(date_end)-datenum(date_start)+1;
%last_date='2010/11/03';
for j=1:total_days;
    R=addtodate(datenum(date_start), j, 'day');
    date=datestr(R,26);
    year=date(1:4);
    month=date(6:7);
    day=date(9:10);
    
    directory=(['Z:\GROUP\SAT\BEHR\BEHR_files\NEW\']);
    cd(directory)    
    load(['OMI_BEHR_',year,month,day,'.mat'])
    
    Longitude=[];
    Latitude=[];
    Loncorn1=[];
    Loncorn2=[];
    Loncorn3=[];
    Loncorn4=[];
    Latcorn1=[];
    Latcorn2=[];
    Latcorn3=[];
    Latcorn4=[];
    Time=[];
    ViewingZenithAngle=[];
    SolarZenithAngle=[];
    ViewingAzimuthAngle=[];
    SolarAzimuthAngle=[];
    CloudFraction=[];
    CloudRadianceFraction=[];
    TerrainHeight=[];
    TerrainPressure=[];
    TerrainReflectivity=[];
    vcdQualityFlags=[];
    CloudPressure=[];
    ColumnAmountNO2TropStd=[];
    Row=[];
    Swath=[];
    AMFStrat=[];
    SlantColumnAmountNO2=[];
    ColumnAmountNO2=[];
    AMFTrop=[];
    ColumnAmountNO2Strat=[];
    ColumnAmountNO2Initial=[];
    RelativeAzimuthAngle=[];
    ColumnAmountNO2Trop=[];
    XTrackQualityFlags=[];
    MODISCloud=[];
    MODISAlbedo=[];
    GLOBETerpres=[];
    BEHRAMFTrop=[];
    BEHRColumnAmountNO2Trop=[];
    
    for i=1:length(Data);
        a=squeeze(Data(i).Loncorn(1,:,:)); b=find(a~=0);
        Longitude=[Longitude;Data(i).Longitude(b)];
        Latitude=[Latitude;Data(i).Latitude(b)];
        Loncorn1=[Loncorn1;Data(i).Loncorn(1,b)'];
        Loncorn2=[Loncorn2;Data(i).Loncorn(2,b)'];
        Loncorn3=[Loncorn3;Data(i).Loncorn(3,b)'];
        Loncorn4=[Loncorn4;Data(i).Loncorn(4,b)'];
        Latcorn1=[Latcorn1;Data(i).Latcorn(1,b)'];
        Latcorn2=[Latcorn2;Data(i).Latcorn(2,b)'];
        Latcorn3=[Latcorn3;Data(i).Latcorn(3,b)'];
        Latcorn4=[Latcorn4;Data(i).Latcorn(4,b)'];
        Time=[Time;Data(i).Time(b)];
        ViewingZenithAngle=[ViewingZenithAngle;Data(i).ViewingZenithAngle(b)];
        SolarZenithAngle=[SolarZenithAngle;Data(i).SolarZenithAngle(b)];
        ViewingAzimuthAngle=[ViewingAzimuthAngle;Data(i).ViewingAzimuthAngle(b)];
        SolarAzimuthAngle=[SolarAzimuthAngle;Data(i).SolarAzimuthAngle(b)];
        CloudFraction=[CloudFraction;Data(i).CloudFraction(b)];
        CloudRadianceFraction=[CloudRadianceFraction;Data(i).CloudRadianceFraction(b)];
        TerrainHeight=[TerrainHeight;Data(i).TerrainHeight(b)];
        TerrainPressure=[TerrainPressure;Data(i).TerrainPressure(b)];
        TerrainReflectivity=[TerrainReflectivity;Data(i).TerrainReflectivity(b)];
        vcdQualityFlags=[vcdQualityFlags;Data(i).vcdQualityFlags(b)];
        CloudPressure=[CloudPressure;Data(i).CloudPressure(b)];
        ColumnAmountNO2TropStd=[ColumnAmountNO2TropStd;Data(i).ColumnAmountNO2TropStd(b)];
        Row=[Row;Data(i).Row(b)];
        Swath=[Swath;Data(i).Swath(b)];
        AMFStrat=[AMFStrat;Data(i).AMFStrat(b)];
        SlantColumnAmountNO2=[SlantColumnAmountNO2;Data(i).SlantColumnAmountNO2(b)];
        ColumnAmountNO2=[ColumnAmountNO2;Data(i).ColumnAmountNO2(b)];
        AMFTrop=[AMFTrop;Data(i).AMFTrop(b)];
        ColumnAmountNO2Strat=[ColumnAmountNO2Strat;Data(i).ColumnAmountNO2Strat(b)];
        ColumnAmountNO2Initial=[ColumnAmountNO2Initial;Data(i).ColumnAmountNO2Initial(b)];
        RelativeAzimuthAngle=[RelativeAzimuthAngle;Data(i).RelativeAzimuthAngle(b)];
        ColumnAmountNO2Trop=[ColumnAmountNO2Trop;Data(i).ColumnAmountNO2Trop(b)];
        XTrackQualityFlags=[XTrackQualityFlags;Data(i).XTrackQualityFlags(b)];
        MODISCloud=[MODISCloud;Data(i).MODISCloud(b)];
        MODISAlbedo=[MODISAlbedo;Data(i).MODISAlbedo(b)];
        GLOBETerpres=[GLOBETerpres;Data(i).GLOBETerpres(b)];
        BEHRAMFTrop=[BEHRAMFTrop;Data(i).BEHRAMFTrop(b)];
        BEHRColumnAmountNO2Trop=[BEHRColumnAmountNO2Trop;Data(i).BEHRColumnAmountNO2Trop(b)];
    end
    
    hdf_filename=['OMI_BEHR_',year,month,day,'.hdf'];
    cd Z:\GROUP\SAT\BEHR\BEHR_HDF2 
    hdf5write(hdf_filename,'Data/',[])
    
    sn=length(Data);
    
    for i=1:sn;
        if max(Data(i).Swath)==0;
            continue
        else
            varnms=fieldnames(Data(i));
            %varnms={'CloudFraction';'CloudRadianceFraction';'vcdQualityFlags';'Areaweight';'Latitude';'Longitude';'ColumnAmountNO2Trop';'GLOBETerpres';'MODISAlbedo';'BEHRAMFTrop';'BEHRColumnAmountNO2Trop';'MODISCloud';'Row';'AMFTrop';'XTrackQualityFlags'};
            for k=1:length(varnms);
                local=['/Data/Swath',num2str(max(Data(i).Swath(:))),'/',varnms{k}];
                %var=eval(strcat(['OMI(i).',varnms{k},'(2:end,2:end)']));
                var=eval(strcat(['Data(i).',varnms{k}])); 
                if numel(size(var))==2; 
                    var=single(var'); var(:,end)=[];
                elseif numel(size(var))==3; 
                    var=permute(var,[3 2 1]); var(:,end,:)=[];
                end
                hdf5write(hdf_filename,local,var,'WriteMode','append')
            end
        end
    end
end



   
    
    