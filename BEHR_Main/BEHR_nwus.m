%BEHR_nwus
%arr 2/09/2012
%uses new alb, ter, profs from WRF and MODIS clds

warning off all
tic


addpath C:\Ashley\Mapping\tools
addpath C:\Ashley\Mapping\mapping\m_map1.4
addpath C:\Ashley\Mapping\mapping\m_map1.4\m_map
addpath J:\ARCTAS_Sat_Data\CA_PROFS
addpath C:\Ashley\ARCTAS
addpath I:\MODIS_myd09cmg\myd09cmg
addpath C:\Ashley\ARCTAS\terrain_database
addpath J:\MODIS_8day
addpath C:\Ashley\Aerosol\cld_mats
addpath C:\Ashley\Trends
addpath C:\Ashley\BEHR\AMF_tools
addpath C:\Ashley\Tools\Regrid_tools


fileTmp = 'C:\Ashley\BEHR\AMF_tools\nmcTmpYr.txt';
fileDamf='C:\Ashley\BEHR\AMF_tools\damf.txt';
fileNO2='C:\Ashley\BEHR\AMF_tools\PRFTAV.txt';

%load C:\Ashley\BEHR\US_Terrain_P\us_4km.mat
satellite='OMI';
retrieval='SP';

date_start='2006/04/29';
date_end='2012/04/29';

last_file=dir(fullfile('Z:','GROUP','SAT','BEHR','BEHR_files'));
if length(last_file(end).name)<5;
    last_date='2006/04/29';
else
    last_date=last_file(end).name(10:17); last_date=([last_date(1:4),'/',last_date(5:6),'/',last_date(7:8)]);
end

total_days=datenum(date_end)-datenum(last_date)+1;
%last_date='2004/12/31';
for j=1:total_days;
    R=addtodate(datenum(last_date), j+1, 'day');
    date=datestr(R,26);
    year=date(1:4);
    month=date(6:7);
    day=date(9:10);
    
    if exist('profile_file','file')==1 && strcmp(profile_file(2:3),month)==1;
    else
    profile_file=['m',month,'_NO2_profile'];
    load(['C:\Ashley\BEHR\US_NO2_Profile\m',month,'_NO2_profile'])
    end   
    
    OMI=[];
    omi=[];

    filename = ['OMI_SP_',year,month,day,'.mat'];
    %if isequal(exist(['Z:\GROUP\SAT\BEHR\SP_files\',year,'\',month,'\',filename],'file'),0)
    if isequal(exist(['Z:\GROUP\SAT\BEHR\SP_files\',filename],'file'),0)
    continue
    else
        %load(['Z:\GROUP\SAT\BEHR\SP_files\',year,'\',month,'\',filename])    
        load(['Z:\GROUP\SAT\BEHR\SP_files\',filename])  
        for d=1:length(Data);
            if Data(d).Longitude==0;
            elseif length(Data(d).Longitude)==1; 
                continue 
            else       
                c=numel(Data(d).Longitude);
                %s=size(Data(d).Longitude);
%{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                %get global modis albedo map for specified day
                filename2 = ['mcd43c3_',year,month,day,'.mat'];

                in=[1 -1 2 -2 3 -3 4 -4 5 -5 6 -6 7 -7 8 -8 9 -9 10 -10 11 -11 12 -12 13 -13 14 -14 15 -15 16 -16 17 -17 18 -18 19 -19 20 -20 21 -21];
                for ii=1:length(in);
                    if exist(filename2,'file')==2;
                        continue
                    else
                    date = datestr(addtodate(datenum([month,'/',day,'/',year]),in(ii), 'day'),2);
                    filename2 = ['mcd43c3_20',date(7:8),date(1:2),date(4:5),'.mat'];
                    end
                end

                file2 = ['J:\MODIS_8day\',filename2];
                load(file2,'Albedo_BSA_Band3')

                lat_min=Data(d).Latcorn(:); lat_min(lat_min==0)=[]; lat_min=floor(min(lat_min));
                lat_max=Data(d).Latcorn(:); lat_max(lat_max==0)=[]; lat_max=ceil(max(lat_max));
                lon_min=Data(d).Loncorn(:); lon_min(lon_min==0)=[]; lon_min=floor(min(lon_min));
                lon_max=Data(d).Loncorn(:); lon_max(lon_max==0)=[]; lon_max=ceil(max(lon_max));
                    
                %average albedo from MODIS to OMI pixel
                band3=Albedo_BSA_Band3; clear Albedo_BSA_Band3 
                band3_lat=-90+0.05/2:0.05:90-0.05/2; band3_lats=band3_lat'; band3_lats=repmat(band3_lats,1,7200);
                band3_lon=-180+0.05/2:0.05:180-0.05/2; band3_lons=repmat(band3_lon,3600,1);

                %ai=find(band3_lat>=lat_min & band3_lat<=lat_max);
                %bi=find(band3_lon>=lon_min & band3_lon<=lon_max);
                %band3=band3(ai,bi);
                %band3_lats=band3_lats(ai,bi);
                %band3_lons=band3_lons(ai,bi);
                
                s=size(Data(d).Latitude);
                MODISAlbedo=zeros(s);
                GLOBETerpres=zeros(s);
               tic
                for k=1:c;

                    x = [];
                    x1 = Data(d).Loncorn(1,k);   y1 = Data(d).Latcorn(1,k);
                    x2 = Data(d).Loncorn(2,k);   y2 = Data(d).Latcorn(2,k);
                    x3 = Data(d).Loncorn(3,k);   y3 = Data(d).Latcorn(3,k);
                    x4 = Data(d).Loncorn(4,k);   y4 = Data(d).Latcorn(4,k);


                    xall=[x1;x2;x3;x4;x1];
                    yall=[y1;y2;y3;y4;y1];
                    
                    %%%%SPEED IT UP%%%%
                    ai=find(band3_lat>=min(yall) & band3_lat<=max(yall));
                    bi=find(band3_lon>=min(xall) & band3_lon<=max(xall));
                    band3x=band3(ai,bi);
                    band3_latsx=band3_lats(ai,bi);
                    band3_lonsx=band3_lons(ai,bi);
                    %%%%%%%%%%%%%%%%%%%
                    xx = inpolygon(band3_latsx,band3_lonsx,yall,xall);

                    band3_vals=band3x(xx);  band3_zeros=find(band3_vals==0);
                    band3_vals(band3_zeros)=NaN; band3_vals(isnan(band3_vals))=[];
                    band3_avg=mean(band3_vals);
                    
                    %put in ocean surface albedo from LUT
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if isnan(band3_avg)==1;
                        sza_vec = [5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 89];
                        alb_vec = [0.038 0.038 0.039 0.039 0.040 0.042 0.044 0.046 0.051 0.058 0.068 0.082 0.101 0.125 0.149 0.158 0.123 0.073];
                        alb = interp1(sza_vec,alb_vec,Data(d).SolarZenithAngle(k));
                        band3_avg = alb;
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    
                    MODISAlbedo(k)=band3_avg;

                    clear lat1 lat2 lat3 lat4 xx

                        
                    %get terrain pressure

                    pressure=pres;
                    pressure_lats=pres_lat;
                    pressure_lons=pres_lon;

                   % ai=find(pressure_lats(:,1)'>=lat_min & pressure_lats(:,1)'<=lat_max);
                   % bi=find(pressure_lons(1,:)>=lon_min & pressure_lons(1,:)<=lon_max);
                   % pressure=pressure(ai,bi);
                   % pressure_lat=pressure_lats(ai,bi);
                   % pressure_lon=pressure_lons(ai,bi);
                    
                    %%%%SPEED IT UP%%%%
                    ai=find(pressure_lats(:,1)>=min(yall) & pressure_lats(:,1)<=max(yall));
                    bi=find(pressure_lons(1,:)>=min(xall) & pressure_lons(1,:)<=max(xall));
                    pressurex=pressure(ai,bi);
                    pressure_latx=pressure_lats(ai,bi);
                    pressure_lonx=pressure_lons(ai,bi);
                    %%%%%%%%%%%%%%%%%%%

                    xx = inpolygon(pressure_latx,pressure_lonx,yall,xall); 

                    pres_vals=pressurex(xx);  pres_zeros=find(pres_vals==0);
                    pres_vals(pres_zeros)=NaN; pres_vals(isnan(pres_vals))=[];
                    GLOBETerpres(k)=mean(pres_vals);
                end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%}
                Data(d).MODISAlbedo(isnan(Data(d).MODISAlbedo)==1)=0;
                Data(d).GLOBETerpres(isnan(Data(d).GLOBETerpres)==1)=400.0000;

                lon = Data(d).Longitude; 
                lat = Data(d).Latitude;
                mon = str2double(month)*ones(size(Data(d).Latitude));
                pressure = [1020 1015 1010 1005 1000 990 980 970 960 945 925 900 875 850 825 800 770 740 700 660 610 560 500 450 400 350 280 200];% 100 50 20 5];
                [temperature, tmpSAVE] = rNmcTmp2(fileTmp, pressure, lon, lat, mon);


                terpres = Data(d).GLOBETerpres;
                albedo = Data(d).MODISAlbedo; 

                sza = Data(d).SolarZenithAngle;
                vza = Data(d).ViewingZenithAngle;  
                phi = Data(d).RelativeAzimuthAngle;


                surfPres = terpres;
                cldPres = Data(d).CloudPressure;
                %cldFrac = Data(d).CloudFraction; cldRadFrac = Data(d).CloudRadianceFraction;

                surfPres(surfPres>=1013)=1013;
                    
                %dAmfClr
                %[presSave, szaSave, vzaSave, phiSave, albedoSave, surfPresSave, dAmfSave, dAmfClr] = rDamf2(fileDamf, pressure, sza, vza, phi, albedo, surfPres);
                dAmfClr = rDamf2(fileDamf, pressure, sza, vza, phi, albedo, surfPres);
                cloudalbedo=0.8*ones(size(Data(d).CloudFraction));
                %dAmfCld
                %[presSave, szaSave, vzaSave, phiSave, albedoSave, surfPresSave, dAmfSave, dAmfCld] = rDamf2(fileDamf, pressure, sza, vza, phi, cloudalbedo, cldPres);
                dAmfCld = rDamf2(fileDamf, pressure, sza, vza, phi, cloudalbedo, cldPres);
                
                %no2profile
                %[presSAVE, lonSAVE, latSAVE, monSAVE, no2SAVE, no2Profile1] = rno2prof2(fileNO2, pressure, lon, lat, mon);
                %
                loncorns=Data(d).Loncorn;
                latcorns=Data(d).Latcorn;
                [no2_bins] = rProfile_US(PROFILE, loncorns, latcorns, c);
                no2_bins=reshape(no2_bins,length(pressure),size(vza,1),size(vza,2));
                pp_bins = repmat(pressure,c,1);
                no2_mixratio=no2_bins./(10^6); % NO2 from WRF is in ppm in these files
                no2Profile1=no2_mixratio;
                prof_i=zeros(size(Data(d).Latitude)); prof_i(isnan(squeeze(no2Profile1(1,:,:)))==1)=1;
                clear no2_bins
                %

                aa=find(Data(d).MODISAlbedo==0);
                bb=find(Data(d).GLOBETerpres==400.0000);
                %dAmfClr(:,aa)=NaN; 
                %dAmfClr(:,bb)=NaN;
                %dAmfClr=dAmfClr'; dAmfCld=dAmfCld';
                %nans=find(isnan(dAmfClr(:,1)));
                %dAmfClr(nans,:)=0;
                %dAmfCld(nans,:)=0;
                %no2Profile1(nans,:)=0;
                no2Profile2=no2Profile1;
                pTerr=surfPres; %pTerr(nans)=0; surfPres(nans)=0;
                pCld = cldPres; %pCld(nans)=0;
                cldFrac = Data(d).CloudFraction/1000; %cldFrac(nans)=0;
                cldRadFrac = Data(d).CloudRadianceFraction/1000; %cldRadFrac(nans)=0;
                %cldPres(nans)=0;   vza(nans)=0;    lat(nans)=0;    lon(nans)=0;
                %terpres(nans)=0;   mon(nans)=0;    cloudalbedo(nans)=0;   albedo(nans)=0;
                %temperature(nans,:)=0;   phi(nans)=0; sza(nans)=0;

                noGhost=1; ak=1;
                [amf, amfCld, amfClr] = omiAmfAK2(pTerr, pCld, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile1, no2Profile2, noGhost, ak);
                amf(prof_i==1)=NaN;
                
                %Data(d).MODISAlbedo=MODISAlbedo;
                %Data(d).GLOBETerpres=GLOBETerpres;
                
                varnms = {'Longitude';'Latitude';'Time';...
                'ViewingZenithAngle';'SolarZenithAngle';'ViewingAzimuthAngle';...
                'SolarAzimuthAngle';'AMFPolluted';'AMFUnpolluted';'CloudFraction';...
                'CloudRadianceFraction';'ColumnAmountNO2Polluted';'TerrainHeight';...
                'TerrainPressure';'TerrainReflectivity';'vcdQualityFlags';'CloudPressure';...
                'SlantColumnAmountNO2';'ColumnAmountNO2';'AMFInitial';...
                'ColumnAmountNO2Initial';'RelativeAzimuthAngle';'ColumnAmountNO2Trop';...
                'MODISAlbedo';'GLOBETerpres'};
                %for i_var=1:length(varnms);
                %    if isfield(Data(d),varnms(i_var));
                %        eval(strcat(['Data(d).',varnms{i_var},'(nans) = NaN;']));
                %    end
                %end
                %Data(d).Loncorn(:,nans)=NaN;   Data(d).Latcorn(:,nans)=NaN;
                Data(d).BEHRAMFTrop = amf;            
            end

            %Apply new amf to get new polluted column

            b=length(Data);
            for z=1:b;
                if isfield(Data,'BEHRAMFTrop')==0 || isempty(Data(z).BEHRAMFTrop)==1;
                    continue
                else
                    Data(z).BEHRColumnAmountNO2Trop=Data(z).ColumnAmountNO2Trop.*Data(z).AMFTrop./Data(z).BEHRAMFTrop;
                end
            end
        end

        warning off all
        addpath C:\Ashley
        addpath C:\Ashley\Tools
        addpath C:\Ashley\Tools\Regrid_tools\
        addpath C:\Ashley\Tools\Regrid_tools\Compute_Corner_Pts_ms
        addpath C:\Ashley\Mapping\tools
        addpath C:\Ashley\Mapping\mapping\m_map
        addpath C:\Ashley\Mapping\mapping\m_map1.4\m_map


        %
        %*********************************%
        lonmin = -125;  lonmax = -65;  
        latmin = 25;   latmax = 50;
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
            Latitude=(25+0.025):0.05:(50-0.025); Latitude=Latitude'; Latitude=repmat(Latitude,1,1200);
            OMI(hh).Latitude=Latitude;
            Longitude=(-125+0.025):0.05:(-65-0.025); Longitude=repmat(Longitude,500,1);
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

        savename=['OMI_BEHR_',year,month,day];  
        cd('Z:\GROUP\SAT\BEHR\BEHR_files\');%',year,'\',month,'\'])
        save(savename,'Data','OMI')
        toc
        t=toc;
        if t>900
            quit force
        end
    end
end

  