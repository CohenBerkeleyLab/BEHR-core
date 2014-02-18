%BEHR_SWUS
%arr 5/05/2011
%uses new alb, ter, profs from WRF and MODIS clds

warning off all
tic

addpath K:\Ashley\NewCaRetrieval2\tools
addpath K:\Ashley\NewCaRetrieval2\Regrid_tools
addpath K:\Ashley\Mapping\tools
addpath K:\Ashley\Mapping\mapping\m_map1.4
addpath K:\Ashley\Mapping\mapping\m_map1.4\m_map
addpath J:\ARCTAS_Sat_Data\CA_PROFS
addpath K:\Ashley\ARCTAS
addpath I:\MODIS_myd09cmg\myd09cmg
addpath K:\Ashley\ARCTAS\terrain_database
addpath J:\MODIS_8day
addpath K:\Ashley\Aerosol\cld_mats
addpath K:\Ashley\Trends
addpath K:\Ashley\Tools
addpath J:\New_Retrieval_sp_mats\


fileTmp = 'K:\Ashley\NewCaRetrieval2\tools\nmcTmpYr.txt';
fileDamf='K:\Ashley\NewCaRetrieval2\tools\damf.txt';
fileNO2='K:\Ashley\NewCaRetrieval2\tools\PRFTAV.txt';

load us_4km
satellite='OMI';
retrieval='SP';
region='SWUS';

date_start='2004/10/01';
date_end='2011/06/30';

last_file=dir(fullfile('J:','New_Retrieval_sp_mats','SWUS_wCld','BEHR_new'));
if length(last_file(end).name)<5;
    last_date='2004/09/30';
else
    last_date=last_file(end).name(9:16); last_date=([last_date(1:4),'/',last_date(5:6),'/',last_date(7:8)]);
end

total_days=datenum(date_end,'yyyy/mm/dd')-datenum(last_date,'yyyy/mm/dd')+1;

%last_date='2005/06/10';
for j=1:total_days;
    R=addtodate(datenum(last_date,'yyyy/mm/dd'), j, 'day');
    date=datestr(R,26);
    year=date(1:4);
    month=date(6:7);
    day=date(9:10);
    
    if exist('profile_file')==1 && strcmp(profile_file(2:3),month)==1;
    else
    profile_file=['m',month,'_NO2_profile'];
    load(['K:\Ashley\Trends\US_NO2_Profile\m',month,'_NO2_profile'])
    end 
    
    OMI=[];
    omi=[];

    %cd(['J:\New_Retrieval_sp_mats\',region,'_wCld'])
    filename = ['OMI_SP_',year,month,day,'.mat'];
    if isequal(exist(['J:\New_Retrieval_sp_mats\',region,'_wCld\',filename],'file'),0)
        continue
    else
        load(['J:\New_Retrieval_sp_mats\',region,'_wCld\',filename])
   
        for d=1:length(Data);
            if Data(d).Longitude==0;
            elseif length(Data(d).Longitude)==1; %#ok<ISMT>
                continue 
            else       
                c=length(Data(d).Longitude);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                %get global modis albedo map for specified day
                filename2 = ['mcd43c3_',year,month,day,'.mat'];

                in=[1 -1 2 -2 3 -3 4 -4 5 -5 6 -6 7 -7 8 -8 9 -9 10 -10 11 -11 12 -12 13 -13 14 -14 15 -15 16 -16 17 -17 18 -18 19 -19 20 -20 21 -21];
                for ii=1:length(in);
                    if exist(filename2,'file')==2;
                        %load(filename2,'Albedo_BSA_Band3')
                        continue
                    else
                    date = datestr(addtodate(datenum([month,'/',day,'/',year]),in(ii), 'day'),2);
                    filename2 = ['mcd43c3_20',date(7:8),date(1:2),date(4:5),'.mat'];
                    end
                end

                file2 = ['J:\MODIS_8day\',filename2];
                load(file2,'Albedo_BSA_Band3')

                lat_min=floor(min(Data(d).Latcorn(:)));
                lat_max=ceil(max(Data(d).Latcorn(:)));
                lon_min=floor(min(Data(d).Loncorn(:)));
                lon_max=ceil(max(Data(d).Loncorn(:)));
                    
                %average albedo from MODIS to OMI pixel
                band3=Albedo_BSA_Band3; clear Albedo_BSA_Band3 
                band3_lat=-90+0.05/2:0.05:90-0.05/2; band3_lats=band3_lat'; band3_lats=repmat(band3_lats,1,7200);
                band3_lon=-180+0.05/2:0.05:180-0.05/2; band3_lons=repmat(band3_lon,3600,1);

                ai=find(band3_lat>=lat_min & band3_lat<=lat_max);
                bi=find(band3_lon>=lon_min & band3_lon<=lon_max);
                band3=band3(ai,bi);
                band3_lats=band3_lats(ai,bi);
                band3_lons=band3_lons(ai,bi);

                Data(d).NewAlbedo=zeros(c,1);
                Data(d).NewTerpres=zeros(c,1);

                for j=1:c;

                    x = [];
                    x1 = Data(d).Loncorn(1,j);   y1 = Data(d).Latcorn(1,j);
                    x2 = Data(d).Loncorn(2,j);   y2 = Data(d).Latcorn(2,j);
                    x3 = Data(d).Loncorn(3,j);   y3 = Data(d).Latcorn(3,j);
                    x4 = Data(d).Loncorn(4,j);   y4 = Data(d).Latcorn(4,j);


                    xall=[x1;x2;x3;x4;x1];
                    yall=[y1;y2;y3;y4;y1];

                    xx = inpolygon(band3_lats,band3_lons,yall,xall);

                    band3_vals=band3(xx);  band3_zeros=find(band3_vals==0);
                    band3_vals(band3_zeros)=NaN; band3_vals(isnan(band3_vals))=[];
                    band3_avg=mean(band3_vals);
                    %
                    %put in ocean surface albedo from LUT
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    if isnan(band3_avg)==1;
                        sza_vec = [5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 89];
                        alb_vec = [0.038 0.038 0.039 0.039 0.040 0.042 0.044 0.046 0.051 0.058 0.068 0.082 0.101 0.125 0.149 0.158 0.123 0.073];
                        alb = interp1(sza_vec,alb_vec,Data(d).SolarZenithAngle(j));
                        band3_avg = alb;
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    Data(d).NewAlbedo(j)=band3_avg;

                    clear lat1 lat2 lat3 lat4 xx

                        
                    %get terrain pressure

                    pressure=pres;
                    pressure_lats=pres_lat;
                    pressure_lons=pres_lon;

                    ai=find(pressure_lats(:,1)'>=lat_min & pressure_lats(:,1)'<=lat_max);
                    bi=find(pressure_lons(1,:)>=lon_min & pressure_lons(1,:)<=lon_max);
                    pressure=pressure(ai,bi);
                    pressure_lat=pressure_lats(ai,bi);
                    pressure_lon=pressure_lons(ai,bi);

                    xx = inpolygon(pressure_lat,pressure_lon,yall,xall); 

                    pres_vals=pressure(xx);  pres_zeros=find(pres_vals==0);
                    pres_vals(pres_zeros)=NaN; pres_vals(isnan(pres_vals))=[];
                    Data(d).NewTerpres(j)=mean(pres_vals);
                end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

                Data(d).NewAlbedo(isnan(Data(d).NewAlbedo)==1)=0;
                Data(d).NewTerpres(isnan(Data(d).NewTerpres)==1)=400.0000;

                lon = Data(d).Longitude; 
                lat = Data(d).Latitude;
                mon = str2double(month)*ones(length(Data(d).Latitude),1);
                pressure = [1020 1015 1010 1005 1000 990 980 970 960 945 925 900 875 850 825 800 770 740 700 660 610 560 500 450 400 350 280 200];% 100 50 20 5];
                [temperature, tmpSAVE] = rNmcTmp2(fileTmp, pressure, lon, lat, mon);


                terpres=Data(d).NewTerpres;
                albedo = Data(d).NewAlbedo; 

                sza = Data(d).SolarZenithAngle;
                vza = Data(d).ViewingZenithAngle;  
                phi = Data(d).RelativeAzimuthAngle;


                surfPres = terpres;
                cldPres = Data(d).CloudPressure;
                cldFrac = Data(d).CloudFraction; cldRadFrac = Data(d).CloudRadianceFraction;

                surfPres(surfPres>=1013)=1013;
                    
                %dAmfClr
                [presSave, szaSave, vzaSave, phiSave, albedoSave, surfPresSave, dAmfSave, dAmfClr] = rDamf2(fileDamf, pressure, sza, vza, phi, albedo, surfPres);
                cloudalbedo=0.8*ones(size(Data(d).CloudFraction));
                %dAmfCld
                [presSave, szaSave, vzaSave, phiSave, albedoSave, surfPresSave, dAmfSave, dAmfCld] = rDamf2(fileDamf, pressure, sza, vza, phi, cloudalbedo, cldPres);

                %no2profile
                %[presSAVE, lonSAVE, latSAVE, monSAVE, no2SAVE, no2Profile1] = rno2prof2(fileNO2, pressure, lon, lat, mon);
                %
                loncorns=Data(d).Loncorn;
                latcorns=Data(d).Latcorn;
                [no2_bins] = rProfile_US(PROFILE, loncorns, latcorns, c);
                pp_bins = repmat(pressure,c,1);
                no2_mixratio=no2_bins'./(10^6); % NO2 from WRF is in ppm in these files
                no2Profile1=no2_mixratio;
                clear no2_bins
                %

                aa=find(Data(d).NewAlbedo==0);
                bb=find(Data(d).NewTerpres==400.0000);
                dAmfClr(:,aa)=NaN; 
                dAmfClr(:,bb)=NaN;
                dAmfClr=dAmfClr'; dAmfCld=dAmfCld';
                nans=find(isnan(dAmfClr(:,1)));
                dAmfClr(nans,:)=0;
                dAmfCld(nans,:)=0;
                no2Profile1(nans,:)=0;
                no2Profile2=no2Profile1;
                pTerr=surfPres; pTerr(nans)=0; surfPres(nans)=0;
                pCld = cldPres; pCld(nans)=0;
                cldFrac = Data(d).CloudFraction/1000; cldFrac(nans)=0;
                cldRadFrac = Data(d).CloudRadianceFraction/1000; cldRadFrac(nans)=0;
                cldPres(nans)=0;   vza(nans)=0;    lat(nans)=0;    lon(nans)=0;
                terpres(nans)=0;   mon(nans)=0;    cloudalbedo(nans)=0;   albedo(nans)=0;
                temperature(nans,:)=0;   phi(nans)=0; sza(nans)=0;

                noGhost=1;
                ak=1;
                [amf, amfCld, amfClr, avgKernel, vcd, vcdAvgKernel] = omiAmfAK2(pTerr, pCld, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile1, no2Profile2, noGhost, ak);

                varnms = {'Longitude';'Latitude';'Time';...
                'ViewingZenithAngle';'SolarZenithAngle';'ViewingAzimuthAngle';...
                'SolarAzimuthAngle';'AMFPolluted';'AMFUnpolluted';'CloudFraction';...
                'CloudRadianceFraction';'ColumnAmountNO2Polluted';'TerrainHeight';...
                'TerrainPressure';'TerrainReflectivity';'vcdQualityFlags';'CloudPressure';...
                'SlantColumnAmountNO2';'ColumnAmountNO2';'AMFInitial';...
                'ColumnAmountNO2Initial';'RelativeAzimuthAngle';'ColumnAmountNO2Trop';...
                'NewAlbedo'};
                for i_var=1:length(varnms);
                    if isfield(Data(d),varnms(i_var));
                        eval(strcat(['Data(d).',varnms{i_var},'(nans) = NaN;']));
                    end
                end
                Data(d).Loncorn(:,nans)=NaN;   Data(d).Latcorn(:,nans)=NaN;
                Data(d).NewAmf = amf; 
            end
            OMI=omi;
            %Apply new amf to get new polluted column

            b=length(Data);
            for z=1:b;
                if isfield(Data,'NewAmf')==0 || isempty(Data(z).NewAmf)==1;
                    continue
                else
                    Unpolluted=(Data(z).SlantColumnAmountNO2-Data(z).ColumnAmountNO2Polluted.*Data(z).AMFPolluted)./Data(z).AMFUnpolluted;
                    a=find(Data(z).ColumnAmountNO2Polluted==0);
                    Data(z).NewAmf(a)=0;

                    Data(z).NewColumnAmountNO2Polluted=(Data(z).SlantColumnAmountNO2-Unpolluted.*Data(z).AMFUnpolluted)./Data(z).NewAmf;
                    Data(z).NewColumnAmountNO2Polluted(a)=0;
                    Data(z).NewColumnAmountNO2Polluted(isnan(Data(z).NewColumnAmountNO2Polluted))=0;

                    Data(z).NewColumnAmountNO2Trop=Data(z).NewColumnAmountNO2Polluted+Unpolluted.*Data(z).TropFractionUnpolluted;
                    Data(z).NewColumnAmountNO2Trop(Data(z).NewColumnAmountNO2Trop>1E17)=Unpolluted(Data(z).NewColumnAmountNO2Trop>1E17).*Data(z).TropFractionUnpolluted(Data(z).NewColumnAmountNO2Trop>1E17);
                    Data(z).NewColumnAmountNO2Polluted(Data(z).NewColumnAmountNO2Polluted>1E17)=NaN;
                end
            end
        end

        warning off all
        addpath K:\Ashley
        addpath K:\Ashley\NewCaRetrieval2
        addpath K:\Ashley\NewCaRetrieval2\Regrid_tools\
        addpath K:\Ashley\NewCaRetrieval2\Regrid_tools\Compute_Corner_Pts_ms
        addpath K:\Ashley\Mapping\tools
        addpath K:\Ashley\Mapping\mapping\m_map
        addpath K:\Ashley\Mapping\mapping\m_map1.4\m_map


        %
        %*********************************%
        lonmin = -125;  lonmax = -95;  
        latmin = 25;   latmax = 37.5;
        resolution = 0.05; resolution2 = 0.05;
        %*********************************%
        %

        OMI=[];
        s=size(Data);
        hh=0;
        for d=1:s(2);
            if Data(d).ViewingZenithAngle==0;
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
            Latitude=(25+0.025):0.05:(37.5-0.025); Latitude=Latitude'; Latitude=repmat(Latitude,1,600);
            OMI(hh).Latitude=Latitude;
            Longitude=(-125+0.025):0.05:(-95-0.025); Longitude=repmat(Longitude,250,1);
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
        cd('J:\New_Retrieval_sp_mats\SWUS_wCld\BEHR_new')
        save(savename,'Data','OMI')
        toc
        t=toc;
        if t>2400
            quit force
        end
    end
end
   
