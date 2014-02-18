%%readhdf_modis_cld_us
%%arr 04/22/2011

addpath C:\Ashley\NewCaRetrieval2\Regrid_tools
addpath C:\Ashley\NewCaRetrieval2\Regrid_tools\Compute_Corner_Pts_ms

warning off all
tic

%****************************%
lonmin = -126;  lonmax = -64;
latmin = 24;    latmax = 51;
%****************************%

satellite='MODIS';
retrieval='cld';

%for jj=1:length(years)   
    year='2004';%years{jj};
    days=275:366;   
    for i=1;%:length(days);
        d=days(i);
        x=numel(num2str(d));
        if x==1;
            day=(['00',num2str(d)]);
        elseif x==2;
            day=(['0',num2str(d)]);
        elseif x==3;
            day=num2str(d);
        end
        file=(['MYD06_L2.A',year,day,'*.hdf']);
        files=dir(fullfile('\\128.32.208.19\Elements (G)\MODIS_Cloud_US',year,file));
        n=length(files);
        for ii=1:n;
            directory=['\\128.32.208.19\Elements (G)\MODIS_Cloud_US\',year];
            cd(directory) 
            filename=files(ii).name;
            fileinfo=hdfinfo(filename);

            Latitude=hdfread(fileinfo.Vgroup(1).Vgroup(1).SDS(1)); Latitude=double(Latitude); Latitude=Latitude(:);
            Longitude=hdfread(fileinfo.Vgroup(1).Vgroup(1).SDS(2)); Longitude=double(Longitude); Longitude=Longitude(:);
            CloudFraction=hdfread(fileinfo.Vgroup(1).Vgroup(2).SDS(18)); CloudFraction=double(CloudFraction); CloudFraction=CloudFraction(:); 
            CloudFraction(CloudFraction==127)=100; CloudFraction=CloudFraction*0.009999999776482582;
            ViewingZenithAngle=hdfread(fileinfo.Vgroup(1).Vgroup(2).SDS(4)); ViewingZenithAngle=double(ViewingZenithAngle); ViewingZenithAngle=ViewingZenithAngle(:); 
            ViewingZenithAngle=ViewingZenithAngle*0.009999999776482582;
            %CloudFraction(ViewingZenithAngle>40)=1;
            
            x=find(Longitude>lonmax | Longitude<lonmin);
            y=find(Latitude>latmax | Latitude<latmin);
            Longitude(x)=NaN;           Longitude(y)=NaN;               Longitude(isnan(Longitude))=[]; 
            Latitude(x)=NaN;            Latitude(y)=NaN;                Latitude(isnan(Latitude))=[];
            CloudFraction(x)=NaN;       CloudFraction(y)=NaN;           CloudFraction(isnan(CloudFraction))=[];
            ViewingZenithAngle(x)=NaN;  ViewingZenithAngle(y)=NaN;      ViewingZenithAngle(isnan(ViewingZenithAngle))=[];
    
            if isempty(Longitude)||isempty(Latitude);
            else
                if exist('Data','var')==0;
                    Data(1).Longitude=Longitude;              Data(1).Latitude=Latitude;    
                    Data(1).CloudFraction=CloudFraction;      Data(1).ViewingZenithAngle=ViewingZenithAngle;  
                elseif exist('Data','var')==1;
                    Data(1).Longitude=[Data(1).Longitude;Longitude];
                    Data(1).Latitude=[Data(1).Latitude;Latitude];
                    Data(1).CloudFraction=[Data(1).CloudFraction;CloudFraction];
                    Data(1).ViewingZenithAngle=[Data(1).ViewingZenithAngle;ViewingZenithAngle];
                end
            end
        end
        if exist('Data','var')==0
        else
            XI=24.025:0.05:50.975;
            XI=repmat(XI',1,1240);
            YI=-125.975:0.05:-64.025;
            YI=repmat(YI,540,1);
            ZI = griddata(Data(1).Longitude,Data(1).Latitude,Data(1).CloudFraction,YI,XI);
            Data(1).CloudFractionGrid=ZI;
            R = addtodate(datenum(['01/01/',year]), str2num(day)-1, 'day');
            date=datestr(R,26);
        toc
        cd('\\128.32.208.19\Elements (G)\MODIS_Cloud_US\MODIS_Cloud_mats')
        savename=(['MODIS_Cld_',date(1:4),date(6:7),date(9:10)]);
        save(savename, 'Data')
        clear Data data
        end
    end
%end
