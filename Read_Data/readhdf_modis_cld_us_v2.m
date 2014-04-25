%%readhdf_modis_cld_us
%%arr 04/22/2011
%15 Apr 2014: Updated to current SDS numbers and directory structure

%addpath C:\Ashley\NewCaRetrieval2\Regrid_tools
%addpath C:\Ashley\NewCaRetrieval2\Regrid_tools\Compute_Corner_Pts_ms

warning off all
tic

%****************************%
lonmin = -126;  lonmax = -64;
latmin = 24;    latmax = 51;
%****************************%

satellite='MODIS';
retrieval='cld_v51';

%for jj=1:length(years)   
    year='2007';%years{jj};
    days=32:62;   
    for i=1:length(days);
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
        files=dir(fullfile('/Volumes/share/GROUP/SAT/MODIS/MYD06_L2_Collection5',year,file));
        n=length(files);
        for ii=1:n;
            directory=['/Volumes/share/GROUP/SAT/MODIS/MYD06_L2_Collection5/',year];
            cd(directory) 
            filename=files(ii).name;
            fileinfo=hdfinfo(filename);
            
            fprintf('Reading file %u : %s \n',ii,filename);
            
            Latitude=hdfread(hdf_dsetID(fileinfo,1,1,'Latitude')); Latitude=double(Latitude); Latitude=Latitude(:);
            Longitude=hdfread(hdf_dsetID(fileinfo,1,1,'Longitude')); Longitude=double(Longitude); Longitude=Longitude(:);
            CloudFraction=hdfread(hdf_dsetID(fileinfo,1,2,'Cloud_Fraction')); CloudFraction=double(CloudFraction); CloudFraction=CloudFraction(:); 
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
                fprintf('    Falls within lat/lon boundaries\n');
                if exist('mod_Data','var')==0;
                    mod_Data(1).Longitude=Longitude;              mod_Data(1).Latitude=Latitude;    
                    mod_Data(1).CloudFraction=CloudFraction;      mod_Data(1).ViewingZenithAngle=ViewingZenithAngle;  
                elseif exist('mod_Data','var')==1;
                    mod_Data(1).Longitude=[mod_Data(1).Longitude;Longitude];
                    mod_Data(1).Latitude=[mod_Data(1).Latitude;Latitude];
                    mod_Data(1).CloudFraction=[mod_Data(1).CloudFraction;CloudFraction];
                    mod_Data(1).ViewingZenithAngle=[mod_Data(1).ViewingZenithAngle;ViewingZenithAngle];
                end
            end
        end
        if exist('mod_Data','var')==0
        else
            XI=24.025:0.05:50.975;
            XI=repmat(XI',1,1240);
            YI=-125.975:0.05:-64.025;
            YI=repmat(YI,540,1);
            ZI = griddata(mod_Data(1).Longitude,mod_Data(1).Latitude,mod_Data(1).CloudFraction,YI,XI);
            mod_Data(1).CloudFractionGrid=ZI;
            R = addtodate(datenum(['01/01/',year]), str2num(day)-1, 'day');
            date=datestr(R,26);
        toc
        cd('/Volumes/share/GROUP/SAT/BEHR/MODIS_cld_files_v50')
        savename=(['MODIS_Cld_',date(1:4),date(6:7),date(9:10)]);
        save(savename, 'mod_Data')
        fprintf('\t Saving %s\n',savename);
        clear mod_Data
        end
    end
%end
