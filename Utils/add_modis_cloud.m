function [ Data ] = add_modis_cloud( Data )
%add_modis_cloud(Data) Adds MODIS cloud information to OMI_SP "Data" structures
%   Because MODIS cloud products are so large, I have been running
%   read_omno2_v_aug2012 without the cloud files available to average to
%   pixels.  This function will just average the MODIS cloud fraction to
%   the OMI pixels and save them in the data structure that it returns
%
%   If MODIS cloud files are present in the proper directory when
%   read_omno2 is run, this function is not necessary.
DEBUG_LEVEL = 2;
he5_dir = '/Volumes/share/GROUP/SAT/OMI/OMNO2_32';
modis_myd06_dir = '/Volumes/share-sat/SAT/MODIS/MYD06_L2';

for d=1:numel(Data)
    curr_date = datestr(Data(d).Date,29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Find the swaths %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    short_filename=['OMI-Aura_L2-OMNO2_',year,'m',month,day,'*.he5'];
    file_dir = fullfile(he5_dir,year,month); %Used both here to find all he5 files and in the swath for loop to identify each file.
    file=fullfile(file_dir,short_filename);
    sp_files = dir(file);
    n = length(sp_files);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% Add the cloud data %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Convert the OMI date to a Julian calendar day
    d2=1+datenum(str2double(year),str2double(month),str2double(day))-datenum(str2double(year),1,1);
    x=numel(num2str(d2));
    if x==1;
        julian_day=(['00',num2str(d2)]);
    elseif x==2;
        julian_day=(['0',num2str(d2)]);
    elseif x==3;
        julian_day=num2str(d2);
    end
    
    %Find all MODIS files that occur after the current OMI file
    %but before the next OMI file.
    modis_file=(['MYD06_L2.A',year,julian_day,'*.hdf']);
    modis_files=dir(fullfile(modis_myd06_dir,year,modis_file));
    
    % Calculate the start time for the next OMI swath.
    % Usually there is at least one swath starting after
    % the one overflying the US west coast for that day,
    % but occasionally that swath is not present in the
    % OMNO2 data.  In those cases, we need to calculate the
    % time it should have started, knowing that the Aura
    % orbit period is ~99 min.  If for some reason the
    % calculated start time should end up being in the next
    % day, error out so the user is aware of that.
    if e < n
        next_omi_swath_time = str2double(sp_files(e+1).name(29:32));
    else
        omi_hr = str2double(sp_files(e).name(29:30));
        omi_min = str2double(sp_files(e).name(31:32));
        next_omi_swath_time = 100*(floor((omi_min + 99)/60)+omi_hr) + mod(omi_min + 99,60);
        if next_omi_swath_time > 2359; error('read_omno2:modis_cloud','Next OMI swath time for MODIS cloud binning calculated to be in the next day. \nManual attention recommended'); end
    end
    
    for ii=1:length(modis_files);
        mod_filename=modis_files(ii).name;
        if str2double(mod_filename(19:22))<str2double(sp_files(e).name(29:32));
            continue
        elseif str2double(mod_filename(19:22))>next_omi_swath_time;
            continue
        else
            %For each file that fits the criteria mentioned
            %above, import its latitude, longitude, and cloud
            %fraction.
            if DEBUG_LEVEL > 0; fprintf('  Averaging MODIS cloud file %s\n',mod_filename); end
            mod_filename=fullfile(modis_myd06_dir,year,modis_files(ii).name); %Redefine the filename to have the full path to the file
            mod_fileinfo=hdfinfo(mod_filename);
            Latitude=hdfread(hdf_dsetID(mod_fileinfo,1,1,'Latitude')); Latitude=double(Latitude); Latitude=Latitude(:);
            Longitude=hdfread(hdf_dsetID(mod_fileinfo,1,1,'Longitude')); Longitude=double(Longitude); Longitude=Longitude(:);
            CloudFraction=hdfread(hdf_dsetID(mod_fileinfo,1,2,'Cloud_Fraction')); CloudFraction=double(CloudFraction); CloudFraction=CloudFraction(:);
            CloudFraction(CloudFraction==127)=100; CloudFraction=CloudFraction*0.009999999776482582;
            
            x=find(Longitude>lonmax | Longitude<lonmin);
            y=find(Latitude>latmax | Latitude<latmin);
            Longitude(x)=NaN;           Longitude(y)=NaN;               Longitude(isnan(Longitude))=[];
            Latitude(x)=NaN;            Latitude(y)=NaN;                Latitude(isnan(Latitude))=[];
            CloudFraction(x)=NaN;       CloudFraction(y)=NaN;           CloudFraction(isnan(CloudFraction))=[];
            
            if isempty(Longitude)||isempty(Latitude);
            else
                if exist('mod_Data','var')==0;
                    mod_Data(1).Longitude=Longitude;              mod_Data(1).Latitude=Latitude;
                    mod_Data(1).CloudFraction=CloudFraction;
                elseif exist('mod_Data','var')==1;
                    mod_Data(1).Longitude=[mod_Data(1).Longitude;Longitude];
                    mod_Data(1).Latitude=[mod_Data(1).Latitude;Latitude];
                    mod_Data(1).CloudFraction=[mod_Data(1).CloudFraction;CloudFraction];
                end
            end
        end
    end
    
    %If there is no "mod_Data" variable, fill the regular Data
    %field with -127. Otherwise, find all the MODIS cloud
    %pixels in each OMI pixel and average them together.
    if exist('mod_Data','var')==0;
        Data(d).MODISCloud=-127*ones(length(Data(d).Latitude),1);
    else
        for jj=1:length(Data(d).Latitude);
            x = [];
            x1 = Data(d).Loncorn(1,jj);   y1 = Data(d).Latcorn(1,jj);
            x2 = Data(d).Loncorn(2,jj);   y2 = Data(d).Latcorn(2,jj);
            x3 = Data(d).Loncorn(3,jj);   y3 = Data(d).Latcorn(3,jj);
            x4 = Data(d).Loncorn(4,jj);   y4 = Data(d).Latcorn(4,jj);
            
            xall=[x1;x2;x3;x4;x1];
            yall=[y1;y2;y3;y4;y1];
            xx = inpolygon(mod_Data.Latitude,mod_Data.Longitude,yall,xall);
            
            cld_vals=mod_Data.CloudFraction(xx);
            cld_vals(isnan(cld_vals))=[];
            Data(d).MODISCloud(jj,1)=mean(cld_vals);
            Data(d).MODIS_Cloud_File=mod_filename;
            
            clear xx
        end
    end
    clear mod_Data
end
end

