%BEHR_main
%Josh Laughner <joshlaugh5@gmail.com>
%Based on BEHR_nwus by Ashley Russell (02/09/2012)
%Takes "OMI_SP_yyyymmdd.m" files produced by read_omno2_v_aug2012.m as it's
%main input.

%****************************%
% CONSOLE OUTPUT LEVEL - 0 = none, 1 = minimal, 2 = all messages, 3 = times %
% Allows for quick control over the amount of output to the console.
% Choose a higher level to keep track of what the script is doing.
DEBUG_LEVEL = 2;
%****************************%

%****************************%
%**** FILE LOCATIONS ********%
%****************************%
%This is the directory where the final .mat file will be saved. This will
%need to be changed to match your machine and the files' location. Do not
%include a trailing separator, i.e. '/my/favorite/directory' not
%'my/favorite/directory/
mat_dir = '/Volumes/share-sat/SAT/BEHR/Test_BEHR_Files';

%This is the directory where the "OMI_SP_*.mat" files are saved. This will
%need to be changed to match your machine and the files' location. Do not
%include a trailing separator, i.e. '/my/favorite/directory' not
%'my/favorite/directory/
omi_sp_dir = '/Volumes/share-sat/SAT/BEHR/SP_Files_2014';

%Add the path to the AMF_tools folder which contains rNmcTmp2.m,
%omiAmfAK2.m, integPr2.m and others.  In the Git repository for BEHR, this
%is the 'AMF_tools' folder.
amf_tools_path = '/Users/Josh/Documents/MATLAB/BEHR/AMF_tools';
addpath(amf_tools_path)

%This is the directory where the NO2 profiles are stored. This will
%need to be changed to match your machine and the files' location. Do not
%include a trailing separator, i.e. '/my/favorite/directory' not
%'my/favorite/directory/
no2_profile_path = '/Volumes/share/GROUP/SAT/BEHR/Monthly_NO2_Profiles';

%Store paths to relevant files
fileTmp = fullfile(amf_tools_path,'nmcTmpYr.txt');
fileDamf = fullfile(amf_tools_path,'damf.txt');
fileNO2 = fullfile(amf_tools_path,'PRFTAV.txt');
%****************************%

%****************************%
%Process all files between these dates, in yyyy/mm/dd format
%****************************%
date_start='2011/07/01';
date_end='2011/07/03';
%****************************%

%These will be included in the file name
%****************************%
satellite='OMI';
retrieval='BEHR_omiCloudAMF';
%****************************%

%****************************%
% Which cloud product to use to calculate the AMF: OMI or MODIS
%****************************%
cloud_amf = 'omi';
cloud_rad_amf = 'omi';
%****************************%

%Add the path to the AMF_tools folder which contains rNmcTmp2.m,
%omiAmfAK2.m, integPr2.m and others.  In the Git repository for BEHR, this
%is the 'AMF_tools' folder.
addpath('/Users/Josh/Documents/MATLAB/BEHR/AMF_tools')

%Find the last file completed and set the start date to the next day.  This
%will allow the process to be stopped and started with minimum
%intervention.
file_prefix = [satellite,'_',retrieval,'_']; l = length(file_prefix);
last_file=dir(fullfile(mat_dir,sprintf('%s*.mat',file_prefix)));

if ~isempty(last_file)
    last_datenum = datenum(last_file(end).name(l+1:l+8),'yyyymmdd')+1;
else
    last_datenum = 0;
end

if last_datenum >= datenum(date_start) && last_datenum <= datenum(date_end)
    datenums = last_datenum:datenum(date_end);
else
    datenums = datenum(date_start):datenum(date_end);
end

tic
for j=1:length(datenums)
    %Read the desired year, month, and day
  	R=datenums(j);
    date=datestr(R,26);
    year=date(1:4);
    month=date(6:7);
    day=date(9:10);
    if DEBUG_LEVEL > 0; disp(['Processing data for ', date]); end
    
    filename = ['OMI_SP_',year,month,day,'.mat'];
    if DEBUG_LEVEL > 1; disp(['Looking for SP file ',fullfile(omi_sp_dir,filename),'...']); end
    if isequal(exist(fullfile(omi_sp_dir,filename),'file'),0)
        if DEBUG_LEVEL > 0; disp('No SP file exists for given day'); end
        continue
    else
        if DEBUG_LEVEL > 1; fprintf('\t ...Found.\n'); end
        load(fullfile(omi_sp_dir,filename)); %JLL 17 Mar 2014: Will load the variable 'Data' into the workspace
        if exist('profile_file','file')==1 && strcmp(profile_file(2:3),month)==1; %JLL 20 Mar 2014: 
        else
            profile_file=['m',month,'_NO2_profile'];
            if DEBUG_LEVEL > 1; disp(['Loading ',fullfile(no2_profile_path,profile_file)]); end
            load(fullfile(no2_profile_path,profile_file));
        end
        for d=1:length(Data);
            if Data(d).Longitude==0 | length(Data(d).Longitude)==1;
                if DEBUG_LEVEL > 1; fprintf('  Note: Data(%u) is empty\n',d); end
                continue %JLL 17 Mar 2014: Skip doing anything if there's really no information in this data
            else
                if DEBUG_LEVEL>0; fprintf('  Swath %u of %s \n',d,date); end
                c=numel(Data(d).Longitude);
                
                Data(d).MODISAlbedo(isnan(Data(d).MODISAlbedo)==1)=0; %JLL 17 Mar 2014: replace NaNs with fill values
                Data(d).GLOBETerpres(isnan(Data(d).GLOBETerpres)==1)=400.0000;
                
                %JLL 17 Mar 2014: Load some of the variables from 'Data' to
                %make referencing them less cumbersome.
                lon = Data(d).Longitude;
                lat = Data(d).Latitude;
                loncorns=Data(d).Loncorn;
                latcorns=Data(d).Latcorn;
                
                sza = Data(d).SolarZenithAngle;
                vza = Data(d).ViewingZenithAngle;
                phi = Data(d).RelativeAzimuthAngle;
                
                mon = str2double(month)*ones(size(Data(d).Latitude));
                pressure = [1020 1015 1010 1005 1000 990 980 970 960 945 925 900 875 850 825 800 770 740 700 660 610 560 500 450 400 350 280 200];% 100 50 20 5];
                if DEBUG_LEVEL > 1; disp('   Interpolating temperature data'); end
                [temperature, tmpSAVE] = rNmcTmp2(fileTmp, pressure, lon, lat, mon); %JLL 17 Mar 2014: Interpolates temperature values to the pressures and lat/lon coordinates desired
                
                surfPres = Data(d).GLOBETerpres;
                albedo = Data(d).MODISAlbedo;
                
                surfPres(surfPres>=1013)=1013; %JLL 17 Mar 2014: Clamp surface pressure to sea level or less.
                cldPres = Data(d).CloudPressure;
                
                if DEBUG_LEVEL > 1; disp('   Calculating clear and cloudy AMFs'); end
                dAmfClr = rDamf2(fileDamf, pressure, sza, vza, phi, albedo, surfPres); %JLL 18 Mar 2014: Interpolate the values in dAmf to the albedo and other conditions input
                cloudalbedo=0.8*ones(size(Data(d).CloudFraction)); %JLL 18 Mar 2014: Assume that any cloud has an albedo of 0.8
                dAmfCld = rDamf2(fileDamf, pressure, sza, vza, phi, cloudalbedo, cldPres); %JLL 18 Mar 2014: Interpolate dAmf again, this time taking the cloud top and albedo as the bottom pressure
                
                if DEBUG_LEVEL > 1; disp('   Reading NO2 profiles'); end
                [no2_bins] = rProfile_US(PROFILE, loncorns, latcorns, c); %JLL 18 Mar 2014: Bins the NO2 profiles to the OMI pixels; the profiles are averaged over the pixel
                no2_bins = reshape(no2_bins,length(pressure),size(vza,1),size(vza,2));
                no2Profile1 = no2_bins ./ (10^6); % NO2 from WRF is in ppm in these files
                prof_i = zeros(size(Data(d).Latitude)); prof_i(isnan(squeeze(no2Profile1(1,:,:)))==1)=1; %JLL 18 Mar 2014: prof_i is a matrix of 0 or 1s that is 1 wherever the bottom of NO2 profile is NaN
                clear no2_bins no2_mixratio
                
                no2Profile2 = no2Profile1;
                pTerr = surfPres;
                pCld = cldPres;
                if strcmpi(cloud_amf,'omi')
                    cldFrac = Data(d).CloudFraction; %JLL 18 Mar 2014: Cloud fraction and radiance fraction are scaled by 1e-3 in the OMNO2 he5 file
                else
                    cldFrac = Data(d).MODISCloud;
                end
                if strcmpi(cloud_rad_amf,'omi')
                    cldRadFrac = Data(d).CloudRadianceFraction;
                else
                    cldRadFrac = Data(d).MODISCloud;
                end
                
                if DEBUG_LEVEL > 1; disp('   Calculating BEHR AMF'); end
                noGhost=1; ak=1;
                [amf, amfCld, amfClr] = omiAmfAK2(pTerr, pCld, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile1, no2Profile2, noGhost, ak); %JLl 18 Mar 2014: The meat and potatoes of BEHR, where the TOMRAD AMF is adjusted to use the GLOBE pressure and MODIS cloud fraction
                amf(prof_i==1)=NaN;
                
                Data(d).BEHRAMFTrop = amf; %JLL 18 Mar 2014: Save the resulting AMF of the pixel
            end
        end
        
        b=length(Data);
        for z=1:b;
            if isfield(Data,'BEHRAMFTrop')==0 || isempty(Data(z).BEHRAMFTrop)==1;
                continue
            else
                Data(z).BEHRColumnAmountNO2Trop=Data(z).ColumnAmountNO2Trop.*Data(z).AMFTrop./Data(z).BEHRAMFTrop;
                if DEBUG_LEVEL > 0; fprintf('   BEHR [NO2] stored for swath %u\n',z); end
            end
        end
        
        
        addpath('/Users/Josh/Documents/MATLAB/BEHR/Utils/m_map'); %JLL 18 Mar 2014: Adds the path to the m_map toolkit, needed for hdf_quadrangle_5km_new
        
        %*********************************%
        %JLL 19 Mar 2014: These will be used to define the quadrangles -
        %the quads will be smaller than the OMI pixel, and multiple quads
        %will take on the value for the same (closest) OMI pixel.  By
        %keeping the quads' centers the same over all retrievals you wish
        %to average, this will allow easier averaging over multiple OMI
        %swaths. This is a form of oversampling.
        %*********************************%
        lonmin = -125;  lonmax = -65;  
        latmin = 25;   latmax = 50;
        resolution = 0.05; resolution2 = 0.05;
        %*********************************%
        %
        if lonmin > lonmax %Just in case I enter something backwards...
            error('behr_main:lon_maxmin','Lonmin is greater than lonmax')
        elseif latmin > latmax
            error('behr_main:lat_maxmin', 'Latmin is greater than latmax')
        end
        
        %*********************************%
        %JLL 19 Mar 2014: Save all relevant values produced by add2grid to
        %a new structure called 'OMI'
        %*********************************%
        
        if DEBUG_LEVEL > 0; disp('  Preparing OMI structure'); end
        s=size(Data);
        hh=0;
        for d=1:s(2);
            if Data(d).ViewingZenithAngle==0;
            elseif numel(Data(d).ViewingZenithAngle)==1;
                continue
            else
                if DEBUG_LEVEL > 1; fprintf('   Gridding data for swath %u\n',d); end
                hh=hh+1;
                OMI(hh) = add2grid_BEHR(Data(d),resolution,resolution2,[lonmin, lonmax],[latmin, latmax]); %JLL 20 Mar 2014: Superimpose data to a grid determined by lat & lon min/max and resolution above. Default resolution is 0.05 degree
            end
        end

        % Assign these after the rest of the OMI structure to avoid a
        % problem with dissimilar structure assignment
        Latitude=(latmin+0.025):resolution:(latmax-0.025); Latitude=Latitude'; Latitude=repmat(Latitude,1,1200);
        Longitude=(lonmin+0.025):resolution2:(lonmax-0.025); Longitude=repmat(Longitude,500,1);
        for aa=1:hh
            OMI(aa).Latitude = Latitude;
            OMI(aa).Longitude = Longitude;
        end
        savename=[file_prefix,year,month,day];  
        if DEBUG_LEVEL > 0; disp(['   Saving data as',fullfile(mat_dir,savename)]); end
        save(fullfile(mat_dir,savename),'Data','OMI')
        toc
        t=toc;
        %if t>1200
            %error('Time exceeded 20 min. Stopping')
        %end
    end
end
