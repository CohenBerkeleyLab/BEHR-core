function BEHR_InSitu_Reprocessing
%BEHR_InSitu_Reprocessing
%
%   This script will take aircraft data and use it to recalculate AMFs and
%   produce a new satellite column using that AMF.
%
%   Returns a quality flag with each bit representing a specific warning
%   about the data.  These mimic the flags in the spiral verification code.
%       1st: Summary, set to 1 if any flags are set.
%       2nd: Reserved as a second summary bit against future need.
%       3rd: Unused
%       4th: Indicates that < 10% of the data points in the profile had NO2
%           data
%       5-15: Unused
%       16th: Set if the column was skipped due to < 1% valid
%           NO2, pressure, or temperature data
%
%   Josh Laughner <joshlaugh5@gmail.com> 18 Aug 2013

E = JLLErrors;

campaign_name = 'discover-tx';

[Names, merge_dates, merge_dir] = merge_field_names(campaign_name);

start_date = merge_dates{1};
end_date = merge_dates{2};

starttime = '12:00';
endtime = '15:00';

%Which clouds to use for the AMF calculation; 'CloudFraction' for OMI and
%'MODISCloud' for MODIS
cld_field = 'CloudFraction';

%Which NO2 field from the aircraft file to use; for MD options are
%'NO2_NCAR' and 'NO2_LIF'; for CA and TX, 'NO2_MixingRatio' or
%'NO2_MixingRatio_LIF'
no2field = Names.no2_lif;

%The directory where the original BEHR files are located
behr_prefix = 'OMI_BEHR_omi*';
behr_dir = '/Volumes/share-sat/SAT/BEHR/DISCOVER_BEHR/';


%The file prefix and directory to save the resulting files under
save_prefix = 'OMI_BEHR_InSitu_';
save_dir = '/Volumes/share-sat/SAT/BEHR/DISCOVER_BEHR_REPROCESSED/';

amf_tools_path = '/Users/Josh/Documents/MATLAB/BEHR/AMF_tools';
fileTmp = fullfile(amf_tools_path,'nmcTmpYr.txt');
fileDamf = fullfile(amf_tools_path,'damf.txt');

DEBUG_LEVEL = 2;

%%%%% END USER INPUT %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BEGIN FUNCTION %%%%%

fill_val = -9e9; % defined fill value

% Make sure the save prefix ends in an underscore
if ~strcmp(save_prefix(end),'_'); save_prefix(end+1) = '_'; end

dates = datenum(start_date):datenum(end_date);
for d=1:numel(dates)
    curr_date = datestr(dates(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    
    if DEBUG_LEVEL > 0; fprintf('Now on %s\n',curr_date); end
    if DEBUG_LEVEL > 1; fprintf('  Loading data...\n'); end
    
    % Load the merge and BEHR file
    merge_filename = sprintf('*%s_%s_%s.mat',year,month,day);
    behr_filename = sprintf('%s%s%s%s.mat',behr_prefix,year,month,day);
    
    merge_files = dir(fullfile(merge_dir,merge_filename));
    if numel(merge_files)==1
        load(fullfile(merge_dir, merge_files(1).name),'Merge')
    elseif isempty(merge_files)
        if DEBUG_LEVEL > 1; fprintf('No Merge file for %s\n',datestr(dates(d))); end
        continue
    else
        error('run_spiral:tmm','Number of merge files for %s is not 1 or 0',datestr(dates(d)));
    end
    
    behr_files = dir(fullfile(behr_dir,behr_filename));
    if numel(behr_files)==1
        load(fullfile(behr_dir,behr_files(1).name),'Data')
    elseif isempty(behr_files)
        if DEBUG_LEVEL > 1; fprintf('No BEHR file for %s\n',datestr(dates(d))); end
        continue
    else
        error('run_spiral:tmm','Number of BEHR files for %s is not 1 or 0',datestr(dates(d)));
    end
    
    flag_base = uint16(0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% AIRCRAFT DATA PREPROCESSING %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Load the aircraft data and check if greater than 99% of it is
    %nans; if so, fill the new BEHR fields with NaNs and skip this
    %loop.  If >90% of the NO2 data are nans, then set a warning flag
    %but continue processing
    [no2, utc, pres, lon, lat] = remove_merge_fills(Merge,no2field,'alt','PRESSURE');
    no2(no2<0) = NaN; % Because of the exponential shape of the NO2 profile, negative values lead to imaginary values for the AMFs
    
    percent_no2_nans = sum(isnan(no2))/numel(no2); percent_pres_nans = sum(isnan(pres))/numel(pres);
    if percent_no2_nans > 0.99 || percent_pres_nans > 0.99;
        flag_base = bitset(flag_base,16,1);
        if DEBUG_LEVEL > 1; fprintf('NO2 %.2f NaNs, PRESSURE %.2f NaNs, skipping\n', percent_no2_nans, percent_pres_nans); end
        Data = return_null(Data,flag_base);
        saveData(Data,curr_date,save_dir,save_prefix);
        continue
    elseif percent_no2_nans > 0.9 || percent_pres_nans > 0.9
        flag_base = bitset(flag_base,4,1);
        if DEBUG_LEVEL > 1; fprintf('NO2 %.2f NaNs, PRESSURE %.2f NaNs, setting warning flag\n', percent_no2_nans, percent_pres_nans); end
    end
    
    % Load the profile numbers 
    profnum = remove_merge_fills(Merge, Names.profile_numbers);
    
    % Calculate the UTC offset to match flight times to satellite overpass
    tz = round(nanmean(lon)/15);
    
    % Get all unique profile numbers and their start times
    unique_profnums = unique(profnum(profnum~=0));
    start_times = zeros(numel(unique_profnums),1);
    for a=1:numel(unique_profnums)
        xx = profnum == unique_profnums(a);
        start_times(a) = min(utc(xx));
    end
    
    % Remove from consideration any profiles with a start time before 10:45
    % am or after 4:45 pm local standard time
    yy = start_times >= local2utc(starttime,tz) & start_times <= local2utc(endtime,tz);
    unique_profnums = unique_profnums(yy); start_times = start_times(yy);
    
    % Save each profile's NO2, altitude, radar altitude, latitude, and
    % longitude as an entry in a cell array
    s = size(unique_profnums);
    no2_array = cell(s); utc_array = cell(s);
    lat_array = cell(s); lon_array = cell(s); 
    pres_array = cell(s); profnum_array = cell(s);
    for a=1:numel(unique_profnums)
        xx = profnum == unique_profnums(a);
        no2_array{a} = no2(xx);
        lat_array{a} = lat(xx);
        lon_array{a} = lon(xx);
        pres_array{a} = pres(xx);
        utc_array{a} = utc(xx);
        profnum_array{a} = unique_profnums(a);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% RECALCULATE VCDs %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for s=1:numel(Data)
        if DEBUG_LEVEL > 1; fprintf('\tNow on swath %d of %d\n',s,numel(Data)); end
        %Initialize the new fields in the Data structure to hold the AMF
        %calculated with the in-situ column and the reprocesses satellite
        %column itself.
        Data(s).InSituAMF = fill_val * ones(size(Data(s).Latitude));
        Data(s).BEHR_R_ColumnAmountNO2Trop = fill_val * ones(size(Data(s).Latitude));
        Data(s).ProfileCount = fill_val * ones(size(Data(s).Latitude));
        Data(s).InSituFlags = fill_val * ones(size(Data(s).Latitude));
        
        flags = uint16(double(flag_base)*ones(size(Data(s).Latitude)));
        
        %Now iterate through each profile.  Find the pixels it intersects
        %and use the profile to calculate a new AMF for those pixels        
        omi_lat = Data(s).Latitude; omi_lon = Data(s).Longitude;
        corner_lat = Data(s).Latcorn; corner_lon = Data(s).Loncorn;
        omi_no2 = Data(s).ColumnAmountNO2Trop;
        sza = Data(s).SolarZenithAngle; vza = Data(s).ViewingZenithAngle;
        phi = Data(s).RelativeAzimuthAngle;  albedo = Data(s).MODISAlbedo;
        surfPres = Data(s).GLOBETerpres; cloudPres = Data(s).CloudPressure;
        cldFrac = Data(s).(cld_field); cldRadFrac = Data(s).CloudRadianceFraction;
        behr_no2 = Data(s).BEHRColumnAmountNO2Trop; 
        
        
        % This will be where the amfs due to each profile go to be averaged
        % together to find the average amf for the pixel
        tmp_amfs = nan(numel(behr_no2),numel(no2_array));
        tmp_count = zeros(size(behr_no2));
        for p=1:numel(no2_array)
            
            % Find all the pixels that impinge on this profile
            lontest = lon_array{p}; lontest(isnan(lontest)) = [];
            lattest = lat_array{p}; lattest(isnan(lattest)) = [];
            lat_logical = max(lattest) > min(corner_lat,[],1) & min(lattest) < max(corner_lat,[],1);
            lon_logical = max(lontest) > min(corner_lon,[],1) & min(lontest) < max(corner_lon,[],1);
            latlon_logical = lat_logical & lon_logical;
            
            pix_indices = find(latlon_logical);
            
            % Skip this profile if no pixels intersect it
            if sum(latlon_logical) == 0; 
                if DEBUG_LEVEL > 1; fprintf('\t\tNo pixels fall within initial lat/lon boundaries\n'); end
                continue; 
            end
            
            loncorn_p = corner_lon(:,latlon_logical); latcorn_p = corner_lat(:, latlon_logical);
            
            % Finally actually check if the profile falls in the pixel
            % using inpolygon(). Recall that we require there to be 20
            % valid measurements in the lowest 3 km (~ 675 hPa).
            no2_3km = no2_array{p}(pres_array{p}>675);
            lon_3km = lon_array{p}(pres_array{p}>675);
            lat_3km = lat_array{p}(pres_array{p}>675);
            
            pix_coverage = zeros(size(loncorn_p,2),1);
            pix_xx = true(size(loncorn_p,2),1);
            
            for pix=1:size(loncorn_p,2)
                IN_3km = inpolygon(lon_3km, lat_3km, loncorn_p(:,pix), latcorn_p(:,pix));
                if sum(~isnan(no2_3km(IN_3km)))<20
                    pix_xx(pix) = false;
                    continue % Pixel must have 20 valid measurments between 0-3 km altitude (good sampling of boundary layer)
                end
                % Calculate what percentage of the profile actually falls in
                % this pixel; append to all values for this pixel
                IN_all = inpolygon(lon_array{p}, lat_array{p}, loncorn_p(:,pix), latcorn_p(:,pix));
                pix_coverage(pix) = sum(IN_all)/numel(no2_array{p});
            end
            
            %If there are no pixels that have part of the profile below 3
            %km within them, skip this profile
            if sum(pix_indices)==0; 
                if DEBUG_LEVEL > 1; fprintf('\t\tNo pixels fall within strict lat/lon boundaries\n'); end
                continue; 
            end
            
            pix_indices = pix_indices(pix_xx);
                        
            omi_lon_p = omi_lon(pix_indices); omi_lat_p = omi_lat(pix_indices);
            loncorn_p = corner_lon(:,pix_indices); latcorn_p = corner_lat(:, pix_indices);
            omi_no2_p = omi_no2(pix_indices); behr_no2_p = behr_no2(pix_indices);
            sza_p = sza(pix_indices); vza_p = vza(pix_indices);
            phi_p = phi(pix_indices); albedo_p = albedo(pix_indices);
            surfPres_p = surfPres(pix_indices); cloudPres_p = cloudPres(pix_indices);
            cldFrac_p = cldFrac(pix_indices); cldRadFrac_p = cldRadFrac(pix_indices);
            amfs_p = nan(size(omi_lon_p));
            
            if DEBUG_LEVEL > 1; fprintf('\t\tRecalculating AMFs\n'); end
            for pix=1:numel(omi_lon_p)
                % Extrapolate and bin the profile using median extrapolation on
                % the bottom and by inserting a scaled WRF profile on the top,
                % if needed.  Only the surface pressure changes for each
                % pixel.
                [insitu_profile, insitu_pressures] = extrapolate_profile(no2_array{p}, pres_array{p}, 'surfacePressure', surfPres_p(pix),'top', 'wrf-scaled', 'bottom', 'median',...
                    'utc',nanmean(utc_array{p}),'sitenum',profnum_array{p}, 'month',month,'lat',lat_array{p},'lon',lon_array{p},'date',curr_date,'shape','exp');
                
                % The profiles and pressure must be columns for the AMF calculation to
                % function correctly
                if ~iscolumn(insitu_pressures); insitu_pressures = insitu_pressures'; end
                if ~iscolumn(insitu_profile); insitu_profile = insitu_profile'; end
                
                % A surface pressure greater than 1013 will cause the dAmf
                % interpolation to fail.
                this_surfPres = min([surfPres_p(pix),1013]);
                
                % Calculate the new AMF
                mean_lon = nanmean(lon_array{p});
                mean_lat = nanmean(lat_array{p});
                if isnan(mean_lon) || isnan(mean(lat))
                    E.badvar('mean_lon/mean_lat','the mean lon/lat is a NaN');
                end
                [temperature, ~] = rNmcTmp2(fileTmp, insitu_pressures, nanmean(lon_array{p}), nanmean(lat_array{p}), str2double(month));
                dAmfClr2 = rDamf2(fileDamf, insitu_pressures, sza_p(pix), vza_p(pix), phi_p(pix), albedo_p(pix), this_surfPres);
                cloudalbedo=0.8;
                dAmfCld2 = rDamf2(fileDamf, insitu_pressures, sza_p(pix), vza_p(pix), phi_p(pix), cloudalbedo, cloudPres_p(pix));
                
                noGhost = 1; ak = 0;
                amfs_p(pix) = omiAmfAK2(this_surfPres, cloudPres_p(pix), cldFrac_p(pix), cldRadFrac_p(pix), insitu_pressures, dAmfClr2, dAmfCld2, temperature, insitu_profile, insitu_profile, noGhost, ak);
            end
            
            tmp_amfs(pix_indices,p) = amfs_p;
            tmp_count(pix_indices) = tmp_count(pix_indices) + 1;
        end
        
        new_amfs = nanmean(tmp_amfs,2);
        new_amfs = reshape(new_amfs,size(omi_lat));
        new_columns = Data(s).BEHRColumnAmountNO2Trop .* Data(s).BEHRAMFTrop ./ new_amfs;
        
        Data(s).InSituAMF = new_amfs;
        Data(s).BEHR_R_ColumnAmountNO2Trop = new_columns;
        Data(s).ProfileCount = tmp_count;
        Data(s).InSituFlags = flags;
    end
        
    
    
    %Run these checks in a separate for loop so that any continue
    %statements in the first don't cause them to be skipped.
    for s=1:numel(Data)
        % Check if either of the new fields is still a fill value; if so
        % throw a warning to let the user know that for some reason the
        % values were not set at some point in the loop
        if any(Data(s).InSituAMF == fill_val);
            warning('In-situ AMF not set for swath %d on %s',s,curr_date);
        end
        if any(Data(s).BEHR_R_ColumnAmountNO2Trop == fill_val);
            warning('Reprocessed column not set for swath %d on %s',s,curr_date);
        end
        if any(Data(s).InSituFlags == fill_val);
            warning('Flags not set for swath %d on %s',s,curr_date);
        end
    end
    saveData(Data,curr_date,save_dir,save_prefix);
    clear Data
end

function data = return_null(data,flag)
for i=1:numel(data)
    data(i).InSituAMF = nan(size(data(i).Longitude));
    data(i).BEHR_R_ColumnAmountNO2Trop = nan(size(data(i).Longitude));
    data(i).ProfileCount = zeros(size(data(i).Longitude));
    data(i).InSituFlags = uint16(double(flag) * ones(size(data(i).Longitude)));
end

function saveData(Data,curr_date,save_dir,prefix)
year = curr_date(1:4);
month = curr_date(6:7);
day = curr_date(9:10);
savename = sprintf('%s%s%s%s.mat',prefix,year,month,day);
save(fullfile(save_dir,savename),'Data');

