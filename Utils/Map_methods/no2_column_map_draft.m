%no2_column_map_2014(startdate,enddate,varargs)
%NO2 Map Function - Uses the m_map package to draw maps of NO2 column density
%over the US (primarily). Arguments:
%   startdate = a string in yyyy/mm/dd format that represents the starting
%      date of the period to average over. If you want to average multiple,
%      noncontinguous periods, enter a cell array with the start dates for
%      each period.
%   enddate = has the same format and structure as startdate, but is the
%      ending dates of the time period(s)
%   fileprefix = 'OMI_BEHR_' by default. The part of the .mat filenames
%       before the date.  
%   flags = a cell array of flags that change the behavior of the function:
%         'manual_latlon' --> allows input of manual latitude and longitude boundaries (see next 2 parameters)
%         'weekend'/'weekday' --> only average over weekend or weekdays respectively
%   lons = a 1x2 or 1x3 numerical matrix containing lonmin, lonmax, and
%       (optionally) lon resolution. If no resolution specified, set to 0.05.
%   lats = see 'lons'
%   clouds = Use OMI cloud fraction (default) or MODIS cloud fraction
%       ('modis').
%   cloudfraccrit = The maximum allowable cloud fraction. Defaults to 200
%       (0.2 unscaled value) for OMI and 0 for MODIS
%   rowanomaly = Parse mode for row anomaly (see function
%       "omi_rowanomaly"). Default value is 'XTrackFlags', other options
%       are 'AlwaysByRow', 'RowsByTime', and 'XTrackFlagsLight'.
%
%Josh Laughner 20 Mar 2014 <joshlaugh5@gmail.com> 24 Mar 2014

function [] = no2_column_map_2014(varargin)

p =inputParser;
p.addRequired('startdate',@isstr || @iscell);
p.addRequired('enddate',@isstr || @iscell);
p.addParamValue('fileprefix','OMI_BEHR_',@isstr);
p.addParamValue('flags',{},@iscell);
p.addParamValue('lons',[-125 -65], @(x) (length(x)==2 || length(x)==3) & (x(1) < x(2)));
p.addParamValue('lats',[25 50], @(x) (length(x) == 2 || length(x)==3) & (x(1) < x(2)));
p.addParamValue('clouds','omi',@isstr);
p.addParamValue('cloudfraccrit',-1,@isscalar)
p.addParamValue('rowanomaly','XTrackFlags',@(x) strcmpi(x,{'AlwaysByRow','RowsByTime','XTrackFlags','XTrackFlagsLight'})) %Ensure that the rowanomaly value is one of the allowed 4

p.parse; input = p.Results;
startdates_in = input.startdate;
enddates_in = input.enddate;
file_prefix = input.fileprefix;


%****************************%
% CONSOLE OUTPUT LEVEL - 0 = none, 1 = minimal, 2 = all messages, 3 = times %
% Allows for quick control over the amount of output to the console.
% Choose a higher level to keep track of what the script is doing.
DEBUG_LEVEL = 0;
%****************************%

%****************************%
%       FILE PATHS
%****************************%
%The directory where the BEHR output files are stored. Do not include a
%trailing separator.
behr_dir = '/Volumes/share/GROUP/SAT/BEHR/Test_BEHR_files';


%****************************%
%       FLAG PARSING         %
%****************************%
% Sets 'man_latlon' to 1 if the flag 'manual_latlon' is present. This
% should be used if the latitude and longitude boundaries are different for
% any reason.
man_latlon = any(strcmpi('manual_latlon',input.flags));

% Sets to 1 if flags for 'weekend' or 'weekday' are present
weekend_bool = any(strcmpi('weekend',input.flags));
weekday_bool = any(strcmpi('weekday',input.flags));

% Sets to 1 if 'one2one' flag is present
one2one_bool = any(strcmpi('one2one',input.flags));

%Sets lat and lon boundaries to manual input values if the 'man_latlon'
%flag is true.  If only a min and max are supplied, default resolution to
%0.05 degrees. If lat/lon boundaries are to be set automatically, that will
%be done in the main loop when the first file is loaded
if man_latlon
    latmin = input.lats(1) ; latmax = input.lats(2); 
    if length(input.lats) > 2 
        latres = input.lats(3);
    else
        latres = 0.05;
    end
    lonmin = input.lons(1); lonmax = input.lons(2); 
    if length(input.lons) > 2
        lonres = input.lons(3);
    else
        lonres = 0.05;
    end
    lats = latmin:latres:latmax;
    lons = lonmin:lonres:lonmax;
end

%****************************%
%      CLOUD FRACTION        %
%****************************%
cloud_type = input.clouds;
cloud_frac = input.cloudfraccrit;
if strcmpi(cloud_type,'omi') && cloud_frac < 0
    cloud_frac = 200; %Set the cloud fraction criteria to 200 (20% unscaled) if OMI clouds used and no other value given
elseif strcmpi(cloud_type,'modis') && cloud_frac < 0
    cloud_frac = 0;
elseif strcmpi(cloud_type, 'omi') || strcmpi(cloud_type, 'modis') %Check that the cloud type is recognized, if the value of cloud_frac is valid
else
    error('no2_col_map:cloud_type','Cloud type must be "OMI" or "MODIS"')
end

%Add the 'Utils' folder and all subfolders to MATLAB's search path. Within
%the Git repository for BEHR, this is the /Utils folder. This includes the
%m_map package.
addpath(genpath('/Users/Josh/Documents/MATLAB/BEHR/Utils'))

%****************************%
%       MAIN LOOP            %
%****************************%

per = length(startdates_in);
if per ~= length(enddates_in); error('NO2ColMap:TimePeriods','no2_column_map - startdate and enddate are unequal lengths'); end
for period = 1:per %Loop over each temporal period you wish to average
    startdate = startdates_in{period}; enddate = enddates_in{period};
    for a = datenum(startdate,'yyyy/mm/dd'):datenum(enddate,'yyyy/mm/dd')
        date = datestr(a,29); %Convert the current date number to a string in the format yyyy-mm-dd
        year = date(1:4);
        month = date(6:7);
        day = date(9:10);
        
        filepath = fullfile(behr_dir,year,month);
        filename = [file_prefix,year,month,day,'.mat'];
        file = fullfile(filepath,filename);
        
        if exist(file,'file') ~= 2; fprintf(' %s not found\n',filename); 
        else
            load(file,'OMI')
            if period == 1 && a == datenum(startdate,'yyyy/mm/dd') %The first time through, get some information from the BEHR file
                SumWeightedColumn = zeros(size(OMI(1).BEHRColumnAmountNO2Trop));
                SumWeight = zeros(size(OMI(1).BEHRColumnAmountNO2Trop));
                if ~man_latlon %If we want to automatically assign lats and lons, do so now
                    latmin = OMI(1).MapData.LatBdy(1); latmax = OMI(1).MapData.LatBdy(2); latres = OMI(1).MapData.LatRes; 
                    lonmin = OMI(1).MapData.LonBdy(1); lonmax = OMI(1).MapData.LonBdy(2); lonres = OMI(1).MapData.LonRes;
                    lats = latmin:latres:latmax;
                    lons = lonmin:lonres:lonmax;
                end
            elseif ~man_latlon %Test is there is any difference between the OMI lat/lon grid and the existing one if automatic assignment was done.
                if length(OMI(1).Latitude)~=length(lats) || length(OMI(1).Longitude) ~= length(lons)
                    error('no2_col_map:geogrid_mismatch','%s lat/lons do not match existing grid', file);
                end
                lat_test = sum(OMI(1).Latitude - lats); lon_test = sum(OMI(1).Longitude - lons);
                if lat_test ~= 0 && lon_test ~= 0 
                    error('no2_col_map:geogrid_mismatch','%s lat/lons do not match existing grid', file);
                end
                clear lat_test lon_test
            end
            for s=1:length(OMI)
                omi = OMI(s); %We will set the area weight to 0 for any elements that should not contribute to the average
                omi.Areaweight(omi.BEHRColumnAmountNO2Trop<=0) = 0; %Do not average in negative tropospheric column densities
                omi.Areaweight(mod(omi.vcdQualityFlags,2)~=0) = 0; %If the vcdQualityFlags value is not even (least significant bit ~= 0), do not include this element
                
                if strcmpi(cloud_type,'omi'); omi.Areaweight(omi.CloudFraction > cloud_frac) = 0; 
                elseif strcmpi(cloud_type,'modis'); omi.Areaweight(omi.MODISCloud > cloud_frac) = 0;
                end %Do not include the element if the cloud fraction is greater than the allowable criteria
                
                omi.Areaweight(omi.BEHRColumnAmountNO2Trop > 1E17) = 0; %Do not include the element if the NO2 column is too great.
                hh=find(isnan(omi.BEHRColumnAmountNO2Trop)); omi.BEHRColumnAmountNO2Trop(hh)=0; omi.Areaweight(hh)=0; %Set any column NaNs to 0 and do not include them in the average
                
                xx = omi_rowanomaly(omi,input.rowanomaly); %Remove elements affected by the row anomaly.
                omi.Areaweight(xx) = 0;
                
                SumWeightedColumn = SumWeightedColumn + omi.BEHRColumnAmountNO2Trop .* omi.Areaweight;
                SumWeight = SumWeight + omi.Areaweight;
            end %End loop over swaths in OMI
        end %End if statement checking if the BEHR file for the current day exists
    end %End the loop over all days in this time period
end %End the loop over all time periods

%Normalize the areaweight for each pixel to 1.
ColumnData = SumWeightedColumn ./ SumWeight;

%Grid the data. If automatic lat/lon gridding was used, all files were
%tested to have the same lat/lon grid. If not, we will need to match up the
%disperate grids.


end