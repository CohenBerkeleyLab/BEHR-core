function [ output_args ] = BEHR_hdf_v2(  )
%BEHR_hdf_v2 Create the HDF files for BEHR products
%   Detailed explanation goes here

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SET OPTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Set to 'native' to save the native OMI resolution pixels. Set to
% 'gridded' to save the 0.05 x 0.05 gridded data

pixel_type = 'native';

% Make the list of variables to put in the HDF files. Std. variables will
% be added by default; see the "set_variables" function for additional
% options.

[vars, savename] = set_variables();
attr = add_attributes(vars);

% The dates to process, location of the files, and where to save the files.
% If you want to process all files in a directory, set the start and end
% dates to something silly.
start_date = '';
end_date = '';

mat_file_dir = '';
save_dir = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ismember(pixel_type,{'native','gridded'})
    E.badinput('"pixel_type" must be "native" or "gridded"');
end

if datenum(start_date) < datenum('2004-10-01') || datenum(end_date) < datenum('2004-10-01')
    E.badinput('start and end dates must be after Oct 1st, 2004');
elseif datenum(start_date) > datenum(end_date)
    E.badinput('Start date must be earlier than end date');
end

if ~exist(mat_file_dir,'dir')
    E.badinput('mat_file_dir must be a directory');
end

if ~exist(save_dir,'dir')
    E.badinput('save_dir must be a directory');
end

%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN LOOP %%%%%
%%%%%%%%%%%%%%%%%%%%%

FILES = dir(fullfile(mat_file_dir,'OMI_BEHR*.mat'));
for a=1:numel(FILES)
    % Find the date part of the file
    d_ind = regexp(FILES(a).name,'\d\d\d\d_\d\d_\d\d');
    date_string = FILES(a).name(d_ind:d_ind+9);
    if datenum(date_string) >= datenum(start_date) && datenum(date_string) <= datenum(end_date)
        load(fullfile(mat_file_dir,FILES(a).name));
        if strcmpi(pixel_type,'native')
            Data_to_save = Data;
        else
            Data_to_save = OMI;
        end
        
        make_hdf_file(Data_to_save,vars,attr,date_string,save_dir,savename)
    end
end

end

function [vars, savename] = set_variables(varargin)
% Make a list of variables that should be added to the product. All the
% standard variables will be added always. Pass any or all of the following
% strings to add certain variables
%
%   'reprocessed' - fields related to columns that are reprocessed using
%   in-situ profiles
%
% The standard variables to be included (listed in
% http://behr.cchem.berkeley.edu/TheBEHRProduct.aspx)

vars = {'AMFStrat','AMFTrop','BEHRAMFTrop','BEHRColumnAmountNOTrop','CloudFraction',...
    'CloudPressure','CloudRadianceFraction','ColumnAmountNO','ColumnAmountNOInitial',...
    'ColumnAmountNOStrat','ColumnAmountNOTrop','ColumnAmountNOTropStd','GLOBETerpres',...
    'Latcorn','Latitude','Loncorn','Longitude','MODISAlbedo','MODISCloud',...
    'RelativeAzimuthAngle','Row','SlantColumnAmountNO','SolarAzimuthAngle',...
    'SolarZenithAngle','Swath','TerrainHeight','TerrainPressure','TerrainReflectivity',...
    'Time','ViewingAzimuthAngle','ViewingZenithAngle','XTrackQualityFlags','vcdQualityFlags'};



% 

savename = 'OMI_BEHR_';

% The reprocessing related fields
if ismember('reprocessed',varargin)
    repro_vars = {'InSituAMF','BEHR_R_ColumnAmountNO2Trop','ProfileCount','InSituFlags'};
    vars = cat(2,vars,repro_vars);
    savename = strcat(savename,'InSitu_');
end

end

function attr = add_attributes(vars)
E = JLLErrors;

if ~~iscell(vars)
    E.badinput('vars must be a cell array');
end

attr = make_empty_struct_from_cell(vars);

% This cell array will have the variable name, unit, range, fill, product
% (SP or BEHR) and description in that order.
longfill = -1.267650600228229401496703205376e30;
shortfill = -32767;
shortfill2 = shortfill / 1000; % This should be fixed in the next run of BEHR to be the same as shortfill
behrfill = -9e9;
nofill = NaN;

attr_table = {  'AMFStrat', 'unitless', [0, Inf], nofill, 'SP', 'Stratospheric AMF';
                'AMFTrop', 'unitless', [0, Inf], nofill, 'SP', 'Tropospheric AMF (standard product)';
                'BEHRAMFTrop', 'unitless', [0, Inf], shortfill, 'BEHR', 'Tropospheric AMF (BEHR) calculated with MODIS Albedo, GLOBE Terr. Pres., and 12 km NO2 profiles';
                'BEHRColumnAmountNO2Trop', 'molec./cm^2', [0, Inf], longfill, 'BEHR', 'Tropospheric NO2 VCD (BEHR) calculated as SCD_trop / AMF_BEHR';
                'CloudFraction', 'unitless', [0, 1], shortfill2, 'SP', 'OMI geometric cloud fraction';
                'CloudPressure', 'hPa', [0, Inf], shortfill, 'SP', 'OMI cloud top pressure';
                'CloudRadianceFraction', 'unitless', [0, 1], shortfill2, 'SP', 'OMI cloud radiance (top of atmosphere light fraction)';
                'ColumnAmountNO2', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Total NO2 VCD';
                'ColumnAmountNO2Initial', 'molec./cm^2', [0, Inf], longfill, 'SP', '';
                'ColumnAmountNO2Strat', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Stratospheric NO2 VCD';
                'ColumnAmountNO2Trop', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Tropospheric NO2 VCD (standard product)';
                'ColumnAmountNO2TropStd', 'molec./cm^2', [0, Inf], nofill, 'SP', 'Standard deviation of SP NO2 tropospheric VCD';
                'GLOBETerpres', 'hPa', [0, Inf], behrfill, 'BEHR', 'Terrain pressure derived from GLOBE (1 x 1 km) topography using standard scale height relation, avg. to OMI pixel';
                'Latcorn', 'deg', [-90, 90], nofill, 'BEHR', 'Calculated corner latitude of pixels';
                'Latitude', 'deg', [-90, 90], nofill, 'SP', 'Center latitude of pixels';
                'Loncorn', 'deg', [-180, 180], nofill, 'BEHR', 'Calculated corner longitude of pixels';
                'Longitude', 'deg', [-180, 180], nofill, 'SP', 'Center longitude of pixels';
                'MODISAlbedo', 'unitless', [0, 1], behrfill, 'BEHR', 'MODIS MCD43C3 16-day avg. window every 8 days, avg. to OMI pixel';
                'MODISCloud', 'unitless', [0, 1], behrfill, 'BEHR', 'MODIS MYD06_L2 5 x 5 km cloud fraction, avg. to OMI pixel';
                'RelativeAzimuthAngle', 'deg', [0, 180], nofill, 'BEHR', 'Calculated azimuth angle between sun and satellite';
                'Row', 'unitless', [0, 59], nofill, 'SP', 'Across track row number, 0 based';
                'SlantColumnAmountNO2', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Total NO2 SCD';
                'SolarAzimuthAngle', 'deg', [-180, 180], nofill, 'SP', 'Solar azimuth angle';
                'SolarZenithAngle', 'deg', [0 90], nofill, 'SP', 'Solar zenith angle';
                'Swath', 'unitless', [0, Inf], nofill, 'SP', 'Swath number since OMI launch';
                'TerrainHeight', 'm', [-Inf, Inf], nofill, 'SP', 'Terrain height';
                'TerrainPressure', 'hPa', [0, Inf], nofill 'SP', 'Terrain pressure';
                'TerrainReflectivity', 'unitless', [0, 1], shortfill2, 'SP', 'Terrain albedo (OMI albedo product)';
                'Time', 's', [0, Inf], nofill, 'SP', 'Time at start of scan (TAI93: seconds since Jan 1, 1993)';
                'ViewingAzimuthAngle', 'deg', [-180, 180], nofill, 'SP', 'Viewing azimuth angle';
                'ViewingZenithAngle', 'deg', [0, 90], nofill, 'SP', 'Viewing zenith angle';
                'XTrackQualityFlags', 'bit array flag', 'N/A', nofill, 'SP', 'Across track quality flag (for row anomaly)';
                'vcdQualityFlags', 'bit array flag', 'N/A', nofill, 'SP', 'Ground pixel quality flags';
                'InSituAMF', 'unitless', [0, Inf], behrfill, 'BEHR-InSitu', 'AMF calculated using co-located in situ NO2 profile';
                'BEHR_R_ColumnAmountNO2Trop', 'molec./cm^2', [0, Inf], behrfill, 'BEHR-InSitu', 'BEHR Tropospheric NO2 VCD calculated with the in situ AMF';
                'ProfileCount', 'unitless', [0, Inf], behrfill, 'BEHR-InSitu', 'Number of aircraft profiles averaged to create the in situ a priori NO2 profile';
                'InSituFlags', 'bit array flag', 'N/A', nofill, 'BEHR-InSitu', 'In situ profile quality flag';
                };
            
fns = fieldnames(attr);
for a=1:numel(fns)
    xx = strcmp(fns{a}, attr_table(:,1));
    if sum(xx) == 0;
        E.callError('var_attr_def', sprintf('The attributes for the variable %s are not defined in the attr_table',fns{a}));
    elseif sum(xx) > 1
        E.callError('var_attr_def', sprintf('The attributes for the variable %s are defined more than once in the attr_table',fns{a}));
    end
    
    attr.(fns{a}).Unit = attr_table{xx,2};
    attr.(fns{a}).Range = attr_table{xx,3};
    attr.(fns{a}).Fill = attr_table{xx,4};
    attr.(fns{a}).Product = attr_table{xx,5};
    attr.(fns{a}).Description = attr_table{xx,6};
end

end


function make_hdf_file(Data_in, vars, attr, date_string, save_dir, savename)
E = JLLErrors;
global ask_to_overwrite
if isempty(ask_to_overwrite)
    ask_to_overwrite = true;
end

if ~strcmp(savename(end),'_')
    savename = strcat(savename,'_');
end
hdf_filename = strcat(savename, date_string, '.mat');
hdf_fullfilename = fullfile(save_dir, hdf_filename);

% Check if the file exists. Give the user 3 options if it does: abort,
% overwrite, overwrite all.

if exist(hdf_fullfilename,'file') && ask_to_overwrite;
    if ask_to_overwrite
        user_ans = input(sprintf('File %s exists. [O]verwrite, [A]bort, or Overwrite and [d]on''t ask again? ',hdf_fullfilename));
        user_ans = lower(user_ans);
        switch user_ans
            case 'o'
                delete(hdf_fullfilename);
            case 'd'
                delete(hdf_fullfilename);
                ask_to_overwrite = false;
            otherwise
                E.userCancel;
        end
    else
        delete(hdf_fullfilename);
    end
end

for d=1:numel(Data_in)
    swath_id = max(Data_in(d).Swath(:));
    group_name = sprintf('/Data/Swath%05d',swath_id);
    for v=1:numel(vars)
        var_name = sprintf('%s/%s',group_name,vars{v});
        save_data = Data(d).(vars{v});
        sz = size(save_data);
        
        % Make NaNs into fill values - apparently this is better for HDF
        % type files
        nans = isnan(save_data);
        save_data(nans) = attr.(vars{v}).Fill;
        
        if isa(save_data,'double')
            % Convert doubles to singles to save space for the data people
            % will be downloading
            save_data = single(save_data);
        end
        
        % Create the dataset, then write it and add the attributes
        h5create(hdf_fullfilename, var_name, sz, 'Datatype', class(save_data), 'FillValue', attr.(vars{v}).Fill);
        h5write(hdf_fullfilename, var_name, save_data);
        
        atts = fieldnames(attr.(vars{v}));
        for a=1:numel(atts);
            if strcmp(atts{a},'Fill')
                continue % We've already handled the fill value with h5create
            end
            h5writeatt(hdf_fullfilename, var_name, atts{a}, attr.(vars{v}).(atts{a}));
        end
        
    end
end
end