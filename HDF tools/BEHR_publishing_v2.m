function [  ] = BEHR_publishing_v2(output_type, pixel_type, options, start_date, end_date)
%BEHR_publishing_v2 Create the HDF files for BEHR products
%   Detailed explanation goes here
%
%   Non-built in dependencies (updated 22 Sept 2015):
%       Classes/JLLErrors.m
%       Utils/bitopmat.m
%       Utils/iscellcontents.m
%       Utils/make_empty_struct_from_cell.m

global onCluster
if isempty(onCluster)
    onCluster = false;
end

if onCluster
    addpath('~/MATLAB/Utils');
    addpath('~/MATLAB/Classes');
    
    % Cleanup object will safely exit if there's a problem
    cleanupobj = onCleanup(@() mycleanup());
end

E = JLLErrors;
global DEBUG_LEVEL
DEBUG_LEVEL = 1;


%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SET OPTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Output type should be 'txt' or 'hdf'.  Text (csv) files are for native
% resolution only.
if ~exist('output_type','var')
    output_type = 'txt';
else
    allowed_outtype = {'txt','hdf'};
    if ~ismember(output_type,allowed_outtype)
        E.badinput('output_type must be one of %s',strjoin(allowed_outtype,', '));
    end
end

% Set to 'native' to save the native OMI resolution pixels. Set to
% 'gridded' to save the 0.05 x 0.05 gridded data
if ~exist('pixel_type','var')
    pixel_type = 'native';
else
    allowed_pixtype = {'native','gridded'};
    if ~ismember(pixel_type,allowed_pixtype)
        E.badinput('pixel_type must be one of %s',strjoin(allowed_pixtype,', '));
    end
end

% options - add 'reprocessed' here if doing in situ files
if ~exist('options','var')
    options = {};
end

% Make the list of variables to put in the HDF files. Std. variables will
% be added by default; see the "set_variables" function for additional
% options. The pixel type needs to be passed so that it knows whether to
% keep the pixel specific variables or not.

[vars, savename] = set_variables(pixel_type, output_type, options{:});
attr = add_attributes(vars);

% The dates to process, location of the files, and where to save the files.
% If you want to process all files in a directory, set the start and end
% dates to something silly.
if ~exist('start_date','var') || ~exist('end_date','var')
    start_date = '2008-02-18';
    end_date = '2016-01-01';
end

global mat_file_dir
global save_dir
global numThreads
if ~onCluster
    mat_file_dir = BEHR_paths('behr_mat_dir');
    save_subdir = sprintf('behr_%s-%s_%s',pixel_type,output_type,BEHR_version);
    save_dir = fullfile(BEHR_paths('website_staging_dir'),save_subdir);
    if ~exist(save_dir,'dir')
        mkdir(save_dir)
    end
else
    % Check that all global variables are set
    global_unset = {};
    if isempty(mat_file_dir)
        global_unset{end+1} = 'mat_file_dir';
    end
    if isempty(save_dir)
        global_unset{end+1} = 'save_dir';
    end
    if isempty(numThreads)
        global_unset{end+1} = 'numThreads';
    end
    if ~isempty(global_unset)
        E.runscript_error(global_unset);
    end
    
    % Check that both directories exist and that numThreads is the proper
    % type
    if ~isnumeric(numThreads) || ~isscalar(numThreads)
        E.badinput('numThreads should be a scalar number; this is a global setting, check the calling runscript')
    end
    
    dirs_dne = {};
    if ~exist(mat_file_dir,'dir')
        dirs_dne{end+1} = 'mat_file_dir';
    end
    if ~exist(save_dir,'dir')
        dirs_dne{end+1} = 'save_dir';
    end
    if ~isempty(dirs_dne)
        E.dir_dne(dirs_dne);
    end
   
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(output_type,'txt') && strcmpi(pixel_type, 'gridded')
    E.badinput('Gridded output is only intended for HDF files')
end

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

% Split into two different loops: if running on a cluster, it will
% parallelize and assume that you want to overwrite any existing files. If
% running locally, it will not parallelize, and will ask for your decision
% on overwriting files.

FILES = dir(fullfile(mat_file_dir,'OMI_BEHR*.mat'));
if ~onCluster
    ask_to_overwrite = true;
    for a=1:numel(FILES)
        % Find the date part of the file
        d_ind = regexp(FILES(a).name,'\d\d\d\d\d\d\d\d');
        date_string = FILES(a).name(d_ind:d_ind+7);
        if datenum(date_string,'yyyymmdd') >= datenum(start_date) && datenum(date_string,'yyyymmdd') <= datenum(end_date)
            % If the file is less than a MB, it likely has no data (possibly
            % because the OMI swaths needed were not created for that day). If this
            % is true, skip this file.
            if FILES(a).bytes < 1e6
                if DEBUG_LEVEL > 0; fprintf('%s size < 1 MB, skipping due to lack of data\n',FILES(a).name); end
                continue
            end

            load(fullfile(mat_file_dir,FILES(a).name));
            if strcmpi(pixel_type,'native')
                Data_to_save = Data;
            else
                Data_to_save = OMI;
            end

            if DEBUG_LEVEL > 0
                fprintf('Saving %s %s for %s\n', pixel_type, output_type, date_string);
            end

            if strcmpi(output_type,'hdf')
                ask_to_overwrite = make_hdf_file(Data_to_save,vars,attr,date_string,save_dir,savename,pixel_type,ask_to_overwrite);
            elseif strcmpi(output_type,'txt')
                ask_to_overwrite = make_txt_file(Data_to_save,vars,attr,date_string,save_dir,savename,ask_to_overwrite);
            end
        end
    end
else
    if onCluster && isempty(gcp('nocreate'))
        parpool(numThreads);
    end
    ask_to_overwrite = false;
    parfor a=1:numel(FILES)
        % Find the date part of the file
        d_ind = regexp(FILES(a).name,'\d\d\d\d\d\d\d\d');
        date_string = FILES(a).name(d_ind:d_ind+7);
        if datenum(date_string,'yyyymmdd') >= datenum(start_date) && datenum(date_string,'yyyymmdd') <= datenum(end_date)
            % If the file is less than a MB, it likely has no data (possibly
            % because the OMI swaths needed were not created for that day). If this
            % is true, skip this file.
            if FILES(a).bytes < 1e6
                if DEBUG_LEVEL > 0; fprintf('%s size < 1 MB, skipping due to lack of data\n',FILES(a).name); end
                continue
            end

            D = load(fullfile(mat_file_dir,FILES(a).name));
            if strcmpi(pixel_type,'native')
                Data_to_save = D.Data;
            else
                Data_to_save = D.OMI;
            end

            if DEBUG_LEVEL > 0
                fprintf('Saving %s %s for %s\n', pixel_type, output_type, date_string);
            end

            if strcmpi(output_type,'hdf')
                make_hdf_file(Data_to_save,vars,attr,date_string,save_dir,savename,pixel_type,ask_to_overwrite);
            elseif strcmpi(output_type,'txt')
                make_txt_file(Data_to_save,vars,attr,date_string,save_dir,savename,ask_to_overwrite);
            end
        end
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

vars = {'AMFStrat','AMFTrop','BEHRAMFTrop','BEHRColumnAmountNO2Trop',...
    'BEHRPressureLevels','CloudFraction','CloudPressure','CloudRadianceFraction','ColumnAmountNO2',...
    'ColumnAmountNO2Trop','ColumnAmountNO2TropStd','ColumnAmountNO2Strat','GLOBETerpres',...
    'Latcorn','Latitude','Loncorn','Longitude','MODISAlbedo','MODISCloud',...
    'RelativeAzimuthAngle','Row','SlantColumnAmountNO2','SolarAzimuthAngle',...
    'SolarZenithAngle','Swath','TerrainHeight','TerrainPressure','TerrainReflectivity',...
    'Time','ViewingAzimuthAngle','ViewingZenithAngle','XTrackQualityFlags','vcdQualityFlags',...
    'BEHRScatteringWeights','BEHRAvgKernels','BEHRNO2apriori','BEHRGhostFraction'};



% 

savename = sprintf('OMI_BEHR_%s_',BEHR_version);

% Add additional variable categories here. You'll need to add the variable
% names to the "vars" variable, and you should consider adding an
% identifier to the save name to make clear what variables are present.
% Note two other places you'll need to add variables: in the subfunction
% "remove_ungridded_variables" if any of the new variables should be
% included in the gridded products and in the "add_attributes" subfunction
% - there you'll want to include information like unit, range, fill,
% product, and description.

% The reprocessing related fields
if ismember('reprocessed',varargin)
    repro_vars = {'InSituAMF','BEHR_R_ColumnAmountNO2Trop','ProfileCount','InSituFlags'};
    vars = cat(2,vars,repro_vars);
    savename = strcat(savename,'InSitu_');
end

% Remove pixel specific variables (like AMF, VZA, etc.) if the pixel type
% is "gridded". Also add in some variables that are only included in the
% gridded product.
if ismember('gridded', varargin)
    vars = remove_ungridded_variables(vars);
    vars{end+1} = 'Areaweight';
end

% Remove variables that cannot be put into a CSV text file because multiple
% values are required per pixel
if ismember('txt', varargin)
    vars = remove_vector_variables(vars);
end


end

function vars = remove_ungridded_variables(vars)
E = JLLErrors;
if ~iscell(vars) || ~all(iscellcontents(vars,'ischar'))
    E.badinput('"vars" must be a cell array of variable names as strings')
end

% Define what variables should be included in gridded products. You'll need
% to edit this is you add new gridded variables.
gridded_vars = {'AMFTrop', 'Areaweight', 'BEHRAMFTrop','BEHRColumnAmountNO2Trop',...
    'CloudFraction', 'CloudRadianceFraction', 'ColumnAmountNO2Trop', 'GLOBETerpres',...
    'Latitude', 'Longitude', 'MODISAlbedo', 'MODISCloud', 'Row', 'XTrackQualityFlags',...
    'vcdQualityFlags', 'InSituAMF', 'BEHR_R_ColumnAmountNO2Trop'};

gg = ismember(vars,gridded_vars);
vars = vars(gg);

end

function vars = remove_vector_variables(vars)
E=JLLErrors;
if ~iscell(vars) || ~all(iscellcontents(vars,'ischar'))
    E.badinput('"vars" must be a cell array of variable names as strings')
end

% Define what variables cannot be saved in a text file and remove them
vector_vars = {'BEHRPressureLevels','BEHRScatteringWeights','BEHRAvgKernels','BEHRNO2apriori','BEHRGhostFraction'};
vv = ismember(vars, vector_vars);
vars = vars(~vv);
end

function attr = add_attributes(vars)
E = JLLErrors;

if ~iscell(vars)
    E.badinput('vars must be a cell array');
end

attr = make_empty_struct_from_cell(vars);

% This cell array will have the variable name, unit, range, fill, product
% (SP or BEHR) and description in that order.
longfill = single(-1.267650600228229401496703205376e30);
shortfill = single(-32767);
behrfill = single(-3.402e38);
nofill = NaN;

attr_table = {  'AMFStrat', 'unitless', [0, Inf], nofill, 'SP', 'Stratospheric AMF';...
                'AMFTrop', 'unitless', [0, Inf], nofill, 'SP', 'Tropospheric AMF (standard product)';...
                'Areaweight', 'unitless', [0, Inf], behrfill, 'BEHR', 'Reciprocal of pixel area; use to weight a temporal average of grids to account for pixel representativeness';...
                'BEHRAMFTrop', 'unitless', [0, Inf], behrfill, 'BEHR', 'Tropospheric AMF (BEHR) calculated with MODIS Albedo, GLOBE Terr. Pres., and 12 km NO2 profiles';...
                'BEHRColumnAmountNO2Trop', 'molec./cm^2', [0, Inf], behrfill, 'BEHR', 'Tropospheric NO2 VCD (BEHR) calculated as SCD_trop / AMF_BEHR';...
                'BEHRScatteringWeights', 'unitless', [0, Inf], behrfill, 'BEHR', 'Scattering weights derived from the MODIS albedo and GLOBE surface pressure. Includes NO2 cross section temperature correction.';...
                'BEHRAvgKernels', 'unitless', [0, Inf], behrfill, 'BEHR', 'Averaging kernels computed for the weighted average of cloudy and clear conditions';...
                'BEHRPressureLevels', 'hPa', [0, Inf], behrfill, 'BEHR', 'Pressure levels that correspond to the scattering weight, averaging kernel, and NO2 a priori vectors';...
                'CloudFraction', 'unitless', [0, 1], shortfill, 'SP', 'OMI geometric cloud fraction';...
                'CloudPressure', 'hPa', [0, Inf], shortfill, 'SP', 'OMI cloud top pressure';...
                'CloudRadianceFraction', 'unitless', [0, 1], shortfill, 'SP', 'OMI cloud radiance (top of atmosphere light fraction)';...
                'ColumnAmountNO2', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Total NO2 VCD';...
                'ColumnAmountNO2Trop', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Tropospheric NO2 VCD (standard product)';...
                'ColumnAmountNO2TropStd', 'molec./cm^2', [0, Inf], nofill, 'SP', 'Standard deviation of SP NO2 tropospheric VCD';...
                'ColumnAmountNO2Strat', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Stratospheric NO2 VCD';...
                'GLOBETerpres', 'hPa', [0, Inf], behrfill, 'BEHR', 'Terrain pressure derived from GLOBE (1 x 1 km) topography using standard scale height relation, avg. to OMI pixel';...
                'Latcorn', 'deg', [-90, 90], nofill, 'BEHR', 'Calculated corner latitude of pixels';...
                'Latitude', 'deg', [-90, 90], nofill, 'SP', 'Center latitude of pixels';...
                'Loncorn', 'deg', [-180, 180], nofill, 'BEHR', 'Calculated corner longitude of pixels';...
                'Longitude', 'deg', [-180, 180], nofill, 'SP', 'Center longitude of pixels';...
                'MODISAlbedo', 'unitless', [0, 1], behrfill, 'BEHR', 'MODIS MCD43C3 16-day avg. window every 8 days, avg. to OMI pixel';...
                'MODISCloud', 'unitless', [0, 1], behrfill, 'BEHR', 'MODIS MYD06_L2 5 x 5 km cloud fraction, avg. to OMI pixel';...
                'RelativeAzimuthAngle', 'deg', [0, 180], nofill, 'BEHR', 'Calculated azimuth angle between sun and satellite';...
                'Row', 'unitless', [0, 59], nofill, 'SP', 'Across track row number, 0 based';...
                'SlantColumnAmountNO2', 'molec./cm^2', [0, Inf], longfill, 'SP', 'Total NO2 SCD';...
                'SolarAzimuthAngle', 'deg', [-180, 180], nofill, 'SP', 'Solar azimuth angle';...
                'SolarZenithAngle', 'deg', [0 90], nofill, 'SP', 'Solar zenith angle';...
                'Swath', 'unitless', [0, Inf], nofill, 'SP', 'Swath number since OMI launch';...
                'TerrainHeight', 'm', [-Inf, Inf], nofill, 'SP', 'Terrain height';...
                'TerrainPressure', 'hPa', [0, Inf], nofill 'SP', 'Terrain pressure';...
                'TerrainReflectivity', 'unitless', [0, 1], shortfill, 'SP', 'Terrain albedo (OMI albedo product)';...
                'Time', 's', [0, Inf], nofill, 'SP', 'Time at start of scan (TAI93: seconds since Jan 1, 1993)';...
                'ViewingAzimuthAngle', 'deg', [-180, 180], nofill, 'SP', 'Viewing azimuth angle';...
                'ViewingZenithAngle', 'deg', [0, 90], nofill, 'SP', 'Viewing zenith angle';...
                'XTrackQualityFlags', 'bit array flag', 'N/A', nofill, 'SP', 'Across track quality flag (for row anomaly)';...
                'vcdQualityFlags', 'bit array flag', 'N/A', nofill, 'SP', 'Ground pixel quality flags';...
                'InSituAMF', 'unitless', [0, Inf], behrfill, 'BEHR-InSitu', 'AMF calculated using co-located in situ NO2 profile';...
                'BEHR_R_ColumnAmountNO2Trop', 'molec./cm^2', [0, Inf], behrfill, 'BEHR-InSitu', 'BEHR Tropospheric NO2 VCD calculated with the in situ AMF';...
                'ProfileCount', 'unitless', [0, Inf], behrfill, 'BEHR-InSitu', 'Number of aircraft profiles averaged to create the in situ a priori NO2 profile';...
                'InSituFlags', 'bit array flag', 'N/A', nofill, 'BEHR-InSitu', 'In situ profile quality flag';...
                'BEHRColumnAmountNO2Trop_L3', 'molec./cm^2', [0, Inf], behrfill, 'BEHR-L3','BEHR tropospheric NO2 VCDs filtered for quality and row anomaly';...
                'BEHRColumnAmountNO2Trop_L3MODISCloud', 'molec./cm^2', [0, Inf], behrfill, 'BEHR-L3','BEHR tropospheric NO2 VCDs additionally filtered for MODIS Cloud < 20%';...
                'BEHRNO2apriori', 'parts-per-part', [-Inf, Inf], behrfill, 'BEHR', 'NO2 a priori profile used for each pixel. Pressure levels given in BEHRPressureLevels.';...
                'BEHRGhostFraction', 'unitless', [0, Inf], behrfill, 'BEHR', 'Multiply the VCD by this value to estimate the full column (visible + ghost).'...
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


function ask_to_overwrite = make_hdf_file(Data_in, vars, attr, date_string, save_dir, savename, pixel_type, ask_to_overwrite)
E = JLLErrors;

global DEBUG_LEVEL
if isempty(DEBUG_LEVEL)
    DEBUG_LEVEL = 0;
end

if ~strcmp(savename(end),'_')
    savename = strcat(savename,'_');
end
hdf_filename = strcat(savename, date_string, '.hdf');
hdf_fullfilename = fullfile(save_dir, hdf_filename);

% Check if the file exists. Give the user 3 options if it does: abort,
% overwrite, overwrite all.

if exist(hdf_fullfilename,'file')
    if ask_to_overwrite
        user_ans = input(sprintf('File %s exists.\n[O]verwrite, [A]bort, or Overwrite and [d]on''t ask again? ',hdf_fullfilename),'s');
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
        %fprintf('%s exists, skipping\n', hdf_fullfilename);
        %return
    end
end



% Iterate through each swath and save it as under the group
% /Data/Swath#####.
for d=1:numel(Data_in)
    swath_id = max(Data_in(d).Swath(:));
    group_name = sprintf('/Data/Swath%d',swath_id);
    
    if swath_id == 0 || isnan(swath_id)
        % A swath ID of 0 or NaN means that no pixels were gridded, so skip
        % this swath since has it has no useful data.
        continue
    end
    
    if DEBUG_LEVEL > 1; fprintf('\t Now writing %s\n',group_name); end
    if DEBUG_LEVEL > 2; tic; end
    for v=1:numel(vars)
        var_name = sprintf('%s/%s',group_name,vars{v});
        save_data = Data_in(d).(vars{v});
        sz = size(save_data);
        
        % Make NaNs into fill values - apparently this is better for HDF
        % type files. Cell arrays are created for bit array flags in
        % gridded products - so we don't want to do any fill values for
        % them. Rather, since each cell contains a matrix of the bit array
        % flags, we'll do a bitwise OR operation on them so that if a flag
        % is set for any pixel used in that grid cell it carries through.
        if ~iscell(save_data)
            nans = isnan(save_data);
            save_data(nans) = attr.(vars{v}).Fill;
        else
            save_data_cell = save_data;
            save_data = uint16(zeros(sz));
            for c=1:numel(save_data_cell)
                if DEBUG_LEVEL > 3 && mod(c,10000)==1; fprintf('Cell %d of %d\n',c,numel(save_data_cell)); end
                flag = uint16(save_data_cell{c}(:));
                % An empty matrix in the cell means there was no data
                % there, so just assign it a value of 0.
                if isempty(flag); 
                    save_data(c) = uint16(0); 
                else
                    save_data(c) = bitopmat(flag,'or');
                end
            end
        end
        
        if isa(save_data,'double')
            % Convert doubles to singles to save space for the data people
            % will be downloading
            save_data = single(save_data);
        end
        % Ensure that the fill value is of the same type as the data
        fill_val = cast(attr.(vars{v}).Fill, 'like', save_data);
        
        % Create the dataset, then write it and add the attributes
        h5create(hdf_fullfilename, var_name, sz, 'Datatype', class(save_data), 'FillValue', fill_val);
        h5write(hdf_fullfilename, var_name, save_data);
        
        atts = fieldnames(attr.(vars{v}));
        for a=1:numel(atts);
            if strcmp(atts{a},'Fill')
                continue % We've already handled the fill value with h5create
            end
            h5writeatt(hdf_fullfilename, var_name, atts{a}, attr.(vars{v}).(atts{a}));
        end
        
    end
    
    % Write an attribute to the swath group describing if it is gridded or
    % native. Also include the version string
    switch lower(pixel_type)
        case 'native'
            swath_attr = 'OMI SP and BEHR data at native OMI resolution';
        case 'gridded'
            swath_attr = 'OMI SP and BEHR data gridded to 0.05 x 0.05 deg';
        otherwise
            E.badinput('"pixel_type" not recognized');
    end
    h5writeatt(hdf_fullfilename,group_name,'Description',swath_attr);
    h5writeatt(hdf_fullfilename,group_name,'Version',BEHR_version());
    if DEBUG_LEVEL > 2; toc; end
end
end



function ask_to_overwrite = make_txt_file(Data_in, vars, attr, date_string, save_dir, savename, ask_to_overwrite)
E = JLLErrors;

global DEBUG_LEVEL
if isempty(DEBUG_LEVEL)
    DEBUG_LEVEL = 0;
end

if ~strcmp(savename(end),'_')
    savename = strcat(savename,'_');
end
txt_filename = strcat(savename, date_string, '.txt');
txt_fullfilename = fullfile(save_dir, txt_filename);

% Check if the file exists. Give the user 3 options if it does: abort,
% overwrite, overwrite all.

if exist(txt_fullfilename,'file')
    if ask_to_overwrite
        user_ans = input(sprintf('File %s exists.\n[O]verwrite, [A]bort, or Overwrite and [d]on''t ask again? ',txt_fullfilename),'s');
        user_ans = lower(user_ans);
        switch user_ans
            case 'o'
                delete(txt_fullfilename);
            case 'd'
                delete(txt_fullfilename);
                ask_to_overwrite = false;
            otherwise
                E.userCancel;
        end
    else
        delete(txt_fullfilename);
    end
end

% For text files, we will not break up by swaths, instead all pixels will
% be in one giant CSV type output.

% First we'll create the format string based on the variables requested.
% Most variables will have 6 significant digits, using %g (so exponential
% or standard form will be chosen for compactness). Some will be specified
% to be integers - either if the class of the value is an integer or it is
% a flag field. Time will be treated specially because we want very high
% precision, and the Lat/Loncorn fields will need to be expanded into four
% individual fields. Next the header - start with lon, lat, and the
% corners. The order of the rest is less important.
n_vars = numel(vars);
header_cell = cell(1,n_vars+6);
header_cell(1:10) = {'Longitude','Latitude','Loncorn1','Loncorn2','Loncorn3','Loncorn4','Latcorn1','Latcorn2','Latcorn3','Latcorn4'};
format_spec = cell(1,n_vars+6);
format_spec(1:10) = repmat({'%.4f'},1,10);
i=11;
for a=1:n_vars
    if ~ismember(vars{a}, {'Longitude','Latitude','Loncorn','Latcorn'});
        header_cell{i} = vars{a};
        if strcmpi(vars{a},'Time')
            format_spec{i} = '%f';
        elseif isinteger(Data_in(1).(vars{a})(1)) || ~isempty(regexpi(vars{a},'Flag')) || any(strcmpi(vars{a},{'Row','Swath'}))
            format_spec{i} = '%d';
        else
            format_spec{i} = '%.4g';
        end
        i=i+1;
    end
end
header_line = strjoin(header_cell,',');

% Open the file and loop through all the swaths and pixels and
% write the values.

fid = fopen(txt_fullfilename,'w');
fprintf(fid,'%s\n',header_line);

for s=1:numel(Data_in)
    for p=1:numel(Data_in(s).Longitude)
        for a=1:numel(header_cell)
            switch header_cell{a}
                case 'Loncorn1'
                    v = Data_in(s).Loncorn(1,p);
                    if isnan(v)
                        v = attr.Loncorn.Fill;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Loncorn2'
                    v = Data_in(s).Loncorn(2,p);
                    if isnan(v)
                        v = attr.Loncorn.Fill;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Loncorn3'
                    v = Data_in(s).Loncorn(3,p);
                    if isnan(v)
                        v = attr.Loncorn.Fill;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Loncorn4'
                    v = Data_in(s).Loncorn(4,p);
                    if isnan(v)
                        v = attr.Loncorn.Fill;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Latcorn1'
                    v = Data_in(s).Latcorn(1,p);
                    if isnan(v)
                        v = attr.Latcorn.Fill;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Latcorn2'
                    v = Data_in(s).Latcorn(2,p);
                    if isnan(v)
                        v = attr.Latcorn.Fill;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Latcorn3'
                    v = Data_in(s).Latcorn(3,p);
                    if isnan(v)
                        v = attr.Latcorn.Fill;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Latcorn4'
                    v = Data_in(s).Latcorn(4,p);
                    if isnan(v)
                        v = attr.Latcorn.Fill;
                    end
                    fprintf(fid,format_spec{a},v);
                otherwise
                    v = Data_in(s).(header_cell{a})(p);
                    if isnan(v)
                        v = attr.(header_cell{a}).Fill;
                    end
                    fprintf(fid,format_spec{a},v);
            end
            
            if a<numel(header_cell)
                fprintf(fid,',');
            else
                fprintf(fid,'\n');
            end
        end
    end
end

fclose(fid);

end

function mycleanup()
err=lasterror;
if ~isempty(err.message)
    fprintf('MATLAB exiting due to problem: %s\n', err.message);
    if ~isempty(gcp('nocreate'))
        delete(gcp)
    end 

    exit(1)
end
end
