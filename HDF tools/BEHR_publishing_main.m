function [  ] = BEHR_publishing_main(varargin)
%BEHR_publishing_v2 Create the HDF files for BEHR products
%   BEHR_Publishing_v2 can accept a number of input parameters to alter its
%   behavior. All of these have default values that are set up so that
%   calling it without parameters will lead to standard behavior. The
%   parameters are:
%
%       'start': the first date to process, as a date number or a string
%       implicitly understood by datenum(). Default is '2005-01-01'
%
%       'end': the last date to process; same format requirements as
%       'start'. Default is today.
%
%       'output_type': one of the strings 'hdf' or 'txt', determines which
%       output format will be used. 'hdf' will produce HDF version 5 files.
%       Default is 'hdf'
%
%       'pixel_type': one of the strings 'native' or 'gridded', determines
%       whether the native pixels (i.e. the 'Data' structure) or the
%       gridded pixel (i.e. the 'OMI' structure) will be saved. Default is
%       'native'.
%
%       'reprocessed': a boolean (true or false). If true, this tells the
%       publishing algorithm to include fields that used in situ
%       measurements from the DISCOVER-AQ campaign as a priori profiles.
%       That is a specialized product that hasn't been updated in years.
%       Default is false.
%
%       'region': a string indicating which region to publish, must match
%       the directory structure in behr_paths.behr_mat_dir. Only used if
%       "mat_dir" is not specified. Default is 'us'.
%
%       'profile_mode': a string which a priori profiles' retrieval to use,
%       must match the directory structure in behr_paths.behr_mat_dir
%       (within each region). Only used if "mat_dir" is not specified.
%       Default is 'monthly'.
%
%       'mat_dir': the directory from which to load the Matlab files with
%       BEHR output saved in the. If not given (or given as an empty
%       string) then behr_paths.BEHRMatSubdir(region, profile_mode) is
%       called with the values of the "region" and "profile_mode"
%       parameters to determine the directory for the given region and
%       profile mode. Otherwise, BEHR .mat files are read from the
%       directory given.
%
%       'save_dir': the directory to which to save the resulting HDF or CSV
%       files. Default is the value returned by
%       behr_paths.website_staging_dir.
%
%       'organize': a boolean that indicates whether the output should go
%       directly in the save directory (false) or in a subdirectory named
%       behr_<pixel_type>-<output_type>-<behr version>, e.g.
%       "behr_native-hdf-v2-1C". Default is true.
%
%       'overwrite': controls the behavior of this function if one of the
%       files it is trying to output already exists. This is a number, a
%       negative value will cause it to ask you on at least the first file,
%       whether it continues to ask on successive files depends on your
%       response. 0 means do not overwrite, and a positive value means
%       always overwrite. Default is -1, i.e. ask the user.
%
%       'DEBUG_LEVEL': a scalar number that controls the verbosity of this
%       function. 0 is minimum verbosity, higher numbers print more to the
%       screen.
%
%   This function can also be parallelized using the global variables
%   numThreads and onCluster.

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

p = inputParser;
p.addParameter('output_type', 'hdf');
p.addParameter('pixel_type', 'native');
p.addParameter('start', '2005-01-01');
p.addParameter('end', datestr(today, 'yyyy-mm-dd'));
p.addParameter('reprocessed', false);
p.addParameter('mat_dir', '');
p.addParameter('region', 'us');
p.addParameter('profile_mode', 'monthly');
p.addParameter('save_dir', behr_paths.website_staging_dir);
p.addParameter('organize', true);
p.addParameter('overwrite', -1);
p.addParameter('DEBUG_LEVEL', 1);

p.parse(varargin{:});
pout = p.Results;

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SET OPTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

% Start and end date
start_date = pout.start;
end_date = pout.end;

% Output type should be 'txt' or 'hdf'.  Text (csv) files are for native
% resolution only.
output_type = pout.output_type;

allowed_outtype = {'txt','hdf'};
if ~ismember(output_type,allowed_outtype)
    E.badinput('output_type must be one of %s',strjoin(allowed_outtype,', '));
end


% Set to 'native' to save the native OMI resolution pixels. Set to
% 'gridded' to save the 0.05 x 0.05 gridded data
pixel_type = pout.pixel_type;

allowed_pixtype = {'native','gridded'};
if ~ismember(pixel_type,allowed_pixtype)
    E.badinput('pixel_type must be one of %s',strjoin(allowed_pixtype,', '));
end

% Other options - reprocessed should be TRUE to include in situ fields
is_reprocessed = pout.reprocessed;
if ~isscalar(is_reprocessed) || ~islogical(is_reprocessed)
    E.badinput('REPROCESSED must be a scalar logical')
end

% Whether subdirectories should be created within the save directory,
% organized by pixel type, output type, and BEHR version.
organized_subdir = pout.organize;
if ~isscalar(organized_subdir) || ~islogical(organized_subdir)
    E.badinput('ORGANIZE must be a scalar logical')
end

% How to handle overwriting. 1 = overwrite, 0 = don't overwrite, -1 = ask.
overwrite = pout.overwrite;

DEBUG_LEVEL = pout.DEBUG_LEVEL;
if ~isscalar(DEBUG_LEVEL) || ~isnumeric(DEBUG_LEVEL)
    E.badinput('DEBUG_LEVEL must be a scalar number')
end

% File locations

mat_file_dir = pout.mat_dir;
save_dir = pout.save_dir;

if isempty(mat_file_dir)
    mat_file_dir = behr_paths.BEHRMatSubdir(pout.region, pout.profile_mode);
end

% Check that the directories exist like this so that a single error message
% describes if both directories don't exist - handy for running on the
% cluster so that you don't wait forever for the job to start, only to have
% it fail b/c you forgot to make the output directory.
dirs_dne = {};
if ~exist(mat_file_dir,'dir')
    dirs_dne{end+1} = sprintf('mat_file_dir (%s)', mat_file_dir);
end
if ~exist(save_dir,'dir')
    dirs_dne{end+1} = sprintf('save_dir (%s)', save_dir);
end
if ~isempty(dirs_dne)
    E.dir_dne(dirs_dne);
end



% Make the list of variables to put in the HDF files. Std. variables will
% be added by default; see the "set_variables" function for additional
% options. The pixel type needs to be passed so that it knows whether to
% keep the pixel specific variables or not.


global numThreads
if onCluster
    global_unset = {};
    % Check that all global variables are set
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

git_heads.core = git_head_hash(behr_paths.behr_core);
git_heads.behr_utils = git_head_hash(behr_paths.behr_utils);
git_heads.gen_utils = git_head_hash(behr_paths.utils);

%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN LOOP %%%%%
%%%%%%%%%%%%%%%%%%%%%

% Split into two different loops: if running on a cluster, it will
% parallelize and assume that you want to overwrite any existing files. If
% running locally, it will not parallelize, and will ask for your decision
% on overwriting files.

FILES = dir(fullfile(mat_file_dir,'OMI_BEHR*.mat'));
if ~onCluster
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

            [vars, attr, savename, full_save_dir] = set_variables(Data_to_save, pixel_type, output_type, is_reprocessed, save_dir, organized_subdir);

            if strcmpi(output_type,'hdf')
                overwrite = make_hdf_file(Data_to_save, vars, attr, full_save_dir, savename, pixel_type, overwrite, git_heads, DEBUG_LEVEL);
            elseif strcmpi(output_type,'txt')
                overwrite = make_txt_file(Data_to_save,vars,attr,date_string,full_save_dir,savename,overwrite,DEBUG_LEVEL);
            end
        end
    end
else
    if onCluster && isempty(gcp('nocreate'))
        parpool(numThreads);
    end
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

            [vars, attr, savename, full_save_dir] = set_variables(Data_to_save, pixel_type, output_type, is_reprocessed, save_dir, organized_subdir);

            if strcmpi(output_type,'hdf')
                make_hdf_file(Data_to_save, vars, attr, full_save_dir, savename, pixel_type, overwrite, git_heads, DEBUG_LEVEL);
            elseif strcmpi(output_type,'txt')
                make_txt_file(Data_to_save, vars, attr, date_string, full_save_dir, savename, overwrite, DEBUG_LEVEL);
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%
% SUBFUNCTIONS %
%%%%%%%%%%%%%%%%

function [vars, attr, savename, save_dir] = set_variables(Data, pixel_type, output_type, reprocessed, save_dir, organized_subdir)
% Make a list of variables that should be added to the product. All the
% standard variables will be added always. Pass any or all of the following
% strings to add certain variables
%
%   'reprocessed' - fields related to columns that are reprocessed using
%   in-situ profiles
%
% The standard variables to be included (listed in
% http://behr.cchem.berkeley.edu/TheBEHRProduct.aspx)
E = JLLErrors;

if organized_subdir
    save_subdir = sprintf('behr-%s-%s_%s-%s_%s', Data(1).BEHRProfileMode, Data(1).BEHRRegion, pixel_type, output_type, BEHR_version);
    save_dir = fullfile(save_dir, save_subdir);
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
end

% The attribute table contains the list of all variables we expect to provide.
% Choose the proper subset.
savename = behr_filename(Data(1).Date, Data(1).BEHRProfileMode, Data(1).BEHRRegion, output_type);
if reprocessed
    attr =  BEHR_publishing_attribute_table('pub-insitu', 'struct');
    savename = strrep(savename, 'BEHR','BEHR-InSitu');
else
    attr = BEHR_publishing_attribute_table('pub', 'struct');
end

vars = fieldnames(attr);

% Remove variables not defined as "gridded" by BEHR_publishing_gridded_fields
if strcmpi('gridded', pixel_type)
    vars = remove_ungridded_variables(vars);
end

% Remove any fields not present in the structure. If one of the expected variables is
% not present, error (because we expect it to be there!) The exception are the PSM weight
% fields, they will not be present if we are gridding with CVM only.
vv = isfield(Data, vars);

if sum(~vv) > 0
    pp = false(size(vars));
    if strcmpi('native', pixel_type)
        % If not doing the gridded data, then none of the weights fields will be present.
        pp = pp | ismember(vars, BEHR_publishing_gridded_fields.cvm_weight_vars);
    end
    if strcmpi('native', pixel_type) || Data(1).Only_CVM
        pp = pp | ismember(vars, BEHR_publishing_gridded_fields.psm_weight_vars);
    end

    missing_vars = vars(~vv & ~pp);
    if numel(missing_vars) > 0
        E.callError('missing_variable', 'The following variables are missing: %s', strjoin(missing_vars, ', '));
    end
end

vars = vars(vv);



% Remove variables that cannot be put into a CSV text file because multiple
% values are required per pixel
if strcmpi('txt', output_type)
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
gridded_vars = BEHR_publishing_gridded_fields.all_gridded_vars;

gg = ismember(vars,gridded_vars);
vars = vars(gg);

end

function vars = remove_vector_variables(vars)
E=JLLErrors;
if ~iscell(vars) || ~all(iscellcontents(vars,'ischar'))
    E.badinput('"vars" must be a cell array of variable names as strings')
end

% Define what variables cannot be saved in a text file and remove them
vector_vars = {'BEHRPressureLevels','BEHRScatteringWeights','BEHRAvgKernels','BEHRNO2apriori'};
vv = ismember(vars, vector_vars);
vars = vars(~vv);
end



function overwrite = make_hdf_file(Data_in, vars, attr, save_dir, savename, pixel_type, overwrite, current_git_heads, DEBUG_LEVEL)
E = JLLErrors;

[~,~,ext] = fileparts(savename);
if ~strcmp(ext,'.hdf')
    E.badinput('SAVENAME must have the extension .hdf, not %s', ext);
end
hdf_fullfilename = fullfile(save_dir, savename);

% Check if the file exists. We may be already set to automatically
% overwrite or not, otherwise we have to ask the user what to do.

if exist(hdf_fullfilename,'file')
    [overwrite, do_i_return] = do_i_overwrite(overwrite, hdf_fullfilename);
    if do_i_return
        return
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
            save_data(nans) = attr.(vars{v}).fillvalue;
        else
            E.notimplemented('Cell array variable');
        end
        
        if isa(save_data,'double')
            if ismember(vars{v}, BEHR_publishing_gridded_fields.flag_vars)
                save_data_uint = uint32(save_data);
                if ~isequal(save_data_uint, save_data)
                    E.callError('flag_conversion_error', 'After conversion to uint32, the %s field changed', vars{v});
                else
                    save_data = save_data_uint;
                end
            else
                % Convert doubles to singles to save space for the data people
                % will be downloading
                save_data = single(save_data);
            end
        end
        % Ensure that the fill value is of the same type as the data
        fill_val = cast(attr.(vars{v}).fillvalue, 'like', save_data);
        
        % Create the dataset, then write it and add the attributes
        h5create(hdf_fullfilename, var_name, sz, 'Datatype', class(save_data), 'FillValue', fill_val);
        h5write(hdf_fullfilename, var_name, save_data);
        
        atts = fieldnames(attr.(vars{v}));
        for a=1:numel(atts)
            if strcmp(atts{a},'fillvalue')
                continue % We've already handled the fill value with h5create
            end
            h5writeatt(hdf_fullfilename, var_name, atts{a}, attr.(vars{v}).(atts{a}));
        end
        
        % If doing gridded data, there's another attribute to set, which is
        % whether the value was gridded with PSM, CVM, neither, or is a
        % weight field.
        if strcmpi(pixel_type, 'gridded')
            grid_type_attr = 'gridding_method';
            if ~isempty(regexpi(vars{v}, 'weight'))
                grid_type = 'weight';
            elseif ismember(vars{v}, BEHR_publishing_gridded_fields.all_psm_vars) && ~Data_in(d).Only_CVM
                grid_type = 'parabolic spline method';
            elseif ismember(vars{v}, BEHR_publishing_gridded_fields.all_cvm_vars) || (ismember(vars{v}, BEHR_publishing_gridded_fields.all_psm_vars) && Data_in(d).Only_CVM)
                grid_type = 'constant value method';
            elseif ismember(vars{v}, BEHR_publishing_gridded_fields.flag_vars)
                grid_type = 'flag, bitwise OR';
            elseif ismember(vars{v}, BEHR_publishing_gridded_fields.publish_only_gridded_vars)
                grid_type = 'grid property';
            else
                grid_type = 'undefined';
            end
            h5writeatt(hdf_fullfilename, var_name, grid_type_attr, grid_type);
        end
        
        % If this is the BEHRQualityFlags field, add the bit meanings to
        % the attributes. We need to give behr_quality_flags input that it
        % can treat as if it were getting actual data in order to create
        % the flags definition cell array.
        if strcmpi(vars{v}, 'BEHRQualityFlags')
            [~,flags_definition] = behr_quality_flags();
            flags_definition = flags_definition(~iscellcontents(flags_definition, 'isempty'));
            h5writeatt(hdf_fullfilename, var_name, 'FlagMeanings', strjoin(flags_definition, ', '));
        end
    end
    
    % Write an attribute to the swath group describing if it is gridded or
    % native. Also include the version string and the Git head hashes
    switch lower(pixel_type)
        case 'native'
            swath_attr = 'OMI SP and BEHR data at native OMI resolution';
        case 'gridded'
            swath_attr = 'OMI SP and BEHR data gridded to 0.05 x 0.05 deg';
        otherwise
            E.badinput('"pixel_type" not recognized');
    end
    
    % General swath attributes
    h5writeatt(hdf_fullfilename, group_name, 'Description', swath_attr);
    h5writeatt(hdf_fullfilename, group_name, 'Version', BEHR_version());
    
    % Files read in in the process of generating the product
    swath_attr_fields = BEHR_publishing_gridded_fields.swath_attr_vars;
    for a=1:numel(swath_attr_fields)
        attr_val = Data_in(d).(swath_attr_fields{a});
        if iscellstr(attr_val)
            attr_val = strjoin(attr_val, ', ');
        elseif isnumeric(attr_val)
            if isscalar(attr_val)
                attr_val = num2str(attr_val);
            else
                attr_val = mat2str(attr_val);
            end
        elseif ~ischar(attr_val)
            E.notimplemented('Attribute values must be a number, string, or cell array of strings');
        end
        h5writeatt(hdf_fullfilename, group_name, swath_attr_fields{a}, attr_val);
    end

    h5writeatt(hdf_fullfilename, group_name, 'GitHead-Core_Pub', current_git_heads.core);
    h5writeatt(hdf_fullfilename, group_name, 'GitHead-BEHRUtils_Pub', current_git_heads.behr_utils);
    h5writeatt(hdf_fullfilename, group_name, 'GitHead-GenUtils_Pub', current_git_heads.gen_utils);

    
    if DEBUG_LEVEL > 2; toc; end
end
end



function ask_to_overwrite = make_txt_file(Data_in, vars, attr, date_string, save_dir, savename, ask_to_overwrite, DEBUG_LEVEL) %#ok<INUSD>

warning('make_txt_file has not been validated for BEHR > v2.1C, check the resulting files carefully before use')

if ~strcmp(savename(end),'_')
    savename = strcat(savename,'_');
end
txt_filename = strcat(savename, date_string, '.txt');
txt_fullfilename = fullfile(save_dir, txt_filename);

% Check if the file exists. Give the user 3 options if it does: abort,
% overwrite, overwrite all.

if exist(txt_fullfilename,'file')
    [ask_to_overwrite, do_i_return] = do_i_overwrite(ask_to_overwrite, txt_fullfilename);
    if do_i_return
        return;
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
                        v = attr.Loncorn.fillvalue;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Loncorn2'
                    v = Data_in(s).Loncorn(2,p);
                    if isnan(v)
                        v = attr.Loncorn.fillvalue;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Loncorn3'
                    v = Data_in(s).Loncorn(3,p);
                    if isnan(v)
                        v = attr.Loncorn.fillvalue;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Loncorn4'
                    v = Data_in(s).Loncorn(4,p);
                    if isnan(v)
                        v = attr.Loncorn.fillvalue;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Latcorn1'
                    v = Data_in(s).Latcorn(1,p);
                    if isnan(v)
                        v = attr.Latcorn.fillvalue;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Latcorn2'
                    v = Data_in(s).Latcorn(2,p);
                    if isnan(v)
                        v = attr.Latcorn.fillvalue;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Latcorn3'
                    v = Data_in(s).Latcorn(3,p);
                    if isnan(v)
                        v = attr.Latcorn.fillvalue;
                    end
                    fprintf(fid,format_spec{a},v);
                case 'Latcorn4'
                    v = Data_in(s).Latcorn(4,p);
                    if isnan(v)
                        v = attr.Latcorn.fillvalue;
                    end
                    fprintf(fid,format_spec{a},v);
                otherwise
                    v = Data_in(s).(header_cell{a})(p);
                    if isnan(v)
                        v = attr.(header_cell{a}).fillvalue;
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

function [overwrite, do_return] = do_i_overwrite(overwrite, full_file_name)
E = JLLErrors;

do_return = false;
if overwrite < 0
    user_ans = ask_multichoice(sprintf('File %s exists.\nOverwrite, skip, abort, overwrite all, or overwrite none?', full_file_name), {'o','s','a','oa','on'});
    user_ans = lower(user_ans);
    switch user_ans
        case 'o'
            delete(full_file_name);
        case 'oa'
            delete(full_file_name);
            overwrite = 1;
        case 's'
            do_return = true;
        case 'on'
            overwrite = 0;
            do_return = true;
        case 'a'
            E.userCancel;
        otherwise
            E.notimplemented('User answer %s', user_ans);
    end
elseif overwrite > 0
    delete(full_file_name);
else
    fprintf('%s exists, skipping\n', full_file_name);
    do_return = true;
end
end

function mycleanup()
err=lasterror; %#ok<LERR>
if ~isempty(err.message)
    fprintf('MATLAB exiting due to problem: %s\n', err.message);
    if ~isempty(gcp('nocreate'))
        delete(gcp)
    end 

    exit(1)
end
end
