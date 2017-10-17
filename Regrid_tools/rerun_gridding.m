function [  ] = rerun_gridding( start_date, end_date, prof_mode, varargin )
%RERUN_GRIDDING Runs psm_wrapper on a series of BEHR files, saving the new gridded output
%   RERUN_GRIDDING( START_DATE, END_DATE, PROF_MODE ) will load each BEHR
%   file generated using daily or monthly profiles (specified by the string
%   PROF_MODE, which must match the profile mode in the files to be loaded)
%   in the date range given from the standard BEHR directory, regrid the
%   Data struct into the OMI struct, and save the results to the current
%   directory. START_DATE and END_DATE may be strings understood as dates
%   by Matlab or date numbers.
%
%   This behavior can be modified with the following parameters:
%       'behr_mat_dir' - set what directory the BEHR .mat files should be
%       loaded from. By default, the path given as behr_paths.behr_mat_dir.
%
%       'save_dir' - the directory the new files should be saved to. By
%       default, the current directory. This directory must already exist.
%
%       'region' - which region's files to regrid. Default is 'us'.
%
%       'overwrite' - default false, whether to overwrite existing files in
%       the save directory. Must be a scalar logical (true or false).
%
%       'DEBUG_LEVEL' - run time verbosity, also passed to psm_wrapper to
%       control its verbosity.

E = JLLErrors;
p = inputParser;
p.addParameter('behr_mat_dir', behr_paths.behr_mat_dir);
p.addParameter('save_dir', '.');
p.addParameter('region', 'us');
p.addParameter('overwrite', false);
p.addParameter('DEBUG_LEVEL', 2);

p.parse(varargin{:});
pout = p.Results;

load_dir = pout.behr_mat_dir;
save_dir = pout.save_dir;
region = pout.region;
overwrite = pout.overwrite;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

start_date = validate_date(start_date);
end_date = validate_date(end_date);

if ~ischar(load_dir)
    E.badinput('The value for "behr_mat_dir" must be a string');
elseif ~exist(load_dir, 'dir')
    E.badinput('The directory given for "behr_mat_dir" (%s) does not exist', save_dir);
end

if ~ischar(save_dir)
    E.badinput('The value for "save_dir" must be a string');
elseif ~exist(save_dir, 'dir')
    E.badinput('The directory given for "save_dir" (%s) does not exist', save_dir);
end

if ~isscalar(overwrite) || ~islogical(overwrite)
    E.badinput('The value for "overwrite" must be a scalar logical');
end

if ~isscalar(DEBUG_LEVEL) || ~isnumeric(DEBUG_LEVEL)
    E.badinput('The value for "DEBUG_LEVEL" must be a scalar number');
end

for this_dnum=start_date:end_date
    behr_file = behr_filename(this_dnum, prof_mode, region);
    save_name = fullfile(save_dir, behr_file);
    if ~overwrite && exist(save_name, 'file')
        if DEBUG_LEVEL > 0
            fprintf('Output file %s exists, skipping\n', save_name);
        end
        continue
    else
        if DEBUG_LEVEL > 0
            fprintf('Now doing gridding for %s\n', datestr(this_dnum));
        end
    end
    
    D = load(fullfile(load_dir, behr_file), 'Data');
    Data = D.Data;
    OMI = psm_wrapper(Data, Data(1).Grid, DEBUG_LEVEL); %#ok<NASGU>
    
    save(save_name, 'Data', 'OMI');
end

end

