function [  ] = unit_test_driver( self_test )
%UNIT_TEST_DRIVER Driver function for BEHR unit test
%   This function, when called, asks a series of questions interactively to
%   determine how the unit tests should proceed. It is capable of
%   automatically generating OMI_SP and OMI_BEHR files using the current
%   versions of read_main and BEHR_main, if this is requested,
%   it saves the resulting files in a subdirectory of "UnitTestData" which
%   will be created in the same directory as this function. The
%   subdirectory will be named "ProducedYYYYMMDD". It will also contain a
%   text file that describes the status of the BEHR git repository at the
%   time of production, including both the commit hash and one line
%   description of HEAD and the diff against HEAD.
%
%   Whether you produce the data with this function or not, it will then
%   called both BEHR_UNIT_TEST and (if testing read_main)
%   READING_PRIORI_TESTS. BEHR_UNIT_TEST takes a pair of Data or OMI
%   structures and attempts to verify that they are the same. If they are
%   not the same, the unit test will fail overall for that date, but will
%   also print out information about why it failed. Since changes to the
%   algorithm intended to change the output will cause it to fail, it is up
%   to the user to determine if the changes are the expected ones.
%   READING_PRIORI_TESTS does not test the a priori profiles; rather it
%   checks that certain elements of the OMI_SP files make sense a priori,
%   that is, on their own without needing a previous version of the file to
%   compare against.
%
%   At a minimum, this should be used before producing a new version of
%   BEHR to test selected dates (see below) against the existing version.
%   This would be done by allowing it to generate the new data and compare
%   against the files in the directories specified by BEHR_paths().
%
%   This does not run the entire OMI data record. Within the code is a cell
%   array of test dates which will be run. A few days are included in this
%   array to test normal operation, the rest are there to test behavior
%   under weird cases that have caused issues in the past. Consequently,
%   you should add days to this as you find days that cause the BEHR
%   algorithm to error or behave strangely, but you should not remove
%   existing days.
%
%   Internally, the ranges of dates that have daily profiles are specified,
%   so as more daily profiles become available, one should modify this
%   function to update that. It will only generate and test daily profile
%   files within those specified ranges.
%
%   UNIT_TEST_DRIVER( true ) runs a self test, so it only tries to do one
%   day. This is useful if you've made changes to the unit test code itself
%   and just want to make sure it works.
%
%   Josh Laughner <joshlaugh5@gmail.com> 8 May 2017

E = JLLErrors;
DEBUG_LEVEL = 2;

if ~exist('self_test','var')
    self_test = false;
end

% Test these dates. It's a good idea to check at least one regular day
% before the row anomaly started (2005-2006), after it was at its worst
% (after July 2011), plus zoom mode operation in both time periods.
% Additional days should be added that have caused or illuminated bugs

test_region = 'US';

test_dates = {'2005-06-02';... % pre-row anomaly summertime day, at least one day after zoom mode finishes
              '2006-01-01';... % pre-row anomaly wintertime day, at least one day after zoom mode finishes
              '2012-06-03';... % post-row anomaly summertime day, at least one day after zoom mode finishes
              '2013-01-01';... % post-row anomaly wintertime day, at least one day after zoom mode finishes
              '2014-07-08';... % post-row anomaly zoom mode, found by looking for days where OMPIXCORZ is produced for BEHR-relevant orbits for US region
              '2006-09-20';... % pre-row anomaly zoom mode, found by looking for days where OMPIXCORZ is produced for BEHR-relevant orbits for US region
              '2005-07-13';... % day mentioned in the 2.1Arev1 changelog with no NO2 data
              '2010-01-29';... % day mentioned in the 2.1Arev1 changelog with no NO2 data
              '2005-05-04';... % the center lon/lat for the OMPIXCOR product just different enough that cutting down by bounds results in different size arrays, so switched to find_submatrix2
              '2005-05-14';... % Has both a row that is only partially fill values in lon/lat and the OMPIXCOR corners are mostly 0
              '2010-05-09';... % day where the first swath is only 4 long in the along track dimension, which previously caused an error in the gridding algorithm, since it matches the length of the corner dimension
              '2006-04-26';... % day where along and across track dimensions in the first orbit are both 60, which previously screwed up the dimension checking in the gridding algorithm
              '2009-02-15'...  % day where the along track dimension is cut all the way down to one, which exposed a bug in the matlab-python interface
              };

% These are dates that the algorithm should be run for, but for which it is
% okay if no data is produced. This allows the unit tests to skip them
test_dates_no_data = {'2016-05-30';... % OMI was in safe mode; algorithm should gracefully handle the lack of data
                      '2007-12-19'};   % Should read in an empty structure from read_main.m, BEHR_main should just skip it and not try to produce anything
              
% These are the date ranges for which daily profiles are available. When
% running in daily mode, only test dates within these ranges will be
% produced. The first column represents the beginning of each range; the
% second column the end.
daily_profiles_available = {'2005-01-01', '2005-12-31';...
                            '2007-01-01', '2007-12-31';...
                            '2012-01-01', '2013-12-31'};
                  
test_dates = unique(cat(1,test_dates,test_dates_no_data));

if self_test
    test_dates = test_dates(1);
end

my_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(my_dir,'SubTests'));
what_to_test = ask_multichoice('Which step should be tested?', {'all', 'reading', 'behrmain', 'publishing'});

use_behrpaths = ask_yn('Use the paths specified by BEHR_paths() for the old data?');
generate_new_data = ask_yn('Generate the new files? If not you will be asked to choose the directories to load new files from');
if generate_new_data
    root_save_folder = make_data_folder();
    read_save_folder = fullfile(root_save_folder, 'Reading');
    main_root_save_folder = fullfile(root_save_folder, 'Main');
    pub_root_save_folder = fullfile(root_save_folder, 'Publishing');
end

save_results_to_file = ask_yn('Save results to file? (If not, will be printed to screen).');
if save_results_to_file
    results_file = make_results_file(what_to_test);
    fid = fopen(results_file,'w');
else 
    % An fid of 1 will cause fprintf to print to the command window, as if
    % no fid was given
    fid = 1;
end

try
    prompt_str = ['\nSpecify any fields to ignore in unit testing, separated by a space.\n',...
                  'Regular expressions can be used. By default, fields beginning with\n',...
                  '"GitHead" are ignored because they are expected to be different if the\n',...
                  'algorithm has changed. To override this, make one of the strings input\n',...
                  '"keepgit" (without the quotes), i.e. entering "keepgit .*File" will\n',...
                  'ignore any field ending in "File" but do compare the GitHead fields: '];
    fields_to_ignore = input(prompt_str, 's');
    fields_to_ignore = strsplit(fields_to_ignore);

    xx = strcmpi('keepgit', fields_to_ignore);
    if ~any(xx)
        fields_to_ignore = veccat({'GitHead.*'},fields_to_ignore);
    else
        fields_to_ignore(xx) = [];
    end

    if generate_new_data
        make_git_report(behr_paths.behr_core, 'GitReport-Core.txt');
        make_git_report(behr_paths.behr_utils, 'GitReport-BEHRUtils.txt');
        make_git_report(behr_paths.utils, 'GitReport-GenUtils.txt');
    end
    switch what_to_test
        % Each of the testing subfunctions allows paths to be given to them 
        % by test_all() to minimized user interaction if all three steps are
        % to be run. I've set it up so that if empty strings are passed, it
        % considers those paths to not be given, but something has to be passed.
        case 'reading'
            success = test_reading('', '');
        case 'behrmain'
            success_m = test_behr_main('monthly', '', '');
            success_d = test_behr_main('daily', '', '');
            success = success_m & success_d;
        case 'publishing'
            success_m = test_publishing('monthly', '', '', '', '');
            success_d = test_publishing('daily', '', '', '', '');
            success = success_m * success_d;
        case 'all'
            success = test_all();
        otherwise
            E.notimplemented(what_to_test);
    end
    
    for a=1:numel(success)
        fprintf(fid, '%s: %s\n', datestr(test_dates{a}), passfail(success(a)));
    end
    fprintf(fid, 'Overall: %s\n', passfail(all(success)));
    
    msg = sprintf('BEHR unit test completed on %s step(s): %s', what_to_test, datestr(now));
    border = repmat('*', 1, numel(msg));
    fprintf(fid, '\n%s\n', border);
    fprintf(fid, '%s\n', msg);
    fprintf(fid, '%s\n\n', border);
catch err
    if fid > 2
        fclose(fid);
    end
    rethrow(err);
end
if fid > 2
    fclose(fid);
end

if save_results_to_file
    fprintf('Results saved to %s\n', results_file);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NESTED FUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dfolder = make_data_folder()
        dfolder = fullfile(my_dir, 'UnitTestData', sprintf('Produced%s', datestr(today, 'yyyymmdd')));
        if exist(dfolder, 'dir')
            if ~ask_yn(sprintf('Directory\n %s\n exists, it will be cleared before continuing. Proceed?', dfolder))
                E.userCancel()
            else
                % This is set up to delete all files in one of the
                % "ProducedYYYYMMDD" folders b/c all files in those folders
                % should be produced by the code in the state represented
                % by the GitReport.txt file. If we allow ourselves to only
                % remove a subset of files, that is no longer guaranteed to
                % be true.
                
                remove_contents(dfolder);
            end
        else
            mkdir(dfolder);
        end
        fprintf('Unit test data will be stored in %s\n', dfolder);
    end

    function rfile = make_results_file(test_steps)
        rfolder = fullfile(my_dir, 'UnitTestResults');
        if ~exist(rfolder, 'dir')
            mkdir(rfolder);
        end
       rfilename = sprintf('BEHR_%s_Unit_Test_Results_%s.txt', test_steps, datestr(now, 'yyyymmdd_HHMMSS'));
       rfile = fullfile(rfolder, rfilename);
    end

    function make_git_report(repo_dir, report_name)
        currdir = cd(repo_dir);
        try
            % Overall status (current branch, modified/added/deleted/untracked
            % files). Use --porcelain to remove formatting (bold, color, etc)
            % and --branch to force it to show the branch. --no-pager means it
            % won't try to put the output through "less" and avoids a "terminal
            % not fully functional" warning.
            [~, gitstat] = system('git --no-pager status --porcelain --branch');
            
            % Get the one line commit message for the last, decorated with any
            % tags or branch heads. Remove color to avoid introducing special
            % characters into the text.
            [~, githead] = system('git --no-pager log -1 --pretty=oneline --decorate --no-color');
            
            % Get the differenced since the last commit, sans color and pager
            % for the same reasons as above. By specifically diffing against
            % HEAD, we get staged and unstaged changes.
            [~, gitdiff] = system('git --no-pager diff --no-color HEAD');
        catch err
            cd(currdir)
            rethrow(err);
        end
        cd(currdir);
        
        % Extract the branch from the status - with "--porcelain --branch"
        % it is on its own line prefaced by ##
        [i,j] = regexp(gitstat, '##.*?\n', 'once');
        gitbranch = gitstat(i:j);
        gitbranch = strrep(gitbranch, '#', '');
        gitbranch = strtrim(gitbranch);
        gitstat(i:j) = [];
        
        % Also add space before each file in the diff (usually the "diff
        % --git" is bolded so it's easier to see but we removed formatting)
        gitdiff = strrep(gitdiff, 'diff --git', sprintf('\ndiff --git'));
        
        gfid = fopen(fullfile(root_save_folder, report_name), 'w');
        begin_msg = sprintf('Git report on %s for unit test data generated on %s', repo_dir, datestr(now));
        gborder = repmat('*', size(begin_msg));
        
        fprintf(gfid, '%s\n%s\n%s\n\n', gborder, begin_msg, gborder); 
        fprintf(gfid, 'Current branch: %s\n\n', gitbranch);
        fprintf(gfid, 'HEAD at time of generation:\n%s\n%s\n\n', githead, gborder);
        fprintf(gfid, 'Git status at time of generation\n  (M = modified, A = added, D = deleted, ?? = untracked):\n\n%s\n%s\n\n', gitstat, gborder);
        fprintf(gfid, 'Git diff (working dir against HEAD) at time of generation:\n\n%s', gitdiff);
        fclose(gfid);
    end

    function successes = test_all()
        % If doing all three tests, give the option to set up all the directories now
        if ~use_behrpaths
            if ask_yn('Compare against an old unit test?')
                old_root = getdir('Choose the root ProducedYYYYMMDD directory', {});
                old_sp_dir = fullfile(old_root, 'Reading');
                old_daily_behr_dir = fullfile(old_root, 'Main', 'daily');
                old_monthly_behr_dir = fullfile(old_root, 'Main', 'monthly');
                old_daily_native_dir = fullfile(old_root, 'Publishing', 'native_hdf', 'daily');
                old_monthly_native_dir = fullfile(old_root, 'Publishing', 'native_hdf', 'monthly');
                old_daily_gridded_dir = fullfile(old_root, 'Publishing', 'gridded_hdf', 'daily');
                old_monthly_gridded_dir = fullfile(old_root, 'Publishing', 'gridded_hdf', 'monthly');
            else
                old_sp_dir = getdir('You''ll need to choose the directory with the old OMI_SP files', test_dates);
                old_monthly_behr_dir = getdir('You''ll need to choose the directory with the old monthly OMI_BEHR files', test_dates);
                old_daily_behr_dir = getdir('You''ll need to choose the directory with the old daily OMI_BEHR files', test_dates);
                old_daily_native_dir = getdir('You''ll need to choose the directory with the old daily native pixel HDF files', test_dates);
                old_monthly_native_dir = getdir('You''ll need to choose the directory with the old monthly native pixel HDF files', test_dates);
                old_daily_gridded_dir = getdir('You''ll need to choose the directory with the old daily gridded HDF files', test_dates);
                old_monthly_gridded_dir = getdir('You''ll need to choose the directory with the old monthly gridded HDF files', test_dates);
            end
        else
            % if using behr_paths, these will be set automatically
            old_sp_dir = '';
            old_monthly_behr_dir = '';
            old_daily_behr_dir = '';
            old_daily_native_dir = '';
            old_monthly_native_dir = '';
            old_daily_gridded_dir = '';
            old_monthly_gridded_dir = '';
        end

        if ~generate_new_data
            new_sp_dir = getdir('You''ll need to choose the directory with the new OMI_SP files', test_dates);
            new_behr_dir = getdir('You''ll need to choose the directory with the new OMI_BEHR files with "daily" and "monthly" subfolders', test_dates);
            new_native_dir = getdir('You''ll need to choose the directory containing the new native HDF files with "daily" and "monthly" subfolders', test_dates);
            new_gridded_dir = getdir('You''ll need to choose the directory containing the new gridded HDF files with "daily" and "monthly" subfolders', test_dates);
        else
            % if generating new data, these are set from the save_folder automatically
            new_sp_dir = '';
            new_behr_dir = '';
            new_native_dir = '';
            new_gridded_dir = '';
        end
        read_success = test_reading(old_sp_dir, new_sp_dir);
        behr_monthly_success = test_behr_main('monthly',old_monthly_behr_dir, new_behr_dir, read_save_folder);
        behr_daily_success = test_behr_main('daily',old_daily_behr_dir, new_behr_dir, read_save_folder);
        pub_monthly_success = test_publishing('monthly',old_monthly_native_dir, old_monthly_gridded_dir, new_native_dir, new_gridded_dir, fullfile(main_root_save_folder,'monthly'));
        pub_daily_success = test_publishing('daily',old_daily_native_dir, old_daily_gridded_dir, new_native_dir, new_gridded_dir, fullfile(main_root_save_folder,'daily'));
        successes = read_success & behr_monthly_success & behr_daily_success & pub_monthly_success & pub_daily_success;
    end

    function successes = test_reading(old_dir, new_dir)
        if generate_new_data
            % If generating new data, then our new_dir will always be the location where we generate the new data.
            new_dir = read_save_folder;
            mkdir(new_dir);
            for i=1:numel(test_dates)
                read_main('start', test_dates{i}, 'end', test_dates{i}, 'sp_mat_dir', new_dir, 'overwrite', true, 'region', test_region);
            end
        else
            % Otherwise, a directory for new files may have already been passed (if running all tests, generally). 
            % if not, ask now which directory contains the new SP files.
            if isempty(new_dir)
                new_dir = getdir('You''ll need to choose the directory with the new OMI_SP files', test_dates);
            end

            % This fixed some weird bug where "read_save_folder" wasn't set because we weren't generating data, but
            % it got used later. That probably shouldn't happen, so that bug should be fixed eventually and this
            % removed.
            if ~exist('save_folder', 'var')
                read_save_folder = new_dir;
            end
        end
        
        % Since a common use case of this function is to test a new version against the prior version, we gave
        % the user the option at the beginning of using the standard paths for old data. If that's not what they
        % chose, or an old directory wasn't already given, we need to ask now.
        if use_behrpaths
            old_dir = behr_paths.SPMatSubdir(test_region);
        elseif isempty(old_dir)
            old_dir = getdir('You''ll need to choose the directory with the old OMI_SP files', test_dates);
        end
        
        successes = true(size(test_dates));
        for i=1:numel(test_dates)
            if DEBUG_LEVEL > 0
                fprintf(fid, '\n');
            end
            filepat = sp_savename(test_dates{i}, test_region, '.mat', true);
            try
                [old_data, old_file] = load_by_glob(fullfile(old_dir, filepat));
                [new_data, new_file] = load_by_glob(fullfile(new_dir, filepat));
            catch err
                if strcmp(err.identifier, 'load_by_glob:file_not_found')
                    if ismember(test_dates{i}, test_dates_no_data)
                        if DEBUG_LEVEL > 0
                            fprintf(fid, 'No data for %s as expected\n', test_dates{i});
                        end
                    else
                        if DEBUG_LEVEL > 0
                            fprintf(fid, 'FAIL: No data produced for %s!!!\n', test_dates{i});
                        end
                        successes(i) = false;
                    end
                    continue
                else
                    rethrow(err);
                end
            end
            if DEBUG_LEVEL > 0
                fprintf(fid, '\nChecking %s\n', test_dates{i});
                fprintf(fid, 'Loaded old file: %s\n', old_file{1});
                fprintf(fid, 'Loaded new file: %s\n', new_file{1});
            end
            
            if DEBUG_LEVEL > 0
                header_msg = '***** Running priori tests on data read in ****';
                header_border = repmat('*', 1, length(header_msg));
                fprintf(fid, '\n%1$s\n%2$s\n%1$s\n', header_border, header_msg);
            end
            
            successes(i) = reading_priori_tests(new_data.Data, DEBUG_LEVEL, fid) && successes(i);
            
            if DEBUG_LEVEL > 0
                header_msg = '***** Running reading unit tests, comparing to previous data ****';
                header_border = repmat('*', 1, length(header_msg));
                fprintf(fid, '\n%1$s\n%2$s\n%1$s\n', header_border, header_msg);
            end
            
            successes(i) = behr_unit_test(new_data.Data, old_data.Data, DEBUG_LEVEL, fid, fields_to_ignore) && successes(i);
        end
    end

    
    %%%%%%%%%%%%%%%%%%%
    % BEHR MAIN TESTS %
    %%%%%%%%%%%%%%%%%%%
    
    function successes = test_behr_main(prof_mode, old_dir, new_dir, sp_data_dir)
        if generate_new_data
            new_dir = fullfile(main_root_save_folder, lower(prof_mode));
            mkdir(new_dir)
            % If sp_data_dir not already given, give the choice of using behr_paths.sp_mat_dir or a user-specified dir
            if ~exist('sp_data_dir', 'var')
                if ask_yn('Use the paths specified by behr_paths for the SP files to be read into BEHR_main?')
                    sp_data_dir = behr_paths.SPMatSubdir(test_region);
                else
                    sp_data_dir = getdir('You''ll need to choose the directory with existing OMI_SP files', test_dates);
                end
            end
            
            for i=1:numel(test_dates)
                if strcmpi(prof_mode, 'daily') && ~can_do_daily(test_dates{i})
                    continue
                end
                
                BEHR_main('start', test_dates{i}, 'end', test_dates{i}, 'behr_mat_dir', new_dir, 'sp_mat_dir', sp_data_dir, 'profile_mode', prof_mode, 'overwrite', true);
            end
        else
            % If we're not generating data, then check if new_dir is not empty (i.e. already given)
            if isempty(new_dir)
                new_dir = getdir(sprintf('You''ll need to choose the directory with the new %s OMI_BEHR files', prof_mode), test_dates);
            end
        end
        
        % Since a common use case of this function is to test a new version against the prior version, we gave
        % the user the option at the beginning of using the standard paths for old data. If that's not what they
        % chose, or an old directory wasn't already given, we need to ask now.
        if use_behrpaths
                old_dir = behr_paths.BEHRMatSubdir(test_region, prof_mode);
        elseif isempty(old_dir)
            old_dir = getdir(sprintf('You''ll need to choose the directory with the old %s OMI_BEHR files', prof_mode), test_dates);
        end
        
        successes_data = true(size(test_dates));
        successes_grid = true(size(test_dates));
        for i=1:numel(test_dates)
            if strcmpi(prof_mode, 'daily') && ~can_do_daily(test_dates{i})
                successes_data(i) = true;
                successes_grid(i) = true;
                continue
            end
            
            if DEBUG_LEVEL > 0
                fprintf(fid, '\n');
            end
            filepat = behr_filename(test_dates{i}, prof_mode, test_region, '.mat', true);
            try
                [old_data, old_file] = load_by_glob(fullfile(old_dir, filepat));
                [new_data, new_file] = load_by_glob(fullfile(new_dir, filepat));
            catch err
                if strcmp(err.identifier, 'load_by_glob:file_not_found')
                    if ismember(test_dates{i}, test_dates_no_data)
                        if DEBUG_LEVEL > 0
                            fprintf(fid, 'No data for %s as expected\n', test_dates{i});
                        end
                    else
                        if DEBUG_LEVEL > 0
                            fprintf(fid, 'FAIL: No data produced for %s!!!\n', test_dates{i});
                        end
                        successes_data(i) = false;
                        successes_grid(i) = false;
                    end
                    continue
                else
                    rethrow(err);
                end
            end
            if DEBUG_LEVEL > 0
                fprintf(fid, '\nChecking %s\n', test_dates{i});
                fprintf(fid, 'Loaded old file: %s\n', old_file{1});
                fprintf(fid, 'Loaded new file: %s\n', new_file{1});
            end
            
            % Only the native pixel struct has all of the necessary fields
            % to assess the accuracy of the quality flags. We must assume
            % that the gridding algorithm does its job properly... at least
            % until we write a unit test for that.
            if DEBUG_LEVEL > 0
                header_msg = '***** Running priori tests on result of main algorithm, Data struct ****';
                header_border = repmat('*', 1, length(header_msg));
                fprintf(fid, '\n%1$s\n%2$s\n%1$s\n', header_border, header_msg);
            end
            
            successes_data(i) = main_priori_tests(new_data.Data, DEBUG_LEVEL, fid) && successes_data(i);
            
            if DEBUG_LEVEL > 0
                header_msg = '***** Running BEHR_main unit tests on Data struct ****';
                header_border = repmat('*', 1, length(header_msg));
                fprintf(fid, '\n%1$s\n%2$s\n%1$s\n', header_border, header_msg);
            end
            
            successes_data(i) = behr_unit_test(new_data.Data, old_data.Data, DEBUG_LEVEL, fid, fields_to_ignore) && successes_data(i);
            
            if DEBUG_LEVEL > 0
                header_msg = '***** Running BEHR_main unit tests on OMI struct ****';
                header_border = repmat('*', 1, length(header_msg));
                fprintf(fid, '\n%1$s\n%2$s\n%1$s\n', header_border, header_msg);
            end
            
            successes_grid(i) = behr_unit_test(new_data.OMI, old_data.OMI, DEBUG_LEVEL, fid, fields_to_ignore) && successes_grid(i);
        end
        
        successes = successes_data & successes_grid;
    end


    %%%%%%%%%%%%%%%%%%%%
    % PUBLISHING TESTS %
    %%%%%%%%%%%%%%%%%%%%
    
    function successes = test_publishing(prof_mode, old_native_dir, old_gridded_dir, new_native_dir, new_gridded_dir, behr_data_dir)
        if generate_new_data
            % If we're generating new data, then our new file directories will always be in the save folder
            % (where we generate the new data)
            new_native_dir = fullfile(pub_root_save_folder, 'native_hdf', lower(prof_mode));
            new_gridded_dir = fullfile(pub_root_save_folder, 'gridded_hdf', lower(prof_mode));
            if ~exist('behr_data_dir', 'var')
                if ask_yn('Use the paths specified by behr_paths for the BEHR files to be read into BEHR_publishing_main?')
                    behr_data_dir = behr_paths.BEHRMatSubdir(test_region, prof_mode);
                else
                    behr_data_dir = getdir(sprintf('You''ll need to choose the directory with existing %s OMI_BEHR files', prof_mode), test_dates);
                end
            end
            
            if ~exist(new_native_dir, 'dir')
                mkdir(new_native_dir);
            end
            if ~exist(new_gridded_dir, 'dir')
                mkdir(new_gridded_dir)
            end
            
            for i=1:numel(test_dates)
                if strcmpi(prof_mode, 'daily') && ~can_do_daily(test_dates{i})
                    continue
                end
                
                BEHR_publishing_main('start', test_dates{i}, 'end', test_dates{i}, 'output_type', 'hdf', 'pixel_type', 'native', 'mat_dir', behr_data_dir, 'save_dir', new_native_dir,...
                    'organize', false, 'overwrite', true);
                BEHR_publishing_main('start', test_dates{i}, 'end', test_dates{i}, 'output_type', 'hdf', 'pixel_type', 'gridded', 'mat_dir', behr_data_dir, 'save_dir', new_gridded_dir,...
                    'organize', false, 'overwrite', true);
            end
        else
            % Otherwise, this may be already given
            if isempty(new_native_dir)
                new_native_dir = getdir(sprintf('You''ll need to choose the directory containing the new %s native HDF files', prof_mode), test_dates);
            end
            if isempty(new_gridded_dir)
                new_gridded_dir = getdir(sprintf('You''ll need to choose the directory containing the new %s gridded HDF files', prof_mode), test_dates);
            end
        end
        
        % Since a common use case of this function is to test a new version against the prior version, we gave
        % the user the option at the beginning of using the standard paths for old data. If that's not what they
        % chose, or an old directory wasn't already given, we need to ask now.
        if use_behrpaths
            old_root_dir = behr_paths.website_staging_dir;
            % Assume that the staging directory is adjacent to the webData
            % directory where the files are actually moved to show up on
            % the website
            old_native_dir = fullfile(old_root_dir, '..', 'webData', sprintf('behr_%s_hdf', lower(prof_mode)));
            old_gridded_dir = fullfile(old_root_dir, '..', 'webData', sprintf('behr_%s_regridded_hdf', lower(prof_mode)));
        else
            if isempty(old_native_dir)
                old_native_dir = getdir(sprintf('You''ll need to choose the directory with the old %s native pixel HDF files', prof_mode), test_dates);
            end
            if isempty(old_gridded_dir)
                old_gridded_dir = getdir(sprintf('You''ll need to choose the directory with the old %s gridded HDF files', prof_mode), test_dates);
            end
        end
        
        successes_native = test_publishing_subfunc(prof_mode, old_native_dir, new_native_dir);
        successes_grid = test_publishing_subfunc(prof_mode, old_gridded_dir, new_gridded_dir);
        
        successes = successes_native & successes_grid;
    end

    function successes = test_publishing_subfunc(prof_mode, old_dir, new_dir)
        native_or_gridded = regexp(new_dir, '(native)|(gridded)', 'match', 'once');
        
        successes = false(size(test_dates));
        for i=1:numel(test_dates)
            if strcmpi(prof_mode, 'daily') && ~can_do_daily(test_dates{i})
                continue
            end
            
            if DEBUG_LEVEL > 0
                fprintf(fid, '\n');
            end
            filepat = behr_filename(test_dates{i}, prof_mode, test_region, '.hdf', true);
            try
                [old_data, old_file] = load_hdf_by_glob(fullfile(old_dir, filepat));
                [new_data, new_file] = load_hdf_by_glob(fullfile(new_dir, filepat));
            catch err
                if strcmp(err.identifier, 'load_hdf_by_glob:file_not_found')
                    if ismember(test_dates{i}, test_dates_no_data)
                        if DEBUG_LEVEL > 0
                            fprintf(fid, 'No data for %s as expected\n', test_dates{i});
                        end
                        successes(i) = true;
                    else
                        if DEBUG_LEVEL > 0
                            fprintf(fid, 'FAIL: No data produced for %s!!!\n', test_dates{i});
                        end
                    end
                    continue
                else
                    rethrow(err);
                end
            end
            if DEBUG_LEVEL > 0
                fprintf(fid, '\nChecking %s\n', test_dates{i});
                fprintf(fid, 'Loaded old file: %s\n', old_file);
                fprintf(fid, 'Loaded new file: %s\n', new_file);
            end
            
            if DEBUG_LEVEL > 0
                header_msg = sprintf('***** Running BEHR_publishing unit tests on %s HDFs ****', native_or_gridded);
                header_border = repmat('*', 1, length(header_msg));
                fprintf(fid, '\n%1$s\n%2$s\n%1$s\n', header_border, header_msg);
            end
            
            successes(i) = behr_unit_test(new_data, old_data, DEBUG_LEVEL, fid, fields_to_ignore);
        end
    end

    % NESTED UTILITY FUNCTIONS %
    function b = can_do_daily(date_in)
        dnums = cellfun(@datenum, daily_profiles_available);
        date_check = datenum(date_in) >= dnums(:,1) & datenum(date_in) <= dnums(:,2);
        b = any(date_check);
    end
end


%%%%%%%%%%%%%%%%
% SUBFUNCTIONS %
%%%%%%%%%%%%%%%%

function s = passfail(b)
if b
    s = 'PASS';
else
    s = 'FAIL';
end
end

function d = getdir(prompt, test_dates)
if ~isempty(test_dates)
    fprintf('%s for the following dates:\n  %s\n (press ENTER)\n', prompt, strjoin(test_dates, ', '));
else
    fprintf('%s\n (press ENTER)\n', prompt);
end
input('','s'); %wait for the user

E=JLLErrors;
if isDisplay
    d = uigetdir;
else
    while true
        d = input('Enter the directory: ', 's');
        if strcmpi(d,'q')
            E.userCancel;
        elseif exist(d, 'dir')
            break
        else
            fprintf('That directory does not exist.\n')
        end
    end
end
end

function [Data, file_name] = load_hdf_by_glob(file_pattern)
F = dir(file_pattern);
if numel(F) < 1
    error('load_hdf_by_glob:file_not_found', 'No file matching "%s" found', file_pattern);
elseif numel(F) > 1
    error('load_hdf_by_glob:too_many_files', 'Multiple files matching "%s" found', file_pattern);
end

hdf_dir = fileparts(file_pattern);
file_name = fullfile(hdf_dir, F(1).name);
hdfi = h5info(file_name);
Data = behrhdf2struct(hdfi);


end

function remove_contents(directory)
F = dir(fullfile(directory, '*'));
for i=1:numel(F)
    if ~strcmp(F(i).name, '.') && ~strcmp(F(i).name, '..')
        if F(i).isdir
            remove_contents(fullfile(directory, F(i).name));
            rmdir(fullfile(directory, F(i).name));
        else
            delete(fullfile(directory, F(i).name));
        end
    end
end
end


