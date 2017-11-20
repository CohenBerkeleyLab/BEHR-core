function [  ] = integration_test( varargin )
%INTEGRATION_TEST Quick check that all the BEHR code executes
%   INTEGRATION_TEST() runs read_main, BEHR_main, and
%   BEHR_publishing_v2 for one day. Data is output to ./UnitTestData/tmp.
%   The only purpose of this function is a rapid verification that the
%   primary components of BEHR work; it does not check for correctness.
%
%   INTEGRATION_TEST( 'free-branch' ) overrides the requirement that BEHR-core and
%   BEHR-core-utils be on the same branch.
%
%   INTEGRATION_TEST( 'reading' ) only does the reading code.
%
%   INTEGRATION_TEST( 'main' ) only does the main code.
%
%   INTEGRATION_TEST( 'pub' ) only does the publishing code.
%
%   Any of these three can be combined to do two of the three pieces: e.g.
%   INTEGRATION_TEST( 'reading', 'main' ) will run the reading and main
%   code, but not the publishing code.
E = JLLErrors;

xx = strcmpi('free-branch', varargin);
if any(xx)
    req_same_branch = false;
else
    req_same_branch = true;
end

do_read = ismember('reading', varargin);
do_main = ismember('main', varargin);
do_pub = ismember('pub', varargin);

if ~do_read && ~do_main && ~do_pub
    % If no step specified, assume we do all of them
    do_read = true;
    do_main = true;
    do_pub = true;
end

if req_same_branch
    core_branch = git_branch(behr_paths.behr_core);
    utils_branch = git_branch(behr_paths.behr_utils);
    if ~strcmp(core_branch, utils_branch)
        E.callError('git', 'Cannot proceed; core and core-utils directories on different branches');
    end
end

mydir = fileparts(mfilename('fullpath'));
out_dir = fullfile(mydir, 'UnitTestData', 'tmp');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
pub_native_dir = fullfile(out_dir, 'native');
if ~exist(pub_native_dir, 'dir')
    mkdir(pub_native_dir);
end
pub_gridded_dir = fullfile(out_dir, 'gridded');
if ~exist(pub_gridded_dir, 'dir')
    mkdir(pub_gridded_dir);
end

test_date = '2012-02-01';

if do_read
    timer_read = tic;
    read_main('start', test_date, 'end', test_date, 'sp_mat_dir', out_dir, 'overwrite', true);
    t_read = toc(timer_read);
else
    t_read = NaN;
end
if do_main
    timer_main = tic;
    BEHR_main('start', test_date, 'end', test_date, 'sp_mat_dir', out_dir, 'behr_mat_dir', out_dir, 'profile_mode', 'monthly', 'overwrite', true);
    BEHR_main('start', test_date, 'end', test_date, 'sp_mat_dir', out_dir, 'behr_mat_dir', out_dir, 'profile_mode', 'daily', 'overwrite', true);
    t_main = toc(timer_main);
else 
    t_main = NaN;
end
if do_pub
    timer_pub = tic;
    BEHR_publishing_main('start', test_date, 'end', test_date, 'output_type', 'hdf', 'pixel_type', 'native', 'mat_dir', out_dir, 'save_dir', pub_native_dir,...
        'organize', false, 'overwrite', true);
    BEHR_publishing_main('start', test_date, 'end', test_date, 'output_type', 'hdf', 'pixel_type', 'gridded', 'mat_dir', out_dir, 'save_dir', pub_gridded_dir,...
        'organize', false, 'overwrite', true);
    t_pub = toc(timer_pub);
else
    t_pub = NaN;
end

fprintf('All major BEHR components finished successfully.\n');
fprintf('Times:\n\t Reading: %f s \n\t Main: %f s \n\t Publishing: %f s \n', t_read, t_main, t_pub);

end
