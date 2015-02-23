% runscript_template.m - Template for a script to execute Matlab on the BRC
% cluster. Follow the comments below to edit this script.
%
% Josh Laughner <joshlaugh5@gmail.com> 14 Jan 2015

% We enclose the entire script in a try-catch loop so that if an error
% occurs, Matlab does not stop execution (which would potentially leave
% it awaiting interactive input on a compute node). Instead, we handle
% the error ourselve by printing the error message to the console and
% exiting with exit code > 0. If desired, the user can check the
% identifier on the error and edit the exit code to give more info.

try

% Most (or all) of my code meant to be run on the cluster should check the
% global variable onCluster to determine whether to execute code that needs
% to be run only on the cluster or not.  This includes calls to addpath()
% necessary (since I don't know if the cluster installation of Matlab
% allows additional folders to be added to the default path) and opening
% parallel pools.  
global onCluster;
onCluster = true;

% Note that I have prevented Matlab from automatically opening a parallel
% pool when parallel code is used (e.g. parfor). This way, by keeping calls
% to parpool() (which opens a parallel pool) inside if statements
% controlled by onCluster, when I run code on my local machine, no parallel
% pool will be opened, and the parfor loop will execute as a normal for
% loop.



% This will retrieve the variable MATLAB_NUM_THREADS from the bash shell
% that started this instance of Matlab.  By setting this variable (using
% the syntax: export MATLAB_NUM_THREADS=...) this script will automatically
% read in that environmental variable.  This is useful because it allows
% you to get the number of cores available from the SLURM scheduler (or
% whatever scheduler the cluster is using) via something like:
%   export MATLAB_NUM_THREADS=$SLURM_NTASKS_PER_NODE
% assuming you had previously set SLURM to use a specific number of tasks
% per node.

global numThreads;
tmp = getenv('MATLAB_NUM_THREADS');
if ~isempty(tmp); % getenv returns an empty string if a variable is not set in the shell
    numThreads = str2double(tmp);
    clear('tmp');
end

% As an aside, I'm guess that given the way Matlab is set up on the cluster
% (i.e. you run an actual Matlab instance rather than accessing Matlab
% workers on the cluster from your lab computer) that it might only be
% possible to assign one node per Matlab instance - that is, I don't know
% if Matlab will be able to communicate across nodes. Note that in R2014a,
% the limitation on the number of local workers was removed (was 12
% previously)

% Now execute your Matlab code here.  Keep in mind you may need to specify
% whole paths to custom function files, or add the paths.

% ==================== CODE TO RUN ==================== %



% ================== END CODE TO RUN ================== %


% Checks for an existing parallel pool and closes it, then exits Matlab
if ~isempty(gcp('nocreate'))
    delete(gcp);
end
exit;

catch err
    errmsg = getReport(err);
    fprintf('MATLAB exiting due to problem: \n\n%s',errmsg);
    
    % Check that there isn't an active parallel pool - if there is, delete it
    if ~isempty(gcp('nocreate'))
        delete(gcp)
    end

    % Exit with non-zero status code
    exit(1);
end
