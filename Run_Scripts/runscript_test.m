% runscript_template.m - Template for a script to execute Matlab on the BRC
% cluster. Follow the comments below to edit this script.
%
% Josh Laughner <joshlaugh5@gmail.com> 14 Jan 2015

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

% Set the paths for read_omno2 and BEHR_main
global sp_mat_dir;
sp_mat_dir = '/Volumes/share-sat/SAT/BEHR/SP_Files_2014';
global omi_he5_dir;
omi_he5_dir = '/Volumes/share/GROUP/SAT/OMI/OMNO2_32';
global modis_myd06_dir;
modis_myd06_dir = '/Volumes/share-sat/SAT/MODIS/MYD06_L2';
global modis_mcd43_dir;
modis_mcd43_dir = '/Volumes/share-sat/SAT/MODIS/MCD43C3';
global globe_dir;
globe_dir = '/Volumes/share-sat/SAT/BEHR/GLOBE_Database';

global behr_mat_dir;
behr_mat_dir = '/Volumes/share-sat/SAT/BEHR/BEHR_Files_2014';
global amf_tools_path;
amf_tools_path = '/Users/Josh/Documents/MATLAB/BEHR/AMF_tools';
global no2_profile_path;
no2_profile_path = '/Volumes/share-sat/SAT/BEHR/Monthly_NO2_Profiles';

%Run both
cd('~/Documents/MATLAB/BEHR/Read_Data');
read_omno2_v_aug2012('2013-08-01','2013-08-31');
cd('~/Documents/MATLAB/BEHR/BEHR_Main');
BEHR_main('2013-08-01','2013-08-31');

% ================== END CODE TO RUN ================== %


% Checks for an existing parallel pool and closes it, then exits Matlab
if ~isempty(gcp('nocreate'))
    delete(gcp);
end
exit;
