function satload(method)

% satload - Add file paths to Matlab for this session that are
% useful/needed for certain satellite work.  Takes a single argument
% (currently) that determines which set of files to load.
% Methods:
%   'states' = load the path to have the state outlines
%   'mod04' = load the modis MOD04_L2 and its subfolders

if strcmp(method,'states')
   addpath('/Volumes/share/GROUP/SAT/BEHR/Ashley_tools/Ashley/Mapping/tools');
   disp('/share/GROUP/SAT/BEHR/Ashley_tools/Ashley/Mapping/tools added to path')
elseif strcmpi(method,'mod04')
   addpath(genpath('/Volumes/share/GROUP/SAT/MODIS/MOD04_L2'));
   disp('MOD04_L2 on the server added to the path')
elseif strcmp(method,'globe')
    addpath('/Volumes/share/GROUP/SAT/BEHR/GLOBE_files')
    disp('GLOBE_files on the server added to the path')
end