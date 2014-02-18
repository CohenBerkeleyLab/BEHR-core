function satload(method)

% satload - Add file paths to Matlab for this session that are
% useful/needed for certain satellite work.  Takes a single argument
% (currently) that determines which set of files to load.
% Methods:
%   'states' = load the path to have the state outlines

if strcmp(method,'states')
   addpath('/Volumes/share/GROUP/SAT/BEHR/Ashley_tools/Ashley/Mapping/tools')
   disp('/share/GROUP/SAT/BEHR/Ashley_tools/Ashley/Mapping/tools added to path')
end