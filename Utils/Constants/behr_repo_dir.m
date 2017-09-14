function [ bdir ] = behr_repo_dir(  )
%BEHR_REPO_DIR Returns the top directory of the BEHR repository
%   BDIR = BEHR_REPO_DIR() Returns the top directory of the BEHR
%   repository; should be OS independent. This can be used to generate
%   paths relative to the repo root that can then be valid across any
%   machine, at least within the repository.
%
%   It does rely on being within the Utils/Constants folder within the BEHR
%   repository; a warning will be issued if this does not appear to be the
%   case. You should check the path it returns in that case.

warning('BEHR_REPO_DIR has been deprecated in favor of behr_paths.behr_core');

mydir = fileparts(mfilename('fullpath'));
%Assumes this is in the Utils/Constants directory within the BEHR repo
if isempty(regexp(mydir,'Utils/Constants$','once'))
    warning('BEHR_REPO_DIR assumes that it is within Utils/Constants in the BEHR repo, this does not appear to be true (path = %s)',mydir)
end
bdir = fullfile(mydir,'..','..');

end

