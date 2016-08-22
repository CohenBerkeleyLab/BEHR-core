function [  ] = update_behr_version_website( )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

verfile = fullfile(BEHR_paths('website_staging_dir'),'..','webData','behr_version.txt');
fid = fopen(verfile,'w');
fprintf(fid, BEHR_version);
fclose(fid);

end

