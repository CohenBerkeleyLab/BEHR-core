function [  ] = update_behr_version_website( )
%UPDATE_BEHR_VERSION_WEBSITE Writes the current BEHR version to a text file accessible through the website
%   Some of the web tools I've written depend on being able to read the
%   current BEHR version number to download the proper data. That's written
%   to a text file on the file server in a location accessible from the web
%   site. This function writes the current version defined in BEHR_version
%   to that text file. By running this automatically on the satellite
%   download computer, it ensures the version number is always correct.

verfile = fullfile(behr_paths.website_staging_dir,'..','webData','behr_version.txt');
fid = fopen(verfile,'w');
fprintf(fid, BEHR_version);
fclose(fid);

end

