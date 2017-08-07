function [  ] = BEHR_initial_setup( )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[behr_repo_path, myname] = fileparts(mfilename('fullpath'));

make_paths_file();

    function make_paths_file
        % Check that behr_paths.m does not already exist anywhere. If there's one
        % in the right place, then ask if we want to overwrite it. If there's one
        % somewhere else, the user will need to move or delete it so that we don't
        % get two behr_paths.m files.
        paths_filename = 'behr_paths.m';
        template_filename = 'behr_paths_template.m';
        paths_file = fullfile(behr_repo_path,'Utils','Constants',paths_filename);
        template_file = fullfile(behr_repo_path,'Utils','Constants',template_filename);
        
        if exist(paths_file,'file')
            user_ans = input(sprintf('A %s file already exists in BEHR/Utils/Constants. Replace it? (y to replace, any other key to abort): ', paths_filename), 's');
            if ~(strcmpi(user_ans,'Yes') || strcmpi(user_ans, 'y'))
                error('BEHR_setup:user_cancel','A %s file already exists and you have chosen not to overwrite it.', paths_filename);
            end
        elseif ~isempty(which(paths_filename))
            error('BEHR_setup:file_exists','%s exists on your Matlab search path, but at %s, not at BEHR/Utils/Constants. Please delete or move that version to BEHR/Utils/Constants.',paths_filename, which(paths_filename));
        end
        
        % This defines all the paths that should exist in that file. The
        % field name is the name the property will be given, the 'comment'
        % subfield is the explanatory comment, and the 'default' subfield
        % is what the default will be. If either field is omitted, that is
        % fine. A generic comment will be inserted and a blank path will be
        % inserted. Including the subfield "no_quote" at all will cause it
        % not to automatically quote the default path, which is useful if
        % you need to allow for multiple paths in a cell array.
        
        sat_file_server = '128.32.208.13';
        wrf_file_server = 'cohenwrfnas.dyn.berkeley.edu';
        
        if ismac
            sat_folder = fullfile('/Volumes','share-sat');
            wrf1_folder = fullfile('/Volumes', 'share-wrf1');
            wrf2_folder = fullfile('/Volumes', 'share-wrf2');
        elseif isunix
            sat_folder = fullfile('/mnt','share-sat');
            wrf1_folder = fullfile('/mnt', 'share-wrf1');
            wrf2_folder = fullfile('/mnt', 'share-wrf2');
        elseif ispc
            drive_letter_request = 'Enter the drive letter (just the letter, not the ":\" that "%s" is mounted as';
            sat_folder = strcat(input(sprintf(drive_letter_request, 'share-sat'), 's'), ':');
            wrf1_folder = strcat(input(sprintf(drive_letter_request, 'share-wrf1'), 's'), ':');
            wrf2_folder = strcat(input(sprintf(drive_letter_request, 'share-wrf2'), 's'), ':');
        end
            
        
        % Local repos/folders
        paths.classes.comment = sprintf('The Matlab classes repository; it should contain at least the class JLLErrors. May be cloned from share-sat/SAT/BEHR/BEHR_MatlabClasses_GitRepo.git on the file server at %s.', sat_file_server);
        paths.classes.default = fullfile(behr_repo_path, '..', 'Classes');
        paths.utils.comment = sprintf('The directory of the general Matlab Utils repository (not the one inside the BEHR repo). May be cloned from share-sat/SAT/BEHR/MiscUtils.git on the file server at %s.', sat_file_server);
        paths.utils.default = fullfile(behr_repo_path, '..', 'Utils');
        paths.amf_tools_dir.comment = 'The AMF_tools directory in the BEHR repository on your computer. It should contain the files damf.txt and nmcTmpYr.txt';
        paths.amf_tools_dir.default = fullfile(behr_repo_path, 'AMF_tools');
        paths.psm_dir.comment = 'The PSM_Gridding directory in the BEHR repository on your computer. It should contain the files PSM_Main.py and psm_wrapper.m';
        paths.psm_dir.default = fullfile(behr_repo_path, 'PSM_Gridding');
        
        % Matlab file folders
        paths.sp_mat_dir.comment = sprintf('The mounted path to the directory on the file server at %s containing OMI_SP_vX-YZ_yyyymmdd.mat files. The file server should be mounted on your computer.',sat_file_server);
        paths.sp_mat_dir.default = fullfile(sat_folder, 'SAT', 'BEHR', 'SP_Files_2014');
        paths.behr_mat_dir.comment = sprintf('The mounted path to the directory on the file server at %s containing OMI_BEHR_vX-YZ_yyyymmdd.mat files. The file server should be mounted on your computer.',sat_file_server);
        paths.behr_mat_dir.default = fullfile(sat_folder, 'SAT', 'BEHR', 'BEHR_Files_2014');
        
        % OMI and ancillary data folders
        paths.omno2_dir.comment = sprintf('The OMNO2 directory on the file server at %s. It should contain folders organized by year and month with OMI-Aura_L2-OMNO2 files in them. The file server should be mounted on your computer.',sat_file_server);
        paths.omno2_dir.default = fullfile(sat_folder, 'SAT', 'OMI', 'OMNO2', 'version_3_3_0');
        paths.ompixcor_dir.comment = sprintf('The OMPIXCOR directory on the file server at %s. It should contain folders organized by year and month with OMI-Aura_L2-OMNPIXCOR files in them. The file server should be mounted on your computer.',sat_file_server);
        paths.ompixcor_dir.default = fullfile(sat_folder, 'SAT', 'OMI', 'OMPIXCOR', 'version_003');
        paths.myd06_dir.comment = sprintf('Please find the MYD06_L2 directory on the file server at %s. It should contain folders for each year with MYD06_L2 files in them. The file server should be mounted on your computer.',sat_file_server);
        paths.myd06_dir.default = fullfile(sat_folder, 'SAT', 'MODIS', 'MYD06_L2');
        paths.mcd43c1_dir.comment = sprintf('The MCD43C1 directory on the file server at %s. It should contain folders for each year with MCD43C1 files in them. The file server should be mounted on your computer.',sat_file_server);
        paths.mcd43c1_dir.default = fullfile(sat_folder, 'SAT', 'MODIS', 'MCD43C1');
        paths.globe_dir.comment = sprintf('The GLOBE database directory on the file server at %s. It should contain files a10g through p10g and their .hdr files. The file server should be mounted on your computer.',sat_file_server);
        paths.globe_dir.default = fullfile(sat_folder, 'SAT', 'BEHR','GLOBE_Database');

        paths.website_staging_dir.comment = sprintf('The directory where data should be staged before being put on the folders visible to the website. Also on the file server at %s which should be mounted on your computer.', sat_file_server);
        paths.website_staging_dir.default = fullfile(sat_folder,'SAT','BEHR','WEBSITE','staging');
        
        % WRF data is spread across multiple volumes locally. It's just too
        % big.
        paths.wrf_monthly_profiles.comment = sprintf('The path that contains the WRF_BEHR*.nc monthly profile files. This should be on the file server at %s.', wrf_file_server);
        paths.wrf_monthly_profiles.default = fullfile(wrf1_folder,'MonthlyProfiles');
        
        paths.wrf_profiles.comment = sprintf('Add all paths that contain WRF profiles. These should be folders that are organized by year and month, with wrfout_d01 files in them. These will be found on the file server at %s, all volumes must be mounted on your computer.', wrf_file_server);
        paths.wrf_profiles.default = sprintf('{''%s'', ''%s''}', fullfile(wrf1_folder, 'Outputs'), fullfile(wrf2_folder, 'Outputs'));
        paths.wrf_profiles.no_quote = true;
        
        %%%%%%%%%%%%%%%%%%
        % Write the file %
        %%%%%%%%%%%%%%%%%%
        fid_template = fopen(template_file, 'r');
        fid_new = fopen(paths_file, 'w');
        
        tline = fgetl(fid_template);
        [~,paths_classname] = fileparts(paths_filename);
        [~,template_classname] = fileparts(template_filename);
        while ischar(tline)
            tline = strrep(tline, template_classname, paths_classname);
            if strcmp(strtrim(tline), '%#PATHS')
                write_paths(paths, fid_new);
            else
                fprintf(fid_new, '%s\n', tline);
            end
            tline = fgetl(fid_template);
        end
        
        fclose(fid_template);
        fclose(fid_new);
        
        fprintf('\n!!! Defaults %s file created at %s. Review it and edit the paths as needed. !!!\n\n', paths_filename, paths_file);
    end

    
end

function str_out = wrap_comment(comment, indent_level)
% Assume each line is 76 characters long, allow 4 spaces for tabs and 2
% characters for the "% " at the beginning
nchars = 76 - 4*indent_level - 2;
tabs = repmat('\t', 1, indent_level);
i = 1;
str_out = '';
while i < length(comment)
    % Find the last space before the end of the line
    if length(comment) - i > nchars
    i2 = min(length(comment), i+nchars);
    j = strfind(comment(i:i2), ' ');
    j = j(end);
    else
        j = length(comment) - i + 1;
    end
    str_out = [str_out, tabs, '%% ', strtrim(comment(i:i+j-1)), '\n'];
    %str_out{end+1} = sprintf('%% %s\n',strtrim(comment(i:i+j-1)));
    i = i+j;
end
end

function write_paths(paths, fid)
fns = fieldnames(paths);
for a=1:numel(fns)
    if isfield(paths.(fns{a}), 'comment');
        this_comment = wrap_comment(paths.(fns{a}).comment, 2);
    else
        this_comment = wrap_comment(fns{a}, 2);
    end
    fprintf(fid, this_comment);
    
    if isfield(paths.(fns{a}), 'default')
        this_default = paths.(fns{a}).default;
    else
        this_default = '';
    end
    
    if isfield(paths.(fns{a}), 'no_quote')
        fprintf(fid, '\t\t%s = %s;\n\n', fns{a}, this_default);
    else
        fprintf(fid, '\t\t%s = ''%s'';\n\n', fns{a}, this_default);
    end
end
end
