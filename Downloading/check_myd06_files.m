function check_myd06_files
% CHECK_MYD06_FILES Checks which MODIS cloud files cannot be opened
%
% I've been having trouble with certain MODIS cloud files being corrupted and
% unopenable using MATLAB's hdfinfo command, so this function will go through
% all the MYD06 files, try to open them, and if that fails, write them out to a
% file in the home directory.
%
% Josh Laughner <joshlaugh5@gmail.com> 16 Aug 2015

moddir = '/mnt/sat/SAT/MODIS/MYD06_L2/';
startyear = 2011;
startday = 248;
firsttime = true;


fid = fopen('~/bad_myd06.txt','w');
cleanupObj = onCleanup(@() closeFile(fid));
D = dir(fullfile(moddir,'20*'));
for a=1:numel(D)
    yr = D(a).name;
    if str2double(yr) < startyear
        continue
    end
    fprintf('Now checking folder %s\n',yr);
    F = dir(fullfile(moddir,yr,'*.hdf'));
    for b=1:numel(F)
        if str2double(F(b).name(15:17)) < startday && firsttime
            continue
        end
        firsttime = false;
        fprintf('\tFile = %s\n',F(b).name);
        try
            hdfinfo(fullfile(moddir,yr,F(b).name));
        catch err
%            if strcmpi(err.identifier, 'MATLAB:imagesci:validate:fileOpen')
                fprintf('\t%s bad, saving\n',F(b).name);
                fprintf(fid, 'MYD06_L2/%s/%s\n', yr, F(b).name);
%            else
%                rethrow(err)
%            end
        end
    end
end
end

function closeFile(fid)
fclose(fid);
end
