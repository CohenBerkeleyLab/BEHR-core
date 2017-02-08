function [  ] = behr_staging_check( stage_path, stage_pattern, old_path, old_pattern )
%BEHR_staging_check Verifies that all expected files are present in a staging directory
%   BEHR_STAGING_CHECK( STAGE_PATH, STAGE_PATTERN, OLD_PATH, OLD_PATTERN )
%   Looks for files matching STAGE_PATTERN in STAGE_PATH and compares the
%   days present to files matching OLD_PATTERN in OLD_PATH. Will report on
%   duplicate days in each folder or days present in one folder but not the
%   other. Useful for verifying that all days are present in a staging
%   folder compared to a published folder before uploading.

E = JLLErrors;

if nargin < 4
    stage_path = '/Volumes/share-sat/SAT/BEHR/WEBSITE/staging/behr_gridded-hdf_v2-1Arev1';
    stage_pattern = 'OMI_BEHR_*';
    old_path = '/Volumes/share-sat/SAT/BEHR/WEBSITE/webData/behr_regridded_hdf';
    old_pattern = 'OMI_BEHR_*';
    if nargin ~= 0
        E.badinput('Must provide 4 inputs or none, the latter case will compare %s files to %s.', fullfile(stage_path, stage_pattern), fullfile(old_path, old_pattern));
    end
end

F_new = dir(fullfile(stage_path, stage_pattern));
F_old = dir(fullfile(old_path, old_pattern));

dnums_new = nan(size(F_new));
dnums_old = nan(size(F_old));

for a=1:numel(F_new)
    [s,e] = regexp(F_new(a).name,'\d\d\d\d\d\d\d\d');
    dnums_new(a) = datenum(F_new(a).name(s:e),'yyyymmdd');
end
for a=1:numel(F_old)
    [s,e] = regexp(F_old(a).name,'\d\d\d\d\d\d\d\d');
    dnums_old(a) = datenum(F_old(a).name(s:e),'yyyymmdd');
end

% First check that there are no duplicates within either set
fprintf('Checking for duplicates...\n');
nu = numel(unique(dnums_new));
n = numel(dnums_new);
if nu < n
    fprintf('There are %d duplicates in the new directory:\n', n-nu);
    for a=1:numel(dnums_new)
        if sum(dnums_new == dnums_new(a)) > 1
            fprintf('\t%s\n',F_new(a).name);
        end
    end
end

nu = numel(unique(dnums_old));
n = numel(dnums_old);
if nu < n
    fprintf('There are %d duplicates in the old directory:\n', n-nu);
    for a=1:numel(dnums_old)
        if sum(dnums_old == dnums_old(a)) > 1
            fprintf('\t%s\n',F_old(a).name);
        end
    end
end

% Second, check for which dates are present in one but not the other
xx = ismember(dnums_new, dnums_old);
if any(~xx)
    fprintf('Some days are present in the new directory but not the old:\')
    xxf = find(~xx);
    for a=1:numel(xxf)
        fprintf('\t%s\n',F_new(xxf(a)).name);
    end
end

xx = ismember(dnums_old, dnums_new);
if any(~xx)
    fprintf('Some days are present in the old directory but not the new:\n')
    xxf = find(~xx);
    for a=1:numel(xxf)
        fprintf('\t%s\n',F_old(xxf(a)).name);
    end
end

