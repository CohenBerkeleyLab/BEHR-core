function [ indiv_stats, overall_stats  ] = behr_prod_test( varargin  )
%[ INDIV_STATS, OVERALL_STATS] = BEHR_PROD_TEST() 
%   Tests a sample of OMI_BEHR files for differences. Whenever making a new
%   version of BEHR, it's good to do some basic checking to make sure that
%   the differences are what you expect. This function will choose a random
%   selection of BEHR files to compare and report on individual, day-by-day
%   stats and the overall differences. This version uses the hard coded
%   directories and file patterns.
%
%[ INDIV_STATS, OVERALL_STATS] = BEHR_PROD_TEST( NEW_DIR, NEW_PATTERN, OLD_DIR, OLD_PATTERN )
%   This version looks for new files matching the glob pattern NEW_PATTERN
%   in NEW_DIR and old files in OLD_DIR matching OLD_PATTERN.
%
% Additional parameters:
%
%   'nfiles' - how many files to load and test. Default is 100.
%
%   'checkvar' - which variable in .mat files to test. Must be the string
%   'Data' or 'OMI'. Default is 'Data', but has no effect if not loading
%   .mat files.
%
%   'fields' - a cell array of string indicating which fields in the files
%   to check. Default is {'BEHRColumnAmountNO2Trop', 'BEHRAMFTrop',
%   'BEHRColumnAmountNO2TropVisOnly', 'BEHRAMFTropVisOnly'}.
%
%   'start', 'end' - start and end dates of the period to draw the data
%   from. Must be either a date number or a date string that Matlab
%   recognizes automatically. Default is 2005-01-01 to today.

% Enter the directory and file name pattern for the new and old versions.
% The pattern must be a string that, when used in dir(), uniquely returns
% only BEHR files of the desired version.

%%%% USER OPTIONS %%%%
p = inputParser;
p.addOptional('new_dir', '.', @ischar);
p.addOptional('new_pattern', 'OMI_BEHR*', @ischar);
p.addOptional('old_dir', behr_paths.behr_mat_dir, @ischar);
p.addOptional('old_pattern', 'OMI_BEHR*', @ischar);

p.addParameter('nfiles', 100);
p.addParameter('checkvar', 'Data');
p.addParameter('fields', {'BEHRColumnAmountNO2Trop','BEHRAMFTrop','BEHRColumnAmountNO2TropVisOnly','BEHRAMFTropVisOnly'});
p.addParameter('start', '2005-01-01');
p.addParameter('end', today);

p.parse(varargin{:});
pout = p.Results;

new_dir = pout.new_dir;
new_pattern = pout.new_pattern;
old_dir = pout.old_dir;
old_pattern = pout.old_pattern;

n_files = pout.nfiles;
checkvar = pout.checkvar;
fields_to_check = pout.fields;
start_date = pout.start;
end_date = pout.end;

% Validation
if ~exist(new_dir, 'dir')
    E.badinput('new_dir "%s" does not exist', new_dir);
elseif ~exist(old_dir, 'dir')
    E.badinput('old_dir "%s" does not exist', old_dir);
end

if ~ischar(new_pattern)
    E.badinput('NEW_PATTERN must be a string')
elseif ~ischar(old_pattern)
    E.badinput('OLD_PATTERN must be a string')
end

if ~isnumeric(n_files) || ~isscalar(n_files) || n_files < 1 || mod(n_files, 1) ~= 0
    E.badinput('The value for "nfiles" must be a scalar, positive, whole number')
end

if ~ismember(checkvar, {'Data', 'OMI'})
    E.badinput('The value for "checkvar" must be the string "Data" or "OMI"');
end

if ~iscellstr(fields_to_check)
    E.badinput('The value for "fields" must be a cell array of strings');
end

validate_date(start_date);
validate_date(end_date);
    

%%%% END USER OPTIONS %%%%


F_new = dir(fullfile(new_dir,new_pattern));
F_new = cut_down_by_date(F_new, start_date, end_date);
F_old = dir(fullfile(old_dir,old_pattern));
F_old = cut_down_by_date(F_old, start_date, end_date);
n_files = min([n_files, numel(F_new), numel(F_old)]);

[~,~,fileext_new] = fileparts(F_new(1).name);
[~,~,fileext_old] = fileparts(F_old(1).name);

dnums_new = get_file_datenums(F_new);
dnums_old = get_file_datenums(F_old);

diff_struct = struct('num_dif_vals',0,'mean_difference',0,'mean_absolute_difference',0,'mean_percent_difference',0,'mean_absolute_percent_difference',0,...
    'median_difference',0,'median_absolute_difference',0,'median_percent_difference',0,'median_absolute_percent_difference',0,'differences',[],...
    'percent_differences',[],'value_pairs',[]);
mat_hdf_comp_bool = xor(strcmpi(fileext_new,'.mat'), strcmpi(fileext_old, '.mat'));
if mat_hdf_comp_bool;
    % If comparing a .mat and a .hdf file, we will do a different
    % comparison of fill values because the .hdf will have fill values
    % where the .mat has NaNs in some cases
    fills_struct = struct('num_new_nans_or_fills', 0, 'values_that_became_nans_or_fills', [],...
        'num_old_nans_or_fills', 0, 'values_that_replaced_nans_or_fills', []);
else
    fills_struct = struct('num_new_nans', 0, 'values_that_became_nans', [], 'lon_for_became_nans', [], 'lat_for_became_nans', [],...
        'num_new_fills', 0, 'values_that_became_fills', [], 'lon_for_became_fills', [], 'lat_for_became_fills', [],...
        'num_old_nans', 0, 'values_that_replaced_nans', [], 'lon_for_replaced_nans', [],  'lat_for_replaced_nans', [],...
        'num_old_fills', [], 'values_that_replaced_fills', [], 'lon_for_replaced_fills', [], 'lat_for_replaced_fills', []);
end
%substruct = struct('date','','num_dif_vals',0,'num_new_nans',0,'values_that_became_nans',[],'num_new_fills',0,'values_that_became_fills',[],...
%    'num_old_nans',0,'values_that_replaced_nans',[],'num_old_fills',0,'values_that_replaced_fills',[],'differences',[],'percent_differences',[],...
%    'value_pairs',[],'Longitude',[],'Latitude',[]);
substruct = struct('date', '', 'Longitude', [], 'Latitude', [], 'difference_stats', diff_struct, 'fill_and_nan_changes', fills_struct);
overall_stats = make_empty_struct_from_cell(fields_to_check, rmfield(substruct,'date'));
indiv_stats = make_empty_struct_from_cell(fields_to_check, substruct);
indiv_stats = repmat(indiv_stats, n_files, 1);

if isDisplay
    wb = waitbar(0, sprintf('Sampling %d files',n_files));
end

n = 0;
safety = 0;
while n < n_files
    safety = safety+1;
    
    % Choose a random new file, make sure there is a corresponding old one,
    % if so, load both and check the required fields.
    r = ceil(rand * numel(F_new));
    rold = dnums_old == dnums_new(r);
    if sum(rold) < 1
        continue
    end
    
    n = n+1;
    if isDisplay
        waitbar(n/n_files);
    end
    
    % Determine the file type from the extension, that will determine how
    % it is loaded and if fill values are checked.
    fill_vals = nan(numel(fields_to_check)+2, 2);
    [D_new.Data, fill_vals(:,1)] = load_data(fullfile(new_dir, F_new(r).name), [fields_to_check, {'Longitude', 'Latitude'}], checkvar);
    [D_old.Data, fill_vals(:,2)] = load_data(fullfile(old_dir, F_old(rold).name), [fields_to_check, {'Longitude', 'Latitude'}], checkvar);
    
    lon = cat_sat_data(D_new.Data,'Longitude');
    lat = cat_sat_data(D_new.Data,'Latitude');
    
    
    
    for a = 1:numel(fields_to_check)
        data_new = cat_sat_data(D_new.Data,fields_to_check{a}, 'vector', true);
        data_old = cat_sat_data(D_old.Data,fields_to_check{a}, 'vector', true);
        
        num_neq = sum(data_new(:) ~= data_old(:));
        
        if mat_hdf_comp_bool
            is_fill_or_nan_new = isnan(data_new(:)) | data_new(:) == fill_vals(a,1);
            is_fill_or_nan_old = isnan(data_old(:)) | data_old(:) == fill_vals(a,2);
            
            xx_new_nans_fills = is_fill_or_nan_new & ~is_fill_or_nan_old;
            num_new_nans_fills = sum(xx_new_nans_fills);
            values_now_nans_fills = data_old(xx_new_nans_fills);
            
            xx_old_nans_fills = ~is_fill_or_nan_new & is_fill_or_nan_old;
            num_old_nans_fills = sum(xx_old_nans_fills);
            values_replaced_nans_fills = data_new(xx_old_nans_fills);
            
            xx_good = ~(is_fill_or_nan_new | is_fill_or_nan_old);
        else
            xx_newnans = isnan(data_new(:)) & ~isnan(data_old(:));
            num_new_nans = sum(xx_newnans);
            values_now_nans = data_old(xx_newnans);
            lon_for_now_nans = lon(xx_newnans);
            lat_for_now_nans = lat(xx_newnans);
            if isnan(fill_vals(a,1))
                % Cannot test for fill value of NaN using == b/c nan == nan
                % returns false.
                is_new_fill = isnan(data_new(:));
                num_new_fills = num_new_nans;
                values_now_fills = values_now_nans;
                lon_for_now_fills = lon_for_now_nans;
                lat_for_now_fills = lat_for_now_nans;
            else
                is_new_fill = data_new(:) == fill_vals(a,1);
                xx_newfills = data_new(:) == fill_vals(a,1) & data_old(:) ~= fill_vals(a,2);
                num_new_fills = sum(xx_newfills);
                values_now_fills = data_old(xx_newfills);
                lon_for_now_fills = lon(xx_newfills);
                lat_for_now_nans = lat(xx_newfills);
            end
            
            xx_oldnans = ~isnan(data_new(:)) & isnan(data_old(:));
            num_old_nans = sum(xx_oldnans);
            values_replaced_nans = data_new(xx_oldnans);
            lon_for_rep_nans = lon(xx_oldnans);
            lat_for_rep_nans = lat(xx_oldnans);
            if isnan(fill_vals(a,2))
                is_old_fill = isnan(data_old(:));
                num_old_fills = num_old_nans;
                values_replaced_fills = values_replaced_nans;
                lon_for_rep_fills = lon_for_rep_nans;
                lat_for_rep_fills = lat_for_rep_nans;
            else
                is_old_fill = data_old(:) == fill_vals(a,2);
                xx_oldfills = data_new(:) ~= fill_vals(a,1) & data_old(:) == fill_vals(a,2);
                num_old_fills = sum(xx_oldfills);
                values_replaced_fills = data_new(xx_oldfills);
                lon_for_rep_fills = lon(xx_oldfills);
                lat_for_rep_fills = lat(xx_oldfills);
            end
            
            xx_good = ~(is_new_fill | is_old_fill);
        end
        
        
        del = data_new(xx_good) - data_old(xx_good);
        perdel = reldiff(data_new(xx_good), data_old(xx_good))*100;
        mean_diff = nanmean(del);
        mean_absdiff = nanmean(abs(del));
        mean_perdiff = nanmean(perdel);
        mean_absperdiff = nanmean(abs(perdel));
        median_diff = nanmedian(del);
        median_absdiff = nanmedian(abs(del));
        median_perdiff = nanmedian(perdel);
        median_abs_perdiff = nanmedian(abs(perdel));
        
        indiv_stats(n).(fields_to_check{a}).date = datestr(dnums_new(r),'yyyy-mm-dd');
        indiv_stats(n).(fields_to_check{a}).Longitude = lon(xx_good);
        indiv_stats(n).(fields_to_check{a}).Latitude = lat(xx_good);
        
        indiv_stats(n).(fields_to_check{a}).difference_stats.num_dif_vals = num_neq;
        indiv_stats(n).(fields_to_check{a}).difference_stats.mean_difference = mean_diff;
        indiv_stats(n).(fields_to_check{a}).difference_stats.mean_absolute_difference = mean_absdiff;
        indiv_stats(n).(fields_to_check{a}).difference_stats.mean_percent_difference = mean_perdiff;
        indiv_stats(n).(fields_to_check{a}).difference_stats.mean_absolute_percent_difference = mean_absperdiff;
        indiv_stats(n).(fields_to_check{a}).difference_stats.median_difference = median_diff;
        indiv_stats(n).(fields_to_check{a}).difference_stats.median_absolute_difference = median_absdiff;
        indiv_stats(n).(fields_to_check{a}).difference_stats.median_percent_difference = median_perdiff;
        indiv_stats(n).(fields_to_check{a}).difference_stats.median_absolute_percent_difference = median_abs_perdiff;
        indiv_stats(n).(fields_to_check{a}).difference_stats.differences = del;
        indiv_stats(n).(fields_to_check{a}).difference_stats.percent_differences = perdel;
        indiv_stats(n).(fields_to_check{a}).difference_stats.value_pairs = [data_new(xx_good), data_old(xx_good)];
        
        if mat_hdf_comp_bool
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.num_new_nans_or_fills = num_new_nans_fills;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.values_that_became_nans_or_fills = values_now_nans_fills;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.num_old_nans_or_fills = num_old_nans_fills;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.values_that_replaced_nans_or_fills = values_replaced_nans_fills;
        else
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.num_new_nans = num_new_nans;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.values_that_became_nans = values_now_nans;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.lon_for_became_nans = lon_for_now_nans;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.lat_for_became_nans = lat_for_now_nans;
            
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.num_new_fills = num_new_fills;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.values_that_became_fills = values_now_fills;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.lon_for_became_fills = lon_for_now_fills;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.lat_for_became_fills = lat_for_now_fills;
            
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.num_old_nans = num_old_nans;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.values_that_replaced_nans = values_replaced_nans;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.lon_for_replaced_nans = lon_for_rep_nans;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.lat_for_replaced_nans = lat_for_rep_nans;
            
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.num_old_fills = num_old_fills;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.values_that_replaced_fills = values_replaced_fills;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.lon_for_replaced_fills = lon_for_rep_fills;
            indiv_stats(n).(fields_to_check{a}).fill_and_nan_changes.lat_for_replaced_fills = lat_for_rep_fills;
        end
        
        overall_stats.(fields_to_check{a}).Longitude = cat(1, overall_stats.(fields_to_check{a}).Longitude, lon(xx_good));
        overall_stats.(fields_to_check{a}).Latitude = cat(1, overall_stats.(fields_to_check{a}).Latitude, lat(xx_good));
        
        overall_stats.(fields_to_check{a}).difference_stats.num_dif_vals = overall_stats.(fields_to_check{a}).difference_stats.num_dif_vals + num_neq;
        overall_stats.(fields_to_check{a}).difference_stats.differences = cat(1, overall_stats.(fields_to_check{a}).difference_stats.differences, del);
        overall_stats.(fields_to_check{a}).difference_stats.percent_differences = cat(1, overall_stats.(fields_to_check{a}).difference_stats.percent_differences, perdel);
        overall_stats.(fields_to_check{a}).difference_stats.value_pairs = cat(1, overall_stats.(fields_to_check{a}).difference_stats.value_pairs, [data_new(xx_good), data_old(xx_good)]);
        
        if mat_hdf_comp_bool
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.num_new_nans_or_fills = overall_stats.(fields_to_check{a}).fill_and_nan_changes.num_new_nans_or_fills + num_new_nans_fills;
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.values_that_became_nans_or_fills = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.values_that_became_nans_or_fills, values_now_nans_fills);
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.num_old_nans_or_fills = overall_stats.(fields_to_check{a}).fill_and_nan_changes.num_old_nans_or_fills + num_old_nans_fills;
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.values_that_replaced_nans_or_fills = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.values_that_replaced_nans_or_fills, values_replaced_nans_fills);
        else
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.num_new_nans = overall_stats.(fields_to_check{a}).fill_and_nan_changes.num_new_nans + num_new_nans;
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.values_that_became_nans = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.values_that_became_nans, values_now_nans);
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.lon_for_became_nans = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.lon_for_became_nans, lon_for_now_nans);
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.lat_for_became_nans = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.lat_for_became_nans, lat_for_now_nans);
            
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.num_new_fills = overall_stats.(fields_to_check{a}).fill_and_nan_changes.num_new_fills + num_new_fills;
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.values_that_became_fills = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.values_that_became_fills, values_now_fills);
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.lon_for_became_fills = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.lon_for_became_fills, lon_for_now_fills);
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.lat_for_became_fills = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.lat_for_became_fills, lat_for_now_fills);
            
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.num_old_nans = overall_stats.(fields_to_check{a}).fill_and_nan_changes.num_old_nans + num_old_nans;
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.values_that_replaced_nans = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.values_that_replaced_nans, values_replaced_nans);
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.lon_for_replaced_nans = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.lon_for_replaced_nans, lon_for_rep_nans);
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.lat_for_replaced_nans = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.lat_for_replaced_nans, lat_for_rep_nans);
            
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.num_old_fills = overall_stats.(fields_to_check{a}).fill_and_nan_changes.num_old_fills + num_old_fills;
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.values_that_replaced_fills =  cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.values_that_replaced_fills, values_replaced_fills);
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.lon_for_replaced_fills = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.lon_for_replaced_fills, lon_for_rep_fills);
            overall_stats.(fields_to_check{a}).fill_and_nan_changes.lat_for_replaced_fills = cat(1, overall_stats.(fields_to_check{a}).fill_and_nan_changes.lat_for_replaced_fills, lat_for_rep_fills);
        end
    end
    
    if safety > 10*n_files;
        warning('Loop has executed more than 10x the number of requested files, exiting via safety condition');
        break
    end
    
    F_new(r) = [];
    F_old(rold) = [];
end

for a=1:numel(fields_to_check)
    % Overall difference states computed at the end
    all_del = overall_stats.(fields_to_check{a}).difference_stats.differences;
    all_perdel = overall_stats.(fields_to_check{a}).difference_stats.percent_differences;
    overall_stats.(fields_to_check{a}).difference_stats.mean_difference = nanmean(all_del);
    overall_stats.(fields_to_check{a}).difference_stats.mean_absolute_difference = nanmean(abs(all_del));
    overall_stats.(fields_to_check{a}).difference_stats.mean_percent_difference = nanmean(all_perdel);
    overall_stats.(fields_to_check{a}).difference_stats.mean_absolute_percent_difference = nanmean(abs(all_perdel));
    overall_stats.(fields_to_check{a}).difference_stats.median_difference = nanmedian(all_del);
    overall_stats.(fields_to_check{a}).difference_stats.median_absolute_difference = nanmedian(abs(all_del));
    overall_stats.(fields_to_check{a}).difference_stats.median_percent_difference = nanmedian(all_perdel);
    overall_stats.(fields_to_check{a}).difference_stats.median_absolute_percent_difference = nanmedian(abs(all_perdel));
end

if isDisplay
    close(wb)
end

% Put variables in base workspace if no outputs to the function, and plot
% overall histograms
if nargout < 1
    putvar(indiv_stats,overall_stats);
end
for a=1:numel(fields_to_check)
    figure; hist(overall_stats.(fields_to_check{a}).difference_stats.differences, 50);
    title(sprintf('Differences in %s',fields_to_check{a}));
    
    figure; hist(overall_stats.(fields_to_check{a}).difference_stats.percent_differences, 50);
    title(sprintf('Percent differences in %s',fields_to_check{a}));
end

end

function dnums = get_file_datenums(F)
dnums = nan(size(F));
for a=1:numel(F)
    [s,e] = regexp(F(a).name,'\d\d\d\d\d\d\d\d');
    dnums(a) = datenum(F(a).name(s:e),'yyyymmdd');
end
end

function F = cut_down_by_date(F, start_date, end_date)
sdnum = datenum(start_date);
ednum = datenum(end_date);
fdates = get_file_datenums(F);
xx = fdates >= sdnum & fdates <= ednum;
F = F(xx);
end

function [Data, fill_vals] = load_data(filename, fields, varname)
[~,~,fileext] = fileparts(filename);
if strcmp(fileext,'.mat')
    D = load(filename, varname);
    Data = D.(varname);
    fill_vals = nan(numel(fields, 1));
elseif strcmp(fileext,'.hdf')
    [Data, fill_vals] = prod_test_load_hdf(filename, fields);
elseif strcmp(fileext,'.txt')
    [Data, fill_vals] = prod_test_load_txt(filename, fields);
else
    E.notimplemented('The ability to check %s files has not been implemented',fileext)
end
end
