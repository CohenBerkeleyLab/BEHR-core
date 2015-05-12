function [ Data_Struct ] = behr_data_check(  )
%behr_data_check Checks for missing satellite files needed for BEHR
%   BEHR requires 3 different satellite data sets:
%       1) The Level 2 OMNO2 files
%       2) Aqua-MODIS cloud product (MYD06_L2)
%       3) 16-day combined MODIS albedo (MCD43)
%   This function will check each of these in the directories they are
%   supposed to be and return a structure indicating any missing files

% Create the error handling class
E = JLLErrors;
DEBUG_LEVEL = 1;

% Get today's date
todays_d = day(today);
todays_m = month(today);
todays_y = year(today);

% Prepare the output structure, it will have the form satellite/year
n_years = todays_y - 2004 + 1;
struct_base = cell(1,2*n_years);
for a=1:n_years
    y = 2004 + a - 1;
    i = a*2 - 1; % edit the odd indicies of the cell
    struct_base{i} = sprintf('Y_%d',y);
end
Data_Struct = struct('OMNO2',struct(struct_base{:}),'MYD06',struct(struct_base{:}),'MCD43C3',struct(struct_base{:}));

% First, check the OMNO2 data. OMI data starts around Oct 2004.  The
% subfunction checkDaily allows the same code to be used for MODIS Cloud
% and OMI. Look for at least 13 files per day.
if DEBUG_LEVEL > 0; fprintf('\n***** Checking OMNO2 data. *****\n\n'); end
omno2dir = omno2_dir();
omno2pat_fxn = @(s) omiPat(s);
omno2_pathbuilder_fxn = @(sat,yr,mn) omiPathBuilder(sat,yr,mn);
omno2_chk_fxn = @(f) omiDataCheck(f);

Data_Struct = checkDaily(Data_Struct, todays_y, todays_m, omno2pat_fxn, omno2_pathbuilder_fxn, omno2dir, 'OMNO2', omno2_chk_fxn, DEBUG_LEVEL);

% Next, MODIS cloud data. There should be at least 18 files per day, but
% we'll relax that to 17 to be careful.
if DEBUG_LEVEL > 0; fprintf('\n***** Checking MYD06 data. *****\n\n'); end
modclddir = modis_cloud_dir;
modcldpat_fxn = @(s) modisCldPat(s);
modis_pathbuilder_fxn = @(sat,yr,mn) modisPathBuilder(sat,yr,mn);
modis_chk_fxn = @(f) modisCloudDataCheck(f);

Data_Struct = checkDaily(Data_Struct, todays_y, todays_m, modcldpat_fxn, modis_pathbuilder_fxn, modclddir, 'MYD06', modis_chk_fxn, DEBUG_LEVEL);

% Finally, MODIS albedo data. 
if DEBUG_LEVEL > 0; fprintf('\n***** Checking MCD43 data. *****\n\n'); end

modalbdir = modis_albedo_dir;
modalbpath_fxn = @(s) modisAlbPat(s);
modis_pathbuilder_fxn = @(sat,yr) modisPathBuilder(sat,yr);

Data_Struct = check8Daily(Data_Struct, todays_y, todays_m, modalbpath_fxn, modis_pathbuilder_fxn, modalbdir, 'MCD43C3', 1, DEBUG_LEVEL);
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SUB FUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

function Data_Struct = checkDaily(Data_Struct, todays_y, todays_m, sat_pat_fxn, sat_pathbuilder, sat_dir, sat_field, sat_test_fxn, DEBUG_LEVEL)
% Used in the messages
statuses = {'Missing', 'Incomplete'};

for ch_year = 2004:todays_y
    ch_year_str = num2str(ch_year);
    if DEBUG_LEVEL > 0; fprintf('\tNow checking year %s\n',ch_year_str); end
    months = setMonths(ch_year, todays_y, todays_m);
    day_chk = setDays(ch_year, todays_y, todays_m, 1);
    
    
    for ch_month = months;
        ch_month_str = sprintf('%02d',ch_month);
        if DEBUG_LEVEL > 1; fprintf('\t\tMonth: %s\n',ch_month_str); end
        ch_path = sat_pathbuilder(sat_dir,ch_year_str,ch_month_str);
        days = 1:eomday(ch_year,ch_month);
        for ch_day = days;
            ch_day_str = sprintf('%02d',ch_day);
            % Iterate through the days of the month, checking the number of
            % files for each day.
            S.ch_year_str = ch_year_str;
            S.ch_month_str = ch_month_str;
            S.ch_day_str = ch_day_str;
            day_pat = sat_pat_fxn(S);
            FILES = dir(fullfile(ch_path,day_pat));
            
            chk = sat_test_fxn(FILES);

            % Copy this to the year long day vector for today.
            curr_date = sprintf('%s-%s-%s',ch_year_str,ch_month_str,ch_day_str);
            ind = modis_date_to_day(curr_date);
            day_chk(ind) = chk;
        end
    end
    
    % Find the contiguous blocks of values, we'll then parse this into
    % messages for the year to be saved to the structure.
    blocks = findBlock(day_chk);
    messages = cell(size(blocks,1),1);
    mi = 1;
    for m=1:numel(messages)
        st_date = modis_day_to_date(blocks(m,1),ch_year);
        end_date = modis_day_to_date(blocks(m,2),ch_year);
        status = blocks(m,3);
        if status < 2;
            messages{mi} = sprintf('%s: %s to %s',statuses{status+1},st_date,end_date);
            mi=mi+1;
        end
    end
    messages(mi:end) = [];
    if isempty(messages)
        messages = 'All data found';
    end
    
    fname = sprintf('Y_%s',ch_year_str);
    Data_Struct.(sat_field).(fname) = messages;
end
end

function Data_Struct = check8Daily(Data_Struct, todays_y, todays_m, sat_pat_fxn, sat_pathbuilder_fxn, sat_dir, sat_field, req_num_files, DEBUG_LEVEL)
% Used in the messages
statuses = {'Missing', 'Incomplete'};

for ch_year = 2004:todays_y
    ch_year_str = num2str(ch_year);
    if DEBUG_LEVEL > 0; fprintf('\tNow checking year %s\n',ch_year_str); end
    [day_chk, days] = setDays(ch_year, todays_y, todays_m, 8);
    % Here we take advantage of the fact that MATLAB functions don't need
    % as many arguments as they accept as long as the unpassed arguments
    % are not required.
    ch_path = sat_pathbuilder_fxn(sat_dir, ch_year_str);
    % MODIS albedo files are not stored by month, so we just look for each
    % day that should exist in the year folder
    for d = 1:numel(days)
        ch_day_str = sprintf('%03d',days(d));
        S.ch_year_str = ch_year_str;
        S.ch_day_str = ch_day_str;
        day_pat = sat_pat_fxn(S);
        FILES = dir(fullfile(ch_path,day_pat));
        
        % If there are at least req_num_files files, consider the day
        % to be complete. If there are some, but not 13, mark it as
        % incomplete. If there are none, mark it as empty.
        if numel(FILES) >= req_num_files;
            chk = 2;
        elseif numel(FILES) > 0 && numel(FILES) < req_num_files;
            chk = 1;
        else
            chk = 0;
        end
        
        % Copy this to the year long day vector for today.
        day_chk(d) = chk;
    end
    % Find the contiguous blocks of values, we'll then parse this into
    % messages for the year to be saved to the structure.
    blocks = findBlock(day_chk);
    messages = cell(size(blocks,1),1);
    mi = 1;
    for m=1:numel(messages)
        st_date = modis_day_to_date(blocks(m,1),ch_year);
        end_date = modis_day_to_date(blocks(m,2),ch_year);
        status = blocks(m,3);
        if status < 2;
            messages{mi} = sprintf('%s: %s to %s',statuses{status+1},st_date,end_date);
            mi=mi+1;
        end
    end
    messages(mi:end) = [];
    if isempty(messages)
        messages = 'All data found';
    end
    
    fname = sprintf('Y_%s',ch_year_str);
    Data_Struct.(sat_field).(fname) = messages;
end
end

function months = setMonths(ch_year,todays_y, todays_m)
% If we're checking 2004, only do months 10-12. If we're checking this
% year, only go up to the current month. Otherwise do all the months.
if ch_year == 2004;
    months = 10:12;
elseif ch_year == todays_y;
    months = 1:todays_m;
else
    months = 1:12;
end
end

function [day_chk, days] = setDays(ch_year, todays_y, todays_m, ndays)
% Returns day_chk - a vector of zeros and days, a sequential vector. ndays
% determines the step size between days.
if false%leapyear(ch_year)
    days = 1:ndays:366;
else
    days = 1:ndays:365;
end
day_chk = zeros(size(days));

if ch_year == todays_y
    month = sprintf('%02d',todays_m+1);
    test_date = sprintf('%d-%s-01',todays_y,month);
    last_day = modis_date_to_day(test_date);
    day_chk(last_day:end) = 3; % use 3 to indicate data that is not available, period
elseif ch_year == 2004;
    last_day = modis_date_to_day('2004-09-30');
    day_chk(1:last_day) = 3;
end
end

function pat = omiPat(S)
omno2pat = 'OMI-Aura_L2-OMNO2_%sm%s%s*.he5';
pat = sprintf(omno2pat,S.ch_year_str,S.ch_month_str,S.ch_day_str);
end

function pat = modisCldPat(S)
modcldpat = 'MYD06_L2.A%s%03d*.hdf';
modday = modis_date_to_day(sprintf('%s-%s-%s',S.ch_year_str,S.ch_month_str,S.ch_day_str));
pat = sprintf(modcldpat,S.ch_year_str,modday);
end

function pat = modisAlbPat(S)
modalbpat = 'MCD43C3.A%s%s*.hdf';
pat = sprintf(modalbpat, S.ch_year_str,S.ch_day_str);
end

function ch_path = omiPathBuilder(sat_dir,ch_year_str,ch_month_str)
    ch_path = fullfile(sat_dir,ch_year_str,ch_month_str);
end

function ch_path = modisPathBuilder(sat_dir,ch_year_str,~)
    % Because the OMI data is stored by month, and since this function is
    % passed to the checkDaily function, it needs 3 arguments even though
    % we don't use the month for modis.
    ch_path = fullfile(sat_dir,ch_year_str);
end

function chk = omiDataCheck(FILES)
% Check that there are at least 13 files for the given day.
% If so, mark it as complete. If there are some, but not 13,
% mark it as incomplete. If there are none, mark it as missing.
    req_num_files = 13;
    nF = numel(FILES);
    if nF >= req_num_files
        chk = 2;
    elseif nF > 0 && nF < req_num_files
        chk = 1;
    else
        chk = 0;
    end
end

function chk = modisCloudDataCheck(FILES)
    % Check that the expected MODIS MYD06 granules are present.
    % You can look at their borders using draw_modis_cloud_swaths
    % (in BEHR/One-off scripts currently)
    
    % First check that there are any files are given, if not, clearly
    % this is a "missing" day
    if numel(FILES)==0
        chk = 0;
        return;
    end
    
    % Define the expected times manually. These will need to be expanded
    % if BEHR is extended beyond CONUS. This should be organized so that
    % each cell of expected_times contains another cell array that has
    % vectors of times for each orbit. You can figure these out by looking
    % at the MODIS granules using draw_modis_cloud_swaths in BEHR/One off
    % scripts. You're looking for which granules (by the time identifier)
    % are needed to cover the area of interest. 
    expected_times = cell(1,16);
    expected_times{1} = {[1745, 1750], [1925, 1930], [2100, 2105, 2110]};
    expected_times{2} = {[1655], [1825, 1830, 1835], [2005, 2010, 2015], [2145, 2150]};
    expected_times{3} = {[1735, 1740], [1910, 1915], [2050, 2055]};
    expected_times{4} = {[1815, 1820], [1955, 2000], [2130, 2135, 2140]};
    expected_times{5} = {[1720, 1725], [1855, 1900, 1905], [2035, 2040, 2045]};
    expected_times{6} = {[1630], [1800, 1805, 1810], [1940, 1945, 1950], [2120, 2125]};
    expected_times{7} = {[1710, 1715], [1845, 1850, 1855], [2025, 2030], [2205, 2210]};
    expected_times{8} = {[1750, 1755], [1930, 1935], [2105, 2110, 2115]};
    expected_times{9} = {[1700], [1835, 1840], [2010, 2015, 2020], [2150, 2155]};
    expected_times{10} = {[1740, 1745], [1915, 1920, 1925], [2055, 2100]};
    expected_times{11} = {[1650], [1820, 1825, 1830], [2000, 2005], [2140, 2145]};
    expected_times{12} = {[1725, 1730], [1905, 1910], [2045, 2050], [2220, 2225, 2230]};
    expected_times{13} = {[1810, 1815], [1945, 1950, 1955], [2125, 2130]};
    expected_times{14} = {[1715, 1720], [1850, 1855, 1900], [2030, 2035], [2210, 2215]};
    expected_times{15} = {[1625], [1755, 1800, 1805], [1935, 1940], [2115, 2120]};
    expected_times{16} = {[1705], [1840, 1845], [2020, 2025], [2155, 2200, 2205]};

    % Make a list of times in the file names. Look for four numbers with .'s on either side.
    times_list = zeros(1,numel(FILES));
    for a=1:numel(FILES)
        t_ind = regexp(FILES(a).name,'\.\d\d\d\d\.');
        times_list(a) = str2num(FILES(a).name(t_ind+1:t_ind+4));
    end
    
    % Check how many of the times are present in the files. Since
    % the Aqua satellite is on a 16-day repeat cycle, the mod16
    % of the datenum can be used to indicate which set of times to use.
    % Add one to that since MATLAB indexing doesn't start at 0.
    year = FILES(1).name(11:14);
    doy = FILES(1).name(15:17);
    file_date = modis_day_to_date(doy,year);
    mod_ind = mod(datenum(file_date),16);
    
    test_times = expected_times{mod_ind+1};
    cnt = 0;
    tot = 0;
    for t=1:numel(test_times)
        % Add a "buffer" granule to the beginning and end of the swath
        test_times_t = [timemath(test_times{t}(1),-5), test_times{t}, timemath(test_times{t}(end),5)];
        
        % Keep track of how many granules we expect in each orbit
        n = numel(test_times{t});
        tot = tot + n;
        
        % Keep track of how many of the expected times are present
        times_found = ismember(test_times_t, times_list);
        cnt = cnt + sum(times_found);
    end
    
    % The status will be complete if there are as many files as expected,
    % and they are about the right times. If there are too few, it is
    % incomplete, and obviously if there are none, the day is missing.
    if cnt >= tot
        chk = 2;
    elseif cnt < tot && cnt > 0
        chk = 0;
    else
        chk = 1;
    end

end

function t = timemath(t,m)
% Add minutes to time t where the hundreds and thousands place are hours
% and the tens and ones are minutes. 
minutes = mod(t,100);
hours = t - minutes;

minutes = minutes + m;
% what must be added to the hour: >60 = add, <0 = subtract
h = floor(minutes/60);

new_minutes = mod(minutes,60);
new_hour = hours + h*100;

t = new_hour + new_minutes;
end
