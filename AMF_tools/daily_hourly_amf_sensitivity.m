function [ AllDaysStruct, padded_array, start_hours ] = daily_hourly_amf_sensitivity( campaign_name, plot_bool, varargin )
%DAILY_HOURLY_AMF_SENSITIVITY Generate hour-by-hour daily percent AMF differences over a campaign
%   This is essentially a wrapper function to bin_profile_by_start_time and
%   hourly_amf_sensitivity that creates a box-and-whisker plot where the
%   day-to-day variability in the % change AMF is shown for each hour. This
%   will only calculate for a single SZA, VZA, and albedo (to make this not
%   take forever).
%
%   Inputs: 
%       campaign_name - (required) a string of a DISCOVER-AQ campaign that
%       will be recognized by merge_field_names.
%
%       plot_bool - (optional) a boolean that says whether to make the plot
%       (true) or just output the structure with everything you'd need for
%       the plot. If set to false and there are no outputs, this function
%       will error.
%
%       topextrap - (parameter) value is passed on to
%       bin_profile_by_start_time. Allowed values are: 'none', 'fit',
%       'median', 'wrf', 'wrf-scaled'. Defaults to 'wrf' meaning the bins
%       above the aircraft profile will be filled in the the nearest
%       WRF-Chem generated profile.
%
%       bottomextrap - (parameter) value is passed on to
%       bin_profile_by_start_time. Allowed values are: 'none', 'fit',
%       'ground', and 'median'. Defaults to median, since this seems to do
%       alright capturing the ground [NO2].
%
%       Josh Laughner <joshlaugh5@gmail.com> 22 May 2015

E = JLLErrors;
p = inputParser;
DEBUG_LEVEL = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(1,Inf);

if ~ischar(campaign_name) 
    E.badinput('campaign_name must be a string');
elseif isempty(regexpi(campaign_name,'discover'))
    E.badinput('campaign_name must be one of the DISCOVER campaigns - this function is not designed for any other campaigns')
end

if nargin < 2
    plot_bool = true;
elseif ~plot_bool && nargout < 1
    E.callError('no_output','Not plotting and no output variable - need at least one output variable to use plot_bool == false');
end

p.addParameter('topextrap', 'wrf', @(x) ismember(x,{'median','fit','wrf','wrf-scaled','none'}));
p.addParameter('bottomextrap', 'median', @(x) ismember(x, {'median','fit','ground','none'}));
p.parse(varargin{:});
pout = p.Results;
topextrap = pout.topextrap;
bottomextrap = pout.bottomextrap;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

[~,dates] = merge_field_names(campaign_name);

dates_vec = datenum(dates{1}):datenum(dates{2});

D=0;
start_hours = [];
for d=1:numel(dates_vec)
    curr_date = datestr(dates_vec(d),29);
    curr_month = str2double(curr_date(6:7));
    if DEBUG_LEVEL > 0; fprintf('\t Now loading %s\n', curr_date); end
    % We'll only look at the bins 8-9a through 3-4p since this cuts out
    % some profiles at odd times that I'm not sure are being assigned to
    % the right time.
    TodaysStruct = bin_profile_by_start_time(campaign_name, curr_date, curr_date, 8, 15, 'topextrap', topextrap, 'bottomextrap', bottomextrap);
    if isempty(TodaysStruct)
        % If there's no file for today, bin_profile_by_start_time will
        % return an empty matrix, so we need to skip the rest of the loop
        continue
    end
    TodaysStruct.date = curr_date;
    
    % Compare the data from each site against the WRF profile interpolated
    % to its coordinates.
    fns = fieldnames(TodaysStruct);
    for f=1:numel(fns)
        if ~isempty(regexp(fns{f},'Site','ONCE'))
            [perdiff, amf, amf0, vectors] = hourly_amf_sensitivity(TodaysStruct.(fns{f}), 'wrf', curr_month, 'n_sza', 1, 'n_vza', 1, 'n_alb', 1);
            new_fn = strcat(fns{f},'_AMF');
            TodaysStruct.(new_fn).perdiff = perdiff;
            TodaysStruct.(new_fn).amf = amf;
            TodaysStruct.(new_fn).amf0 = amf0;
            TodaysStruct.(new_fn).vectors = vectors;
            start_hours = cat(2,start_hours,TodaysStruct.(fns{f}).bin_start_hours);
        end
    end
    
    D=D+1;
    
    try
        AllDaysStruct(D) = TodaysStruct;
    catch err
        if ~strcmp(err.identifier, 'MATLAB:heterogeneousStrucAssignment')
            rethrow(err)
        else
            % Some campaigns will have sites present only on certain days.
            % In that case, there will be different fields in the
            % AllDaysStruct and TodaysStruct structures. This will identify
            % which structure is missing fields, what fields are missing,
            % and fill them in.
            fns_alldays = fieldnames(AllDaysStruct);
            fns_today = fieldnames(TodaysStruct);
            n_alldays = numel(fns_alldays);
            n_today = numel(fns_today);
            if n_alldays < n_today
                xx = ~ismember(fns_today, fns_alldays);
                AllDaysStruct = add_fields(AllDaysStruct, fns_today(xx));
                AllDaysStruct = orderfields(AllDaysStruct, TodaysStruct);
            elseif n_alldays > n_today
                xx = ~ismember(fns_alldays, fns_today);
                TodaysStruct = add_fields(TodaysStruct, fns_alldays(xx));
                TodaysStruct = orderfields(TodaysStruct, AllDaysStruct);
            else % if there are equal numbers of fields, but the mismatch error occured, that means that each has a field not in the other, which needs special handling
                xx = ~ismember(fns_today, fns_alldays);
                AllDaysStruct = add_fields(AllDaysStruct, fns_today(xx));
                xx = ~ismember(fns_alldays, fns_today);
                TodaysStruct = add_fields(TodaysStruct, fns_alldays(xx));
                fns_all = unique(cat(1,fns_today, fns_alldays));
                AllDaysStruct = orderfields(AllDaysStruct);
                TodaysStruct = orderfields(TodaysStruct);
            end
            AllDaysStruct(D) = TodaysStruct;
        end
    end
end

% Order the fields alphabetically, if they haven't been already.
AllDaysStruct = orderfields(AllDaysStruct);

start_hours = sort(unique(start_hours));
all_perdiffs = cell(size(start_hours));
fns = fieldnames(AllDaysStruct);
xx = ~iscellcontents(regexp(fns,'Site\d\d_AMF'),'isempty');
fns = fns(xx);

% Concatenate all the percent differences from different days for the same
% hour bin
for a=1:numel(AllDaysStruct)
    for f=1:numel(fns)
        if ~isempty(AllDaysStruct(a).(fns{f})) % this will avoid trying to do anything with fields added in the try-catch clause above to make structures match.
            for h=1:numel(start_hours)
                xx_h = AllDaysStruct(a).(fns{f}).vectors.Hours == start_hours(h);
                if sum(xx_h) > 0
                    all_perdiffs{h} = cat(1, all_perdiffs{h}, AllDaysStruct(a).(fns{f}).perdiff(xx_h));
                end
            end
        end
    end
end

% Convert the all_perdiffs cell array to a matrix where each cell becomes a
% column, padding with NaNs as needed.
len_cells = zeros(size(all_perdiffs));
for a=1:numel(len_cells)
    len_cells(a) = length(all_perdiffs{a});
end
max_len = max(len_cells);
padded_array = nan(max_len, numel(all_perdiffs));
for a=1:numel(all_perdiffs)
    padded_array(1:len_cells(a),a) = all_perdiffs{a}(:);
end

if plot_bool
    figure;
    boxplot(padded_array);
    set(gca,'fontsize',16);
    set(gca,'xticklabel',start_hours);
    xlabel('Hour of day (LST)');
    ylabel('% Difference AMF');
    title(upper(campaign_name));
end

end

function S = add_fields(S, fields2add)
E = JLLErrors;
if ~isstruct(S)
    E.badinput('S must be a structure');
elseif ~iscell(fields2add)
    E.badinput('field2add must be a cell array')
end

for f=1:numel(fields2add)
    for a=1:numel(S)
        S(a).(fields2add{f}) = [];
    end
end
end

