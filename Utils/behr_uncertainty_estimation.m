function [ Delta, DeltaGrid ] = behr_uncertainty_estimation( Data, OMI, parameter, percent_change, varargin )
%BEHR_UNCERTAINTY_ESTIMATION Estimate the uncertainty in BEHR NO2
%   [ DELTA, DELTAGRID ] = BEHR_UNCERTAINTY_ESTIMATION( DATA, PARAMETER, PERCENT_CHANGE )
%   This function will run the BEHR retrieval for the structure DATA but
%   with the field PARAMETER changed by PERCENT_CHANGE percent. PARAMETER
%   must match a field in DATA. This will return structures DELTA and
%   DELTAGRID which are the equivalent of Data and OMI except they will
%   have the modified parameter values and the resultant different VCDs and
%   AMFs.
%
%       PERCENT_CHANGE may be either a number or a function handle. If a
%       number, then the field given by PARAMETER is changed to:
%
%           value * (100 + percent_change)/100
%
%       before running the BEHR algorithm. If PERCENT_CHANGE is a function
%       handle instead, then it must take a scalar structure as its sole
%       input and return the value that the field PARAMETER should take on.
%       For example, to use an absolute difference in MODISAlbedo rather
%       than a percent difference, use:
%
%           PERCENT_DIFFERENCE = @(Data) Data.MODISAlbedo + 0.05;
%
%   There are two parameters:
%
%   'remove_unchanged_fields' - is a scalar logical; if true, then all
%   numeric fields except those modified or added by this function are
%   removed from Delta and DeltaGrid before returning. 
%
%   'DEBUG_LEVEL' - a scalar number indicating verbosity. Default is 2; 0
%   means no output.

E = JLLErrors;

% This function requires that rProfile_WRF.m in the BEHR-core-utils repo
% has the commit that allows it to keep the full extrapolated NO2 and
% temperature profiles instead of clipping them just below and above the
% surface pressure and tropopause pressure, respectively.
G = GitChecker;
G.addReqCommits(behr_paths.behr_utils, 'ca2faf3');

p = inputParser;
p.KeepUnmatched = true; % this should avoid errors with extra parameters to be passed through to BEHR_main_one_day
p.addParameter('remove_unchanged_fields', false);
p.addParameter('DEBUG_LEVEL', 2);

p.parse(varargin{:});
pout = p.Results;

remove_unchanged_fields = pout.remove_unchanged_fields;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

if ~isstruct(Data)
    E.badinput('DATA must be a structure')
end
if ~ischar(parameter) || (~isfield(Data, parameter) && ~any(strcmpi(parameter, {'profileloc', 'profiletime'})))
    E.badinput('PARAMETER must be a field name in DATA or the special strings "profileloc" or "profiletime"')
end
if isnumeric(percent_change) 
    if ~isscalar(percent_change)
        E.badinput('If given as a number, PERCENT_CHANGE must be a scalar');
    else
        percent_change = @(Data) Data.(parameter) * (100 + percent_change)/100;
    end
elseif ~isa(percent_change, 'function_handle')
    E.badinput('PERCENT_CHANGE must be a scalar number or a function handle');
end

Delta = Data;

for a=1:numel(Delta)
    % Vary the specified parameter. Later we'll add an ability to vary the
    % NO2 profiles in a realistic way, but for now we'll stick to just
    % varying 2D parameters
    if ~any(strcmpi(parameter, {'profileloc','profiletime'}))
        Delta(a).(parameter) = percent_change(Delta(a));
    end
end

% Now run BEHR but for the modified parameters
if strcmpi(parameter, 'profileloc')
    [Delta, DeltaGrid] = BEHR_main_one_day(Delta, 'profile_mode', Delta(1).BEHRProfileMode, 'lookup_profile', true, 'lookup_sweights', false,...
        'randomize_profile_loc', true, varargin{:});
elseif strcmpi(parameter, 'profiletime')
    [Delta, DeltaGrid] = BEHR_main_one_day(Delta, 'profile_mode', Delta(1).BEHRProfileMode, 'lookup_profile', true, 'lookup_sweights', false,...
        'randomize_profile_time', true, varargin{:});
else
    [Delta, DeltaGrid] = BEHR_main_one_day(Delta, 'profile_mode', Delta(1).BEHRProfileMode, 'lookup_profile', false, 'lookup_sweights', true, 'extra_gridding_fields', {parameter}, varargin{:});
end


% Calculate the percent differences in the NO2 columns and AMFs
for a=1:numel(Delta)
    Delta(a).PercentChangeNO2 = reldiff(Delta(a).BEHRColumnAmountNO2Trop, Data(a).BEHRColumnAmountNO2Trop)*100;
    Delta(a).PercentChangeNO2Vis = reldiff(Delta(a).BEHRColumnAmountNO2TropVisOnly, Data(a).BEHRColumnAmountNO2TropVisOnly)*100;
    Delta(a).PercentChangeAMF = reldiff(Delta(a).BEHRAMFTrop, Data(a).BEHRAMFTrop)*100;
    Delta(a).PercentChangeAMFVis = reldiff(Delta(a).BEHRAMFTropVisOnly, Data(a).BEHRAMFTropVisOnly)*100;
    
    DeltaGrid(a).PercentChangeNO2 = reldiff(DeltaGrid(a).BEHRColumnAmountNO2Trop, OMI(a).BEHRColumnAmountNO2Trop)*100;
    DeltaGrid(a).PercentChangeNO2Vis = reldiff(DeltaGrid(a).BEHRColumnAmountNO2TropVisOnly, OMI(a).BEHRColumnAmountNO2TropVisOnly)*100;
    DeltaGrid(a).PercentChangeAMF = reldiff(DeltaGrid(a).BEHRAMFTrop, OMI(a).BEHRAMFTrop)*100;
    DeltaGrid(a).PercentChangeAMFVis = reldiff(DeltaGrid(a).BEHRAMFTropVisOnly, OMI(a).BEHRAMFTropVisOnly)*100;
end

if remove_unchanged_fields
    % Keep the percent change fields, the NO2 and AMF fields themselves,
    % the quality flags (so we can ID good pixels), and the changed
    % parameter, but remove all other data fields (attribute fields will be
    % kept, i.e. any non-numeric field)
    fields_to_keep = {parameter, 'BEHRColumnAmountNO2Trop', 'BEHRAMFTrop', 'BEHRColumnAmountNO2TropVisOnly', 'BEHRAMFTropVisOnly',...
        'PercentChangeNO2', 'PercentChangeAMF', 'PercentChangeNO2Vis', 'PercentChangeAMFVis', 'BEHRQualityFlags'};
    Delta = cut_down_fields(Delta, fields_to_keep);
    DeltaGrid = cut_down_fields(DeltaGrid, [fields_to_keep, {'Areaweight'}]);
end

end

function Data = cut_down_fields(Data, fields_to_keep)
fns = fieldnames(Data);
numeric_fields = structfun(@isnumeric, Data(1));
keep_fields = ismember(fns, fields_to_keep);
fields_to_remove = fns(~keep_fields & numeric_fields);
Data = rmfield(Data, fields_to_remove);
end
