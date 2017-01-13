function [ varargout ] = cat_sat_data( filepath, datafields, varargin )
%CAT_SAT_DATA( FILEPATH, DATAFIELDS ) Concatenates data from OMI .mat files
%   In some cases, one might wish to use satellite data from multiple days,
%   but we import OMI data and process BEHR data into daily files. This
%   function will load all the .mat files in the directory given by
%   FILEPATH and output a concatenated version of the data in the field or
%   fields given by DATAFIELDS, which should be a string or cell array of
%   strings.
%
%   CAT_SAT_DATA( DATA, DATAFIELDS ) will concatenate all swaths in the
%   structure DATA for the fields specified in DATAFIELDS.
%
%   Parameter arguments are:
%
%       'prefix' - the prefix of the file names (the part before the date).
%       This function will use all files in the directory FILEPATH that
%       match the pattern <PREFIX>*.mat, replacing <PREFIX> with the given
%       prefix. This defaults to an empty string, i.e. by default all .mat
%       files will be matched.
%
%       'startdate' and 'enddate' - any valid representation of dates
%       (datenum or datestr). If only one is given, it will be treated as a
%       start or end point with the other unspecified (e.g. if you specify
%       startdate as 1-Jan-2015, this will operate on all files after
%       1-Jan-2015)
%
%       'newdim' - boolean, defaults to false. When true, each variable will
%       be concatenated along a new dimension (so a 2D variable will be
%       concatenated along the third dimension, a 3D one along the fourth).
%       When false, they will be concatenated in the along track dimension.
%
%       'varname' - must be the string 'Data' or 'OMI'. 'Data' is default.
%       Indicates whether the structure concatenated should be the native
%       pixels ('Data') or the gridded pixels ('OMI').
%
%       'DEBUG_LEVEL' - set to 0 to suppress debugging messages, defaults
%       to 1. Set to 'visual' to use the waitbar dialogue.
%
%   Josh Laughner <joshlaugh5@gmail.com> 10 Sept 2015

E=JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

if isstruct(filepath)
    Data = filepath;
    load_data = false;
else
    load_data = true;
    if ~ischar(filepath) || ~exist(filepath,'dir')
        E.badinput('filepath must be a string specifying a valid directory')
    end
end

if ischar(datafields)
    datafields = {datafields};
elseif ~iscell(datafields) || any(~iscellcontents(datafields,'ischar'))
    E.badinput('datafields must be a string or cell array of strings')
end

p=inputParser;
p.addParameter('prefix','',@ischar);
p.addParameter('startdate','');
p.addParameter('enddate','');
p.addParameter('newdim',false);
p.addParameter('varname','Data');
p.addParameter('DEBUG_LEVEL',1,@(x) (ischar(x) || isnumeric(x) && isscalar(x)));

p.parse(varargin{:});
pout = p.Results;

prefix = pout.prefix;
startdate = pout.startdate;
enddate = pout.enddate;
newdim = pout.newdim;
varname = pout.varname;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

if ~ismember(varname,{'Data','OMI'})
    E.badinput('VARNAME must be either ''Data'' or ''OMI''');
end

wbbool = false;
if ischar(DEBUG_LEVEL)
    if strcmpi(DEBUG_LEVEL, 'visual')
        if isDisplay
            wbbool = true;
            DEBUG_LEVEL = 0;
        end
    else
        warning('Only the string ''visual'' for DEBUG_LEVEL will trigger the use of the waitbar.')
        DEBUG_LEVEL = 1;
    end
end

if isempty(startdate)
    startdate = 0;
else
    try
        startdate = datenum(startdate);
    catch err
        if strcmp(err.identifier,'MATLAB:datenum:ConvertDateString')
            E.badinput('startdate format was not recognized. See datestr and datenum documentation for proper formats')
        else
            rethrow(err)
        end
    end
end

if isempty(enddate)
    enddate = datenum('3000-12-31'); % sufficiently far into the future
else
    try 
        enddate = datenum(enddate);
    catch err
        if strcmp(err.identifier, 'MATLAB:datenum:ConvertDateString')
            E.badinput('enddate format was not recognized. See datestr and datenum documentation for proper formats');
        else
            rethrow(err)
        end
    end
end

if startdate > enddate
    E.badinput('startdate is later than enddate.')
end

if ~isscalar(newdim) || (~islogical(newdim) && ~isnumeric(newdim))
    E.badinput('The parameter newdim must be understood as a scalar logical.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if load_data
    % Get all .mat files in the specified directory
    F = dir(fullfile(filepath, sprintf('%s*.mat',prefix)));
    
    if isempty(F)
        E.filenotfound('satellite .mat file');
    end
else
    F = 0;
end

% Prep output
varargout = cell(1,numel(datafields));

% Loop over all files (within the date limits given). Load the data
% variable, look for the datafields given, and add their data to the output
% which will be one long column vector.

if wbbool && load_data
    wb = waitbar(0,sprintf('Concatenating %s*.mat',strrep(prefix,'_','\_')));
elseif wbbool
    wb = waitbar(0,'Concatenating input structure');
end

for a=1:numel(F)
    if load_data
        [s,e] = regexp(F(a).name, '\d\d\d\d\d\d\d\d');
        filedate = datenum(F(a).name(s:e), 'yyyymmdd');
        if filedate < startdate || filedate > enddate
            continue
        end
        
        D = load(fullfile(filepath, F(a).name),varname);
        if ~isfield(D, varname)
            fprintf('%s does not contain the variable "%s", skipping\n', F(a).name, varname);
            continue
        else
            Data = D.(varname);
        end
    
        if DEBUG_LEVEL > 0
            fprintf('Loading file %s...\n',F(a).name);
        elseif wbbool
            waitbar(a/numel(F));
        end
    end
    
    for b=1:numel(datafields)
        if ~isfield(Data,datafields{b})
            E.callError('fieldnotfound','The field %s is not present in file %s',datafields{b},F(a).name);
        end
        
        for c=1:numel(Data)
            if newdim
                n = ndims(Data(c).(datafields{b}));
                varargout{b} = cat(n+1, varargout{b}, Data(c).(datafields{b}));
            elseif ~newdim && ismatrix(Data(c).(datafields{b}))
                varargout{b} = cat(1, varargout{b}, Data(c).(datafields{b}));
            elseif ~newdim && ~ismatrix(Data(c).(datafields{b}))
                varargout{b} = cat(2, varargout{b}, Data(c).(datafields{b}));
            else
                E.notimplemented(sprintf('concat case: newdim = %d and ndims = %d',newdim,ndims(Data(c).(datafields{b}))));
            end
        end
    end
end

if wbbool
    close(wb);
end

end

