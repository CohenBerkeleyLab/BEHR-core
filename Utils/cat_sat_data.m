function [ varargout ] = cat_sat_data( filepath, datafields, varargin )
%CAT_SAT_DATA(FILEPATH, DATAFIELDS) Concatenates data from OMI .mat files
%   In some cases, one might wish to use satellite data from multiple days,
%   but we import OMI data and process BEHR data into daily files. This
%   function will load all the .mat files in the directory given by
%   FILEPATH and output a concatenated version of the data in the field or
%   fields given by DATAFIELDS, which should be a string or cell array of
%   strings.
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
%       'DEBUG_LEVEL' - set to 0 to suppress debugging messages, defaults
%       to 1.
%
%   Josh Laughner <joshlaugh5@gmail.com> 10 Sept 2015

E=JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ischar(filepath) || ~exist(filepath,'dir')
    E.badinput('filepath must be a string specifying a valid directory')
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
p.addParameter('DEBUG_LEVEL',1,@(x) (isnumeric(x) && isscalar(x)));

p.parse(varargin{:});
pout = p.Results;

prefix = pout.prefix;
startdate = pout.startdate;
enddate = pout.enddate;
newdim = pout.newdim;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

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

% Get all .mat files in the specified directory 
F = dir(fullfile(filepath, sprintf('%s*.mat',prefix)));

if isempty(F)
    E.filenotfound('satellite .mat file');
end

% Prep output
varargout = cell(1,numel(datafields));

% Loop over all files (within the date limits given). Load the data
% variable, look for the datafields given, and add their data to the output
% which will be one long column vector.

for a=1:numel(F)
    [s,e] = regexp(F(a).name, '\d\d\d\d\d\d\d\d');
    filedate = datenum(F(a).name(s:e), 'yyyymmdd');
    if filedate < startdate || filedate > enddate
        continue
    end
    
    load(fullfile(filepath, F(a).name),'Data'); % brings the variable Data into the workspace
    
    if DEBUG_LEVEL > 0
        fprintf('Loading file %s...\n',F(a).name);
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

end

