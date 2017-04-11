function [ fname ] = behr_filename( date_in, ext, any_version )
%BEHR_FILENAME Create the file name for a BEHR file of a given date
%   FNAME = BEHR_FILENAME( DATE_IN ) Given the date, DATE_IN, as a date
%   number or a string which datestr can parse without a specified format,
%   will construct the proper BEHR .mat filename.
%
%   FNAME = BEHR_FILENAME( DATE_IN, EXT ) uses the extension EXT instead of
%   .mat.

if ~exist('ext','var')
    ext = 'mat';
else
    % Remove a leading "."; it is included already
    ext = regexprep(ext, '^\.', '');
end

if ~exist('any_version', 'var')
    any_version = false;
elseif (~islogical(any_version) && ~isnumeric(any_version)) || ~isscalar(any_version)
    E.badinput('ANY_VERSION must be a scalar number or boolean (if given)')
end

if any_version
    ver_str = '*';
else
    ver_str = BEHR_version();
end

fname = sprintf('OMI_BEHR_%s_%s.%s', ver_str, datestr(date_in, 'yyyymmdd'), ext);

end

