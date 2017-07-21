function [ fname ] = sp_savename( date_in, ext, any_version )
%SP_SAVENAME Create the file name for an SP file of a given date
%   FNAME = SP_SAVENAME( DATE_IN ) Given the date, DATE_IN, as a date
%   number or a string which datestr can parse without a specified format,
%   will construct the proper SP .mat filename for the file which contains
%   the SP, pixel corner, modis, and globe data.
%
%   FNAME = SP_SAVENAME( DATE_IN, EXT ) uses the extension EXT instead of
%   .mat.

E = JLLErrors;

validate_date(date_in);

if ~exist('ext','var')
    ext = 'mat';
else
    if ~ischar(ext)
        E.badinput('EXT (if given) must be a string');
    end
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

fname = sprintf('OMI_SP_%s_%s.%s', ver_str, datestr(date_in, 'yyyymmdd'), ext);

end

