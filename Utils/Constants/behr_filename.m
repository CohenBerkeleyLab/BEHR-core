function [ fname ] = behr_filename( date_in, ext )
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
    regexprep(ext, '^\.', '');
end

fname = sprintf('OMI_BEHR_%s_%s.%s', BEHR_version, datestr(date_in, 'yyyymmdd'), ext);

end

