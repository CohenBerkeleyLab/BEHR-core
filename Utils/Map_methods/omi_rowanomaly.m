function [ xx ] = omi_rowanomaly( xtrack_flags, rows, date_in, parse_mode )
%Returns matrix indicies affected by the OMI row anomaly. data is the OMI(d) data structure. Modes = 'AlwaysByRow'; 'RowsByTime'; 'XTrackFlags'; 'XtrackFlagsLight'
%   The OMI row anomaly is a progressive impairment to OMI measurments
%   starting in June 2007.  This function will return matrix indicies of
%   pixels deemed to have been affected by the row anomaly.  The modes
%   available are:
%       'AlwaysByRow' = Will return the indices of any row affected by the
%       anomaly, even before the anomaly began to affect those rows.
%       'RowsByTime' = Returns the indices of rows if the date is later
%       than those rows began to be affected.
%       'XTrackFlags' & 'XTrackFlagsLight' = Uses the standard product value
%       XTrackQualityFlags to determine which pixels are unsafe to use.
%       The first will return any pixel with the flag not equal to 0; the
%       second will only return pixels with severe errors.
%   By default, this function will use 'AlwaysByRows,' which is the most
%   restrictive method.
%
%   See http://www.knmi.nl/omi/research/product/rowanomaly-background.php
%   for information about the OMI row anomaly.

E = JLLErrors;

if ~exist('parse_mode', 'var')
    parse_mode = 'XTrackFlags';
end

switch parse_mode
    case 'AlwaysByRow'
        xx = find((rows >=27 & rows<=45) | (rows >= 53 & rows <= 54)); %JLL 20 Mar 2014: Find the matrix indices of affected rows, return
        
    case 'RowsByTime'
        log = false(size(rows)); %JLL 20 Mar 2014: Based on the date in the file, determine which rows are affected. The logical OR operator ("|") "combines" the logical matrices
        if datenum(date_in, 'yyyy/mm/dd') >= datenum('2007/06/25','yyyy/mm/dd')
            log = log | (rows >= 53 & rows <= 54);
        end
        if datenum(date_in, 'yyyy/mm/dd') >= datenum('2008/12/03', 'yyyy/mm/dd')
            log = log | (rows >= 27 & rows <= 45);
        elseif datenum(date_in, 'yyyy/mm/dd') >= datenum('2008/05/11', 'yyyy/mm/dd')
            log = log | (rows >= 27 & rows <= 42);
        end
        xx = find(log);
        
    case 'XTrackFlags'
        xx = xtrack_flags ~= 0;
    case 'XTrackFlagsLight'
        xx = false(size(xtrack_flags));
        for a=1:numel(xtrack_flags)
            binary_flags = de2bi(xtrack_flags(a)); %JLL 20 Mar 2014: Convert all values to binary arrays.  Each value will be a row in the resulting matrix.
            key_flags = bi2de(binary_flags(:,1:3)); %JLL 20 Mar 2014: Convert the first three (least significant) bits back to decimal
            xx(a) = key_flags == 1 | key_flags == 7;
        end
    otherwise
        allowed_parse_modes = {'AlwaysByRow', 'RowsByTime', 'XTrackFlags', 'XTrackFlagsLight'};
        E.badinput('PARSE_MODE == "%s" is not valid, allowed values are: %s', parse_mode, strjoin(allowed_parse_modes, ', '));
end

end

