function [ xx ] = omi_rowanomaly( data, parse_mode )
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

p = inputParser;
p.addRequired('data', @isstruct);
p.addOptional('mode','AlwaysByRow',@isstr);

p.parse(data,parse_mode);
input = p.Results; omi = input.data; mode = input.mode;

switch mode
    case 'AlwaysByRow'
        xx = find((omi.Row >=27 & omi.Row<=45) | (omi.Row >= 53 & omi.Row <= 54)); %JLL 20 Mar 2014: Find the matrix indices of affected rows, return
        
    case 'RowsByTime'
        log = false(size(omi.Row)); %JLL 20 Mar 2014: Based on the date in the file, determine which rows are affected. The logical OR operator ("|") "combines" the logical matrices
        if datenum(omi.Date, 'yyyy/mm/dd') >= datenum('2007/06/25','yyyy/mm/dd')
            log = log | (omi.Row >= 53 & omi.Row <= 54);
        end
        if datenum(omi.Date, 'yyyy/mm/dd') >= datenum('2008/12/03', 'yyyy/mm/dd')
            log = log | (omi.Row >= 27 & omi.Row <= 45);
        elseif datenum(omi.Date, 'yyyy/mm/dd') >= datenum('2008/05/11', 'yyyy/mm/dd')
            log = log | (omi.Row >= 27 & omi.Row <= 42);
        end
        xx = find(log);
        
    case 'XTrackFlags'
        xx = true(size(omi.XTrackQualityFlags));
        for a=1:numel(omi.XTrackQualityFlags)
            if any(omi.XTrackQualityFlags{a}) ~= 0; %JLL 20 Mar 2014: If the XTrackQualityFlag value is not 0, then the pixel has been affected by the row anomaly.
                xx(a) = 0;
            end
        end
    case 'XTrackFlagsLight'
        xx = true(size(omi.XTrackQualityFlags));
        for a=1:numel(omi.XTrackQualityFlags)
            binary_flags = de2bi(omi.XTrackQualityFlags{a}); %JLL 20 Mar 2014: Convert all values to binary arrays.  Each value will be a row in the resulting matrix.
            key_flags = bi2de(binary_flags(:,1:3)); %JLL 20 Mar 2014: Convert the first three (least significant) bits back to decimal
            xx(a) = key_flags == 1 | key_flags == 7;
        end
end

end

