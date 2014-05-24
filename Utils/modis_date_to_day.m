function [ day ] = modis_date_to_day( date_in )
%modis_day_to_date: Convert a regular date string to the "Julian" day used
% in MODIS filenames
%   Josh Laughner 30 Apr 2014 <joshlaugh5@gmail.com>

%date = datestr(day + datenum(year-1,12,31));
d = datestr(date_in,29); %Convert any input datestr to the same format
year = str2double(d(1:4));
day = datenum(d) - datenum(year-1,12,31);

end

