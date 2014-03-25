function [ date ] = modis_day_to_date( day, year )
%modis_day_to_date: Convert the "Julian" day given in MODIS (MCD43C3) file
%names to a regular date string.
%   Josh Laughner 20 Mar 2014 <joshlaugh5@gmail.com>

date = datestr(day + datenum(year-1,12,31));

end

