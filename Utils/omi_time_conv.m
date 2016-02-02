function [ dnum ] = omi_time_conv( sat_time )
%OMI_TIME_CONV Calculates the UTC time of OMI observations
%   OMI pixels give the time as time at the start of the scan (each row
%   that is) in TAI93, or the number of seconds since midnight, Jan 1st
%   1993 with some leap seconds thrown in. Leap seconds were added on:
%
%       1 Jul 1993
%       1 Jul 1994
%       1 Jan 1996
%       1 Jul 1997
%       1 Jan 1999
%       1 Jan 2006
%       1 Jan 2009
%       1 Jul 2012
%       1 Jul 2015
%
%   This will convert a vector of times from TAI93 to a datenum with time
%   in UTC.

E = JLLErrors;

if ~isnumeric(sat_time)
    E.badinput('sat_time must be a numeric input');
end

% Figure out how many leap seconds there have been
leap_sec = 0;
leap_dates = {  '1993-07-01 00:00:00';...
                '1994-07-01 00:00:00';...
                '1996-01-01 00:00:00';...
                '1997-07-01 00:00:00';...
                '1999-01-01 00:00:00';...
                '2006-01-01 00:00:00';...
                '2009-01-01 00:00:00';...
                '2012-07-01 00:00:00';...
                '2015-07-01 00:00:00'};
leap_tai93 = (datenum(leap_dates) - datenum('1993-01-01 00:00:00'))*3600*24;

% Matlab's datenum function assumes each day is 3600*24 seconds, it does
% not account for leap seconds. So now we will remove the leap seconds from
% the satellite time, working backward.

for a=0:(length(leap_tai93)-1)
    ind = length(leap_tai93) - a;
    sat_time(sat_time > leap_tai93(ind)) = sat_time(sat_time > leap_tai93(ind)) - 1;
end

% Now we can simply take this as seconds since midnight, Jan 1, 1993
% assuming 3600*24 seconds per day.

dnum = sat_time/(3600*24) + datenum('1993-01-01 00:00:00');

end

