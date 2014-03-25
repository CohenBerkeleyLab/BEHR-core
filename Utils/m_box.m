function m_box( Longitude, Latitude, varargin)
%m_box      Josh Laughner <joshlaugh5@gmail.com> 21 Feb 2014
%   A function that will draw a box around a region on a map produced by
%   m_map.  Takes longitude and latitude matrices in one of two formats.
%   By default, longitude and latitude are entered as matrices of corner
%   points. By setting the parameter 'mode' = 'range', the function will
%   instead accept a min/max longitude and latitude.
p = inputParser;
p.addRequired('lon', @(x) length(x)==2 || length(x)==4);
p.addRequired('lat', @(x) length(x)==2 || length(x)==4);
p.addParamValue('mode', 'corners', @isstr);
p.addParamValue('linecolor','r', @isstr);
p.addParamValue('linestyle', '-', @isstr);
p.addParamValue('marker','none', @isstr);
p.addParamValue('linewidth', 1, @isscalar);

p.parse(Longitude, Latitude, varargin{:});
input = p.Results;
lon = input.lon; lat = input.lat;

if strcmpi(input.mode,'corners')
    lon(5) = lon(1); lat(5) = lat(1);
    m_line(lon, lat, 'color', input.linecolor, 'linestyle', input.linestyle, 'linewidth', input.linewidth,'marker',input.marker);
elseif strcmpi(input.mode, 'range');
    lon_b = [lon(1) lon(1) lon(2) lon(2) lon(1)]; lat_b = [lat(1) lat(2) lat(2) lat(1) lat(1)];
        m_line(lon_b, lat_b, 'color', input.linecolor, 'linestyle', input.linestyle, 'linewidth', input.linewidth, 'marker',input.marker);
else
    disp('Input mode not recognized');
end
end

