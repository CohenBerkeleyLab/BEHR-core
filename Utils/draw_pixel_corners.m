function [  ] = draw_pixel_corners( loncorn, latcorn, varargin )
%draw_pixel_corners Draws the corners of satellite pixels
%
%   Takes a 4 x n matrix for longitude and latitude of corners, plus
%   parameter that can be passed to the line command.

if size(loncorn,1)~=4 || size(latcorn,1)~=4;
    error('draw_corners:input','First dimension of input matrices must have length 4');
end

% Create a group to hold the corners
cgroup = hggroup('Tag','Corners');

% Reorder the corners to draw a polygon without lines crossing inside
% (credit to http://stackoverflow.com/questions/13935324/sorting-clockwise-polygon-points-in-matlab)

x = loncorn(:,1); y = latcorn(:,1);
a = atan2(y - mean(y), x - mean(x));
[~,order] = sort(a);

loncorn = loncorn(order,:);
latcorn = latcorn(order,:);

s = size(loncorn);
loncorn2 = [loncorn; loncorn(1,:); nan(1,s(2:end))];
latcorn2 = [latcorn; latcorn(1,:); nan(1,s(2:end))];

line(loncorn2, latcorn2, 'parent',cgroup, varargin{:});

end

