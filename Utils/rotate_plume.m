function [ OMI ] = rotate_plume( Data, center_lon, center_lat, theta, varargin )
%ROTATE_PLUME Code to rotate NO2 plumes to align them by wind direction
%   Required inputs:
%       Data - should be the Data structure created by read_omno2 or
%       BEHR_main.m Only pass one top-level element at a time, i.e. if Data
%       is a 1x4 structure, pass in Data(i) where i is 1, 2, 3, or 4.
%
%       center_lon/center_lat - the center coordinates that the rotation
%       should be carried out around. Must be scalar numbers between -180
%       and 180 and -90 and 90 respectively.
%
%       theta - the angle the plume makes relative to a vector due east
%       from the source. That is, a plume being advected to the northwest
%       will have a theta of ~135 degrees.  This should be degrees (not
%       radians). This might be derived based on the wind direction over
%       the source, the major axis of an ellipse fit to the extent of the
%       plume, or another method.
%
%   Optional inputs:
%       rel_box_corners - a 1x4 vector setting the size of the box to
%       contain the plume. The values are: degrees west of center, degrees
%       east of center, degrees south of center, degrees north of center.
%       Defaults to [2 4 2 2], i.e. the box will be 6 deg E/W by 4 deg N/S
%       offset so that it extends 2 deg further E than west.
%
%   Original code from Luke Valin:
%
%         xBOX = [-2 4 4 -2];
%         yBOX = [-2 -2 2  2];
%         V = hypot(windX(date_i),windY(date_i));
%         T = theta(date_i);
%         for corner = 1 : 4
%          out =   [cos(T) -sin(T); sin(T) cos(T) ]*[xBOX(corner); yBOX(corner)];
%          xBp(corner)= out(1)+pp_lon;
%          yBp(corner)= out(2)+pp_lat;
%         end
%
%         v_ind = find(inpolygon(swath_lon,swath_lat,xBp,yBp))
%
%                          if length(v_ind)<=1
%                             vcd_tot(date_i,:,:)=nan(161,241);
%                             continue
%                          end
%
%
%         swath_X = nan(1,length(v_ind));
%         swath_Y = nan(1,length(v_ind));
%         for t_i =1:length(v_ind)
%
%             out = [cos(-T) -sin(-T); sin(-T) cos(-T) ]*[swath_lon( v_ind(t_i))-pp_lon; swath_lat(v_ind(t_i))-pp_lat];
%             swath_X(t_i)=out(1);
%             swath_Y(t_i)=out(2);
%
%         end
%
%   The idea is to define a 4 deg lat x 6 deg lon box that is assumed to
%   capture a plume.  It is defined along the x-axis, with 2 deg to the
%   west of the center point and 4 deg east. This is then rotated to align
%   with the wind direction defined by theta (which should be degrees CCW
%   from east, i.e. a normal definition of theta in polar coordinates). It
%   then identifies which pixels fall in that box and rotates them back to
%   the x-axis.
%
%   There are two potential ways to implement this with BEHR data. The
%   simplest approach would be to rotate the gridded version in the OMI
%   structure and simply assume that rotating the borders of the grid cells
%   is unecessary.  The potential problem with that is that it will not
%   ensure that the resulting grid lies well on the x-y plane, which would
%   make integration across the plume to get a line density more difficult.
%
%   The better solution I think is to rotate the native pixels and then
%   grid them, thus ensuring that the grid is defined along the x-y axes.
%   This will make the following steps easier, but will require some more
%   careful identification of which pixels fall into the box and rotation
%   of the pixel corners as well.
%
%   To identify which pixels fall into the box, we use another utility
%   function I've written which checks if box A has any overlap with box B.
%   This will capture any case of overlap.  Once the pixels contained in
%   the rotated box are identified, their longitude and latitude centers
%   and corners will be rotated and the gridded within the unrotated box.
%
%   Josh Laughner <joshlaugh5@gmail.com> 4 Feb 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addOptional('rel_box_corners',[2 4 2 2]);
p.addParameter('vza_crit',60);
p.parse(varargin{:});
pout = p.Results;
rel_box_corners = pout.rel_box_corners;
vza_crit = pout.vza_crit;

E = JLLErrors;

if ~isstruct(Data) || ~isscalar(Data) || any(~ismember({'Longitude','Latitude','Loncorn','Latcorn'}, fieldnames(Data)))
    E.badinput('Data must be a scalar structure with the fields Longitude, Latitude, Loncorn, and Latcorn')
end
if ~isscalar(center_lon) || ~isnumeric(center_lon) || center_lon > 180 || center_lon < -180
    E.badinput('center_lon must be a scalar numeric value between -180 and +180')
elseif ~isscalar(center_lat) || ~isnumeric(center_lat) || center_lat > 180 || center_lat < -180
    E.badinput('center_lat must be a scalar numeric value between -90 and +90')
elseif ~isscalar(theta) || ~isnumeric(theta) || theta > 180 || theta < -180
    E.badinput('theta must be a scalar numeric value between -180 and +180')
end


if numel(rel_box_corners) ~= 4 || ~isnumeric(rel_box_corners) || ~isvector(rel_box_corners)
    E.badinput('rel_box_corners (if given) must be a 4 element numeric vector')
end
[x_box, y_box, lonlim, latlim] = convert_rel_box_corners(rel_box_corners, center_lon, center_lat);


if ~isscalar(vza_crit) || ~isnumeric(vza_crit)
    E.badinput('vza_crit must be a numeric scalar')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% First rotate the box to align along theta and find any pixels that
% overlap it at all.
R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
x_box_rot = nan(size(x_box));
y_box_rot = nan(size(y_box));
for corner=1:4
    out = R * [x_box(corner); y_box(corner)];
    x_box_rot(corner) = out(1) + center_lon;
    y_box_rot(corner) = out(2) + center_lat;
end

% Remove pixels with VZA greater than the specified criteria
% which defaults to 60 degrees.
vv = Data.ViewingZenithAngle <= vza_crit;

% Remove pixels outside the box
pp = false(size(Data.Longitude));
for a=1:numel(Data.Longitude)
    pp(a) = test_box_overlap(Data.Loncorn(:,a), Data.Latcorn(:,a), x_box_rot, y_box_rot);
end

xx = pp & vv;

fns = fieldnames(Data);
fns(ismember(fns,{'Longitude','Latitude','Loncorn','Latcorn'})) = []; % we handle this fields separately, so don't include them in the fieldnames loops

Data.Longitude(~xx) = [];
Data.Latitude(~xx) = [];
Data.Loncorn(:,~xx) = [];
Data.Latcorn(:,~xx) = [];
fields_to_remove = {};
for f=1:numel(fns)
    if ~ismatrix(Data.(fns{f}))
        fields_to_remove{end+1} = fns{f}; %#ok<AGROW>
    elseif isnumeric(Data.(fns{f})) && all(size(xx) == size(Data.(fns{f})))
        Data.(fns{f})(~xx) = [];
    end
end

for f=1:numel(fields_to_remove)
    Data=rmfield(Data,fields_to_remove{f});
    fns(strcmp(fns,fields_to_remove{f})) = [];
end
% Rotate the pixels back to the x-axis (only need to change the lon/lat fields)
R = [cosd(-theta), -sind(-theta); sind(-theta), cosd(-theta)];
for a=1:numel(Data.Longitude)
    out = R * [Data.Longitude(a) - center_lon; Data.Latitude(a) - center_lat];
    Data.Longitude(a) = out(1) + center_lon;
    Data.Latitude(a) = out(2) + center_lat;
    for b=1:4
        out = R * [Data.Loncorn(b,a) - center_lon; Data.Latcorn(b,a) - center_lat];
        Data.Loncorn(b,a) = out(1) + center_lon;
        Data.Latcorn(b,a) = out(2) + center_lat;
    end
end


% Finally grid the data to a 0.05 x 0.05 degree grid.

resolution = 0.05; resolution2 = 0.05;
% This will be passed to the gridding function to keep the field order
% correct.
OMI = struct('BEHRColumnAmountNO2Trop', [], 'ViewingZenithAngle', [], 'SolarZenithAngle', [], 'AMFTrop', [], 'CloudFraction', [], 'CloudRadianceFraction', [],...
    'CloudPressure', [], 'ColumnAmountNO2Trop', [], 'RelativeAzimuthAngle', [], 'MODISAlbedo', [], 'GLOBETerpres', [], 'BEHRAMFTrop', [],...
    'Latitude', [], 'Longitude', [], 'MapData', struct, 'Count', [], 'Area', [], 'Areaweight', [], 'vcdQualityFlags', {{}}, 'XTrackQualityFlags', {{}});
OMI = repmat(OMI,1,numel(Data));
hh=0;

if all(Data.ViewingZenithAngle(:) == 0) || numel(Data.ViewingZenithAngle) == 1
    return
else
    OMI = add2grid_BEHR_winds(Data,OMI,resolution,resolution2,lonlim,latlim);
end



end

