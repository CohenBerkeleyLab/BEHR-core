function [windspd, winddir] = calc_avg_wind(xlon, xlat, U, V, clon, clat, radius)
% CALC_AVG_WIND Computes avg wind about a point
%   [WINDSPD, WINDDIR] = CALC_AVG_WIND( XLON, XLAT, U, V, CLON, CLAT )
%   returns vectors WINDSPD and WINDDIR that contain the average wind speed
%   and direction around CLON and CLAT for each day in U and V.
%
%       XLON and XLAT must be 2D matrices that specify the longitude and
%       latitude coordinates for U and V. 
%
%       U and V should be arrays of wind vectors (x and y respectively).
%       The first two dimensions should correspond to longitude and
%       latitude, and have the same lengths as the two dimensions of XLON
%       and XLAT. U and V may have an arbitrary number of dimensions, so
%       long as the last one is the time dimension, that is, if U and V are
%       4D, then U(:,:,:,1) and V(:,:,:,1) are the slices for the first
%       day, U(:,:,:,2) and V(:,:,:,2) the second day and so on. U and V
%       may be input staggered in the WRF coordinate sense, they will be
%       unstaggered internally.
%
%       CLON and CLAT should be scalar numbers giving the center point
%       around which to average wind speed and direction. By default, a 3x3
%       square of grid points around this point is averaged.
%
%   [WINDSPD, WINDDIR] = CALC_AVG_WIND( ____, RADIUS ) allows you to modify
%   how far from CLON and CLAT the wind field is averaged. If not given, it
%   defaults to 1, giving a 3x3 box of grid points.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

if ~exist('radius', 'var')
    radius = 1;
elseif ~isnumeric(radius) || ~isscalar(radius) || mod(radius,1) ~= 0
    E.badinput('RADIUS must be a scalar integer number')
end

if ~isequal(size(xlon), size(xlat))
    E.badinput('XLON and XLAT must be the same size');
elseif ~ismatrix(xlon) || ~ismatrix(xlat) || any(size(xlon) < 2*radius+1)
    % find_square_around relies on these having the points arranged in 2D
    % as they are physically, so there needs
    E.badinput('XLON and XLAT must 2D matrices with at least %d elements in both dimensions for a radius of %d\n(The matrices must have 2*radius+1 elements in each dimension)')
end



% Unstagger U and V if necessary
if size(U,1) == size(xlon,1) + 1
    U = unstagger(U,1);
end
if size(V,2) == size(xlon,2) + 1
    V = unstagger(V,2);
end

sz_U = size(U);
sz_V = size(V);
sz = size(xlon);
if any(sz_U(1:2) ~= sz(1:2)) || any(sz_V(1:2) ~= sz(1:2))
    E.badinput('U and V must be the same size as XLON after being unstaggered')
end
   
if ~isequal(sz_U, sz_V)
    E.badinput('U and V must be the same size as each other after being unstaggered')
end

if ~isnumeric(clon) || ~isscalar(clon) || ~isnumeric(clat) || ~isscalar(clat)
    E.badinput('CLON and CLAT must both be numeric scalar values')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

timedim = ndims(U);
n_times = size(U,timedim);
windspd = nan(1, size(U,timedim));
winddir = nan(1, size(U,timedim));

% Put the time dimension first so that U and V can have an arbitrary number
% of dimensions.
permvec = 1:ndims(U);
permvec(timedim) = [];
permvec = [timedim, permvec];

U = permute(U, permvec);
V = permute(V, permvec);

[xx1,yy1] = find_square_around(xlon, xlat, clon, clat, 1);

for a=1:n_times
    Uslice = U(a, xx1, yy1, :);
    Vslice = V(a, xx1, yy1, :);
    Ubar = nanmean(reshape(Uslice,[],1));
    Vbar = nanmean(reshape(Vslice,[],1));
    windspd(a) = calc_wind_mag(Ubar, Vbar);
    winddir(a) = calc_wind_dir(Ubar, Vbar);
end
end

function mag = calc_wind_mag(U,V)
mag = (U.^2 + V.^2).^0.5;
end

function theta = calc_wind_dir(U,V)
theta = nan(size(U));
for b=1:numel(U)
    theta(b) = atan2d(V(b),U(b));
end
end