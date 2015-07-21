function [ varargout ] = misc_behr_wind_plots( plttype, varargin )
%MISC_BEHR_WIND_PLOTS Various plots studying the effect of wind variation on BEHR
%   Like the other misc. plot functions, this collects several related
%   plotting functions in one place.

E = JLLErrors;


switch lower(plttype)
    case 'windfield'
    case 'windmagangle'
        plot_wind_magnitude_and_angle(varargin{:});
    case 'plotnearprof'
        plot_nearest_profile(varargin{:});
    case 'getnearprof'
        varargout{1} = get_nearest_profiles(varargin{:});
    otherwise
        fprintf('plttype not recognized\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PLOTTING FUNCTIONS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function plot_wind_magnitude_and_angle(lon, lat, U, V, center_lon, center_lat, radius)
        % Plots the average wind magnitude a direction (as degrees
        % counterclockwise from east) over all days given. Requires at
        % least 6 inputs (lon and lat arrays, wind speed arrays, center lon
        % and lat, and the optional radius - defaults to 1). U and V are
        % assumed to have dimensions lon, lat, time. Lon and lat can either
        % be single 2D matrices (dimensions lon x lat) or have the same
        % time dimension. In the former case, the same lat/lon matrix will
        % be used for all times. The radius determines how large an area to
        % average. It defaults to 1 (i.e. a 3x3 square of points to
        % average). The std. dev. of the magnitude and direction will be
        % plotted as an envelope, if radius > 0.  This function will also
        % handle if U and V are still "staggered" and so have extra entries
        % in the lon and lat directions respectively.
        
        % Check input
        if ndims(lon) ~= ndims(lat)
            E.badinput('lon and lat should have the same number of dimensions')
        elseif (ndims(lon) < 2 && ndims(lon) > 3) || (ndims(lat) < 2 && ndims(lat) > 3)
            E.badinput('lon and lat must be 2- or 3-D arrays')
        elseif ndims(U) < 2 || ndims(U) > 3 || ndims(V) < 2 || ndims(V) > 3
            E.badinput('U and V must be 2- or 3-D arrays') % granted 2-D would be boring - 1 day - but may as well
        elseif size(U,3) ~= size(V,3)
            E.badinput('U and V must be the same length in the third dimension')
        elseif ~all(size(lon) == size(lat))
            E.badinput('lon and lat must be the same size')
        elseif ~isscalar(center_lat) || ~isnumeric(center_lat)
            E.badinput('center_lat must be a numeric scalar')
        elseif ~isscalar(center_lon) || ~isnumeric(center_lon)
            E.badinput('center_lon must be a numeric scalar')
        elseif nargin >= 7 && (~isscalar(radius) || ~isnumeric(radius) || radius < 0 || mod(radius,1) ~= 0)
            E.badinput('radius must be a positive scalar integer')
        elseif nargin < 7
            radius = 1;
        end
        
        % Check if U and V are still in staggered coordinates, if so
        % average them to unstagger. Simultaneously, make lon and lat 3D if
        % they are not
        
        sz_U = size(U);
        stag_U = zeros(size(sz_U));
        stag_U(1) = 1;
        sz_V = size(V);
        stag_V = zeros(size(sz_V));
        stag_V(2) = 1;
        
        if ndims(lon) < 3 || size(lon,3) == 1 
            lon = repmat(lon,1,1,sz_U(3));
            lat = repmat(lat,1,1,sz_U(3));
        end
        sz_lonlat = size(lon);
        
        if all(sz_lonlat == sz_U - stag_U)
            U = unstagger(U,1);
        elseif ~all(sz_lonlat == sz_U) 
            E.badinput('U must be the same size as lon & lat, or 1 greater in the first dimension');
        end
        
        if all(sz_lonlat == sz_V - stag_V)
            V = unstagger(V,2);
        elseif ~all(sz_lonlat == sz_V)
            E.badinput('V must be the same size as lon & lat, or 1 greater in the second dimension');
        end
        
        mag_mean = nan(size(U,3),1);
        mag_std = nan(size(U,3),1);
        theta_mean = nan(size(U,3),1);
        theta_std = nan(size(U,3),1);
        
        for a=1:size(U,3)
            
            [xx,yy] = find_square_around(lon(:,:,a), lat(:,:,a), center_lon, center_lat, radius);
            U_a = U(xx,yy,a);
            V_a = V(xx,yy,a);
            mag = (U_a.^2 + V_a.^2).^0.5;
            theta = mod(atand(V_a./U_a),360); % puts the angle in [0, 359] rather than [-180, 180]
        
            mag_mean(a) = nanmean(mag(:));
            if radius > 0
                mag_std(a) = nanstd(mag(:));
            end
            
            theta_mean(a) = nanmean(theta(:));
            if radius > 0
                theta_std(a) = nanstd(theta(:));
            end
        end
        
        figure;
        subplot(2,1,1);
        x = 1:numel(mag_mean);
        
        plot_error_envelope_y(x, mag_mean - mag_std, mag_mean + mag_std, 'r', gcf);
        line(x, mag_mean, 'linewidth', 2, 'color', 'k', 'marker', 'o');
        
        set(gca,'fontsize',16);
        ylabel('Wind magnitude (m/s)');
        set(gca,'ygrid','on');
        
        subplot(2,1,2);
        
        plot_error_envelope_y(x, theta_mean - theta_std, theta_mean + theta_std, 'r', gcf);
        line(x, theta_mean, 'linewidth', 2, 'color', 'k', 'marker', 'o');
        
        set(gca,'fontsize',16);
        xlabel('Day');
        ylabel(sprintf('Wind direction\n(deg. CCW from east)'));
        set(gca,'ygrid','on');
    end

    function plot_nearest_profile(lon, lat, profiles, alt, target_lon, target_lat, multifig, xtext)
        % Function that will plot the profile nearest to the target
        % latitude and longitude. lon and lat should be 2-D or 3-D
        % matrices, and have dimensions longitude, latitude, (time).
        % profiles should be 3-D or 4-D and have dimensions of longitude,
        % latitude, altitude, (time).  alt should be 3- or 4-D as well, and
        % should have dimensions of longitude, latitude, altitude, (time).
        % The target lon and lat should both be scalar. This will reject
        % target lon/lats more than ~12 km outside the lon and lat
        % matrices.
        %
        % Note that having 4-D profiles does not require lon and lat to be
        % 3-D, this will assume that the same longitude and latitude will
        % be used for all times. However, it will error if 3-D lon and lat
        % are passed with a 3-D profiles, since that doesn't make sense.
        % Likewise, alt can be 3-D for multiple days, and it will assume
        % that the same altitudes old true for all days. But if passing a
        % 4-D alt array, profiles must be 4-D.
        %
        % The multifig input can be true or false, if true it plots one
        % profile per figure instead of putting all on one plot with a
        % legend. Defaults to false. Alternatively, pass a cell array of
        % strings to use in the legend when plotting multiple profiles on
        % one figure.
        %
        % The last argument is optional and is the text to be put on the x
        % label. Useful especially if making multiple plots.
        
        % Check input types and dimensions
        if ndims(lon) < 2 || ndims(lon) > 3 || ndims(lat) < 2 || ndims(lon) > 3
            E.badinput('lon and lat must be 2- or 3-D arrays')
        elseif ndims(profiles) < 3 || ndims(profiles) > 4 
            E.badinput('profiles must be a 3- or 4-D array')
        elseif ndims(profiles) == 3 && (ndims(lon) > 2 || ndims(lat) > 2) %#ok<*ISMAT>
            E.badinput('A 3-D lon or lat array implies multiple days, but the dimensionality of profiles does not')
        elseif ~isscalar(target_lat) || ~isnumeric(target_lat) || ~isscalar(target_lon) || ~isnumeric(target_lon)
            E.badinput('The target_lon and target_lat inputs must be scalar numeric values')
        end
        
        sz_lon = size(lon);
        sz_lat = size(lat);
        sz_prof = size(profiles);
        sz_alt = size(alt);
        
        if ndims(lon) ~= ndims(lat)
            E.badinput('lon and lat must have the same number of dimensions')
        elseif ~all(sz_lon == sz_lat)
            E.badinput('lon and lat must be the same size')
        elseif ~all(sz_lon(1:2) == sz_prof(1:2))
            E.badinput('lon, lat, and profiles must have the same first two dimensions')
        elseif ndims(lon) > 2 && (sz_lon(3) ~= sz_prof(4) || sz_lat(3) ~= sz_prof(4))
            E.badinput('A 3-D lon or lat must have the same size in the 3rd dimension as profiles does in the 4th')
        elseif ndims(alt) == 3 && ~all(sz_alt(1:3) == sz_prof(1:3))
            E.badinput('If passing a 3-D alt array, the dimensions must be the same size as the first three dimensions of the profiles array')
        elseif ndims(alt) == 4 && ndims(profiles) ~= 4
            E.badinput('A 4-D alt array implies multiple days, but the dimensionality of profiles does not')
        elseif ndims(alt) == 4 && ~all(sz_alt(1:4) == sz_prof(1:4))
            E.badinput('If passing a 4-D alt array, all dimensions must be the same size as those of profiles')
        elseif ndims(alt) < 3 || ndims(alt) > 4
            E.badinput('alt must be 3- or 4-D')
        end
        
        if nargin < 7 % no value passed, assume we're going to plot everything on one figure
            multifig = false;
        elseif iscell(multifig) % if a cell is passed, it'll be used in the legend. Make sure it has the right number of elements.
            legstr = multifig;
            if numel(legstr) ~= size(profiles,4) % size(X,dim) doesn't break if dim > ndims, just returns 1
                E.badinput('When passing a cell array as multifig, it must have as many elements as days represented in profiles')
            end
            multifig = false;
        elseif ~isscalar(multifig) || ~islogical(multifig) || ~isnumeric(multifig) % a value passed that's not a cell - is it a logical value
            E.badinput('multifig must be a scalar logical or numeric, a cell array with as many cells as days plotted, or be omitted')
        end
        
        if ~multifig && ~exist('legstr','var') % make the default legend entries if needed
            legstr = cell(1,size(profiles,4));
            for a=1:size(profiles,4)
                legstr{a} = sprintf('Day %d',a);
            end
        end
        
        if nargin < 8
            xtext = '';
        elseif ~ischar(xtext)
            E.badinput('xtext must be a string or omitted')
        end
        
        multiday = size(profiles,4) > 1;
        if multiday && size(lon,3) == 1
            lon = repmat(lon,1,1,sz_prof(4));
            lat = repmat(lat,1,1,sz_prof(4));
        end
        if multiday && size(alt,4) == 1
            alt = repmat(alt,1,1,1,sz_prof(4));
        end
        
        figure;
        for a=1:size(profiles,4)
            [x,y] = find_square_around(lon(:,:,a), lat(:,:,a), target_lon, target_lat, 0);
            this_prof = squeeze(profiles(x,y,:,a));
            this_alt = squeeze(alt(x,y,:,a));
            
            % Make a new figure if plotting on separate figures
            if multifig && a > 1
                figure;
            end
            
            plot(this_prof, this_alt)
            set(gca,'fontsize',14);
            xlabel(xtext)
            ylabel('Altitude (m)')
            
            % If doing multiple figures, use the legend text as the title
            if multifig
                title(legstr{a});
            else
                hold on
            end
        end
        
        if ~multifig
            legend(legstr{:});
        end

    end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% OTHER FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [xx, yy] = find_square_around(lon, lat, center_lon, center_lat, radius)
        % Finds indicies for a square of points centered on center_lon and
        % center_lat. The square will have sides of length 2*radius + 1.
        % Slicing lon/lat as lon(xx,yy), lat(xx,yy) will give the points.
        
        % Check input
        if ~ismatrix(lon)
            E.badinput('lon must be a 2-D array')
        elseif ~ismatrix(lat)
            E.badinput('lat must be a 2-D array')
        elseif ~all(size(lon) == size(lat))
            E.badinput('lon and lat must be the same size')
        elseif ~isscalar(center_lat) || ~isnumeric(center_lat)
            E.badinput('center_lat must be a numeric scalar')
        elseif ~isscalar(center_lon) || ~isnumeric(center_lon)
            E.badinput('center_lon must be a numeric scalar')
        elseif ~isscalar(radius) || ~isnumeric(radius) || radius < 0 || mod(radius,1) ~= 0
            E.badinput('radius must be a positive scalar integer')
        end
        
        del = abs(lon - center_lon) + abs(lat - center_lat);
        [~,I] = min(del(:));
        [x,y] = ind2sub(size(del),I);
        xx = (x-radius):(x+radius);
        yy = (y-radius):(y+radius);
    end

    function M = unstagger(M, dim)
        permvec = 1:ndims(M);
        permvec(1) = dim;
        permvec(dim) = 1;
        M = permute(M,permvec);
        M = (M(1:end-1,:,:) + M(2:end,:,:))/2;
        M = permute(M, permvec);
    end

    function profs = get_nearest_profiles(lon, lat, profiles, target_lon, target_lat)
        % Will return a matrix of the profiles closest to the target_lon
        % and target_lat. Follows the same rules about the relationship
        % between the dimensionality of lon, lat, and profile as the
        % function plot_nearest_profile
        if ndims(lon) < 2 || ndims(lon) > 3 || ndims(lat) < 2 || ndims(lon) > 3
            E.badinput('lon and lat must be 2- or 3-D arrays')
        elseif ndims(profiles) < 3 || ndims(profiles) > 4 
            E.badinput('profiles must be a 3- or 4-D array')
        elseif ndims(profiles) == 3 && (ndims(lon) > 2 || ndims(lat) > 2) %#ok<*ISMAT>
            E.badinput('A 3-D lon or lat array implies multiple days, but the dimensionality of profiles does not')
        elseif ~isscalar(target_lat) || ~isnumeric(target_lat) || ~isscalar(target_lon) || ~isnumeric(target_lon)
            E.badinput('The target_lon and target_lat inputs must be scalar numeric values')
        end
        
        sz_lon = size(lon);
        sz_lat = size(lat);
        sz_prof = size(profiles);
        
        if ndims(lon) ~= ndims(lat)
            E.badinput('lon and lat must have the same number of dimensions')
        elseif ~all(sz_lon == sz_lat)
            E.badinput('lon and lat must be the same size')
        elseif ~all(sz_lon(1:2) == sz_prof(1:2))
            E.badinput('lon, lat, and profiles must have the same first two dimensions')
        elseif ndims(lon) > 2 && (sz_lon(3) ~= sz_prof(4) || sz_lat(3) ~= sz_prof(4))
            E.badinput('A 3-D lon or lat must have the same size in the 3rd dimension as profiles does in the 4th')
        end
        
        if size(profiles,4) > 1 && size(lon,3) == 1
            lon = repmat(lon,1,1,sz_prof(4));
            lat = repmat(lat,1,1,sz_prof(4));
        end
        
        profs = nan(size(profiles,3), size(profiles,4));
        for a=1:size(profiles,4)
            [x,y] = find_square_around(lon(:,:,a), lat(:,:,a), target_lon, target_lat, 0);
            profs(:,a) = squeeze(profiles(x,y,:,a));
        end
    end

end

