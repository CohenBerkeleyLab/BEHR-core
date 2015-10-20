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
    case 'plotapriori'
        plot_apriori(varargin{:});
    case 'getnearprof'
        varargout{1} = get_nearest_profiles(varargin{:});
    case 'calcwind'
        varargout{1} = calc_wind_mag(varargin{1}, varargin{2});
        varargout{2} = calc_wind_dir(varargin{1}, varargin{2});
    case 'perdiffvstheta'
        plot_delamf_vs_delangle(varargin{:});
    case 'perrec'
        plot_perrec_vs_distance(varargin{:});
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
            theta = nan(size(U_a));
            for b=1:numel(U_a)
                if U_a(b) >= 0
                    theta(b) = atand(V_a(b)/U_a(b));
                else
                    theta(b) = atand(V_a(b)/U_a(b))+180;
                end
            end
        
            mag_mean(a) = nanmean(mag(:));
            if radius > 0
                mag_std(a) = nanstd(mag(:));
            end
            
            % If all the points are between -90 and 90, we want to average
            % before taking mod(360) b/c that will properly handle points
            % around 0: mean([1,-1]) = 0 vs. mean([1,359]) = 180
            % Otherwise it's okay to take the mod(360) first then average,
            % and in fact we must to correctly handle points around 180:
            % mean([179,181]) = 180 but mean([179,-179]) = 0. Similar
            % logic applies to the std. deviation.
            if all(theta(:)<90) && all(theta(:)>-90)
                theta_mean(a) = mod(nanmean(theta(:)),360);
                if radius > 0
                    theta_std(a) = nanstd(theta(:));
                end
            else
                theta_mean(a) = nanmean(mod(theta(:),360));
                if radius > 0
                    theta_std(a) = nanstd(mod(theta(:),360));
                end
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

    function plot_apriori(DataHourly, DataMonthly, indicies, shape_factor)
        % Function that will plot the a priori profiles from a given
        % satellite pixel to compare the daily and monthly profiles. Takes
        % two Data structures output from BEHR_Main and an n-by-2 matrix of
        % indicies to plot. These are the indicies of the pixels in the
        % matrices in Data. Each row should correspond to a different
        % pixel, and will be plotted in it's own figure.  The profiles will
        % be normalized by their VCDs so that the shape factor is plotted,
        % unless a 0 is given as the optional fourth argument.
        
        % Input checking
        req_fields = {'BEHRNO2apriori','BEHRPressureLevels','GLOBETerpres'};
        if ~isstruct(DataHourly) || ~isstruct(DataMonthly) ||... %continued next line
                any(~isfield(DataHourly, req_fields)) || any(~isfield(DataMonthly, req_fields))
            E.badinput('DataHourly and DataMonthly must be Data structures output from BEHR_Main (must contain the fields %s)',strjoin(req_fields,', '));
        end
        
        if ~isnumeric(indicies) || size(indicies, 2) ~= 2
            E.badinput('indicies must be an n-by-2 matrix of subscript indicies')
        end
        
        if nargin < 4
            shape_factor = 1;
        elseif (~isnumeric(shape_factor) && ~islogical(shape_factor)) || ~isscalar(shape_factor)
            E.badinput('shape_factor must be scalar numeric or logical, or be omitted')
        end
        
        % Plotting
        for a=1:size(indicies,1)
            x = indicies(a,1);
            y = indicies(a,2);
            
            apriori_hr = DataHourly.BEHRNO2apriori(:,x,y);
            pres_hr = DataHourly.BEHRPressureLevels(:,x,y);
            surfP_hr = DataHourly.GLOBETerpres(x,y);
            if shape_factor
                vcd_hr = integPr2(apriori_hr, pres_hr, surfP_hr);
            else 
                vcd_hr = 1;
            end
            
            apriori_mn = DataMonthly.BEHRNO2apriori(:,x,y);
            pres_mn = DataMonthly.BEHRPressureLevels(:,x,y);
            surfP_mn = DataMonthly.GLOBETerpres(x,y);
            if shape_factor
                vcd_mn = integPr2(apriori_mn, pres_mn, surfP_mn);
            else
                vcd_mn = 1;
            end
            
            figure; 
            plot(apriori_hr/vcd_hr, pres_hr);
            hold on
            plot(apriori_mn/vcd_mn, pres_mn);
            
            set(gca,'ydir','reverse');
            set(gca,'fontsize',14);
            if shape_factor
                xlabel('Shape factor')
            else
                xlabel('[NO_2]')
            end
            ylabel('Pressure (hPa)')
            legend('Hourly','Monthly');
            if shape_factor
                title_sf = 'shape factor';
            else
                title_sf = 'conc';
            end
            title(sprintf('Pixel at %.2f W, %.2f N (%s) - %s',DataHourly.Longitude(x,y),DataHourly.Latitude(x,y),mat2str([x,y]), title_sf));
        end
    end

    function varargout = plot_delamf_vs_delangle(Data_D, Data_M, wind_angle, city_lon, city_lat, wind_speed)
        % This function will generate a scatter plot of percent difference
        % in BEHR AMF values between a any two sets of pixels as a function
        % of position of the pixel relative to the city compared to wind
        % direction. Arguments:
        %   Data_D: data structure with the new AMFs. Must contain
        %   Longitude, Latitude, and BEHRAMFTrop as fields.
        %
        %   Data_M: data structure with the old AMFs. Must contain the same
        %   fields.
        %
        %   wind_angle: angle of the wind around the city in degrees CCW
        %   from east, must be between 0 and 360.
        %
        %   city_lon, city_lat: the longitude/latitude coordinates of the
        %   city, using the convention that longitude west of 0 is
        %   negative.
        %
        %   wind_speed: (optional) a wind speed in km/s, used to normalize
        %   pixel distanced to wind speed over the city.  More useful in
        %   conjunction with another function that plots multiple days with
        %   different wind speed and direction on one plot.
        
        % INPUT CHECKING:
        req_fields = {'Longitude','Latitude','BEHRAMFTrop'};
        if ~isstruct(Data_D) || any(~isfield(Data_D,req_fields))
            E.badinput('Data_D must be a structure with fields %s',strjoin(req_fields,', '));
        elseif ~isstruct(Data_M) || any(~isfield(Data_M, req_fields))
            E.badinput('Data_M must be a structure with fields %s',strjoin(req_fields,', '));
        end
        
        % Also check that the pixels are the same in both structures
        if ndims(Data_D.Longitude) ~= ndims(Data_M.Longitude) || any(size(Data_D.Longitude) ~= size(Data_M.Longitude))
            E.sizeMismatch('Data_D.Longitude','Data_M.Longitude');
        elseif any(abs(Data_D.Longitude(:) - Data_M.Longitude(:))>0.01) || any(abs(Data_D.Latitude(:) - Data_M.Latitude(:))>0.01)
            E.badinput('Data_M and Data_D seem to have different pixels (based on lat/lon coordinates');
        end
        
        if ~isscalar(wind_angle) || ~isnumeric(wind_angle) || wind_angle < 0 || wind_angle > 360
            E.badinput('wind_angle must be a scalar number between 0 and 360');
        elseif ~isscalar(city_lon) || ~isnumeric(city_lon) || city_lon < -180 || city_lon > 180
            E.badinput('city_lon must be a scalar number between -180 and 180');
        elseif ~isscalar(city_lat) || ~isnumeric(city_lat) || city_lat < -90 || city_lat > 90
            E.badinput('city_lat must be a scalar number between -90 and 90');
        end
        
        % cbtitle will be used later to label the colorbar, whether it is
        % simply distance to city or is in fact distance normalized by wind
        % speed (in case of multiple days)
        if ~exist('wind_speed','var')
            wind_speed = 1;
            cbtitle = 'Distance to city (km)';
        else
            if ~isscalar(wind_speed) || ~isnumeric(wind_speed) || wind_speed < 0
                E.badinput('wind_speed should be a positive scalar number');
            end
            cbtitle = 'Distance to city normalized by wind speed';
        end
        
        % MAIN FUNCTION:
        % Calculate 3 quantities: the percent difference in AMF from the
        % new and old a priori, the angle of each pixel to the city (in deg
        % CCW from east), and the distance of each pixel to the city
        
        % Per. diff. AMF:
        perdiff = (Data_D.BEHRAMFTrop(:) ./ Data_M.BEHRAMFTrop(:) - 1) * 100;
        
        % Difference between angle to pixel and angle of wind
        city_lon_mat = repmat(city_lon, size(Data_D.Longitude));
        city_lat_mat = repmat(city_lat, size(Data_D.Latitude));
        
        theta_pix = latlon_angle(city_lon_mat, city_lat_mat, Data_D.Longitude, Data_D.Latitude);
        
        % Distance from city to pixels.
        dist_to_city = nan(size(Data_D.Longitude));
        for a=1:numel(dist_to_city)
            dist_to_city(a) = m_lldist([city_lon, Data_D.Longitude(a)], [city_lat, Data_D.Latitude(a)]);
        end
        
        % Figure creation. 36 is the default marker size for a scatter plot
        % in Matlab 2014b.
        figure; 
        scatter(theta_pix(:) - wind_angle, perdiff(:), 36, dist_to_city(:) ./ wind_speed);
        cb = colorbar;
        set(gca,'fontsize',16)
        xlabel('\theta_{pix} - \theta_{wind}')
        ylabel('%\Delta AMF')
        cb.Label.String = cbtitle;
        
    end

    function plot_perrec_vs_distance(zero_lon, zero_lat, Monthly, Daily, Hybrid, plot_mode)
        % This function will plot the percent change recovered in the
        % hybrid a priori versus the standard daily a priori, versus
        % distance from the city center.  This is intended to show that
        % near a city, the hybrid profile (which only changes near surface
        % NO2 with the daily profile, and uses the monthly profile for free
        % troposphere stuff) recoveres a good percent of the change, and
        % further away the change is more dependent on FT variability.
        %
        % Needs 5 required inputs, one optional:
        %  zero_lon, zero_lat - the lon/lat distance is calculated from
        %  Monthly, Daily, Hybrid - the respective structures containing
        %   the Data structure from BEHR output. These can either be the
        %   Data structures themselves (for one day) or a structure
        %   containing multiple Data structures for different days, as the
        %   field Data.  Monthly can be either a single structure that will
        %   be applied for all Daily & Hybrid structures, or have one
        %   structure for each Daily or Hybrid structure.
        %  plot_mode - (optional) determines how multiple days are plotted,
        %   should be entered as one of the strings here:
        %       'allinone' or '' - plots all pixels on the same plot
        %          without any identification by day.
        %       'colored' - will color the points by day.
        %       'separate' - will plot each day on its own figure.
        
        %%%%% INPUT PARSING %%%%%
        if ~isnumeric(zero_lon) || ~isscalar(zero_lon)
            E.badinput('zero_lon must be a numeric scalar')
        elseif ~isnumeric(zero_lat) || ~isscalar(zero_lat)
            E.badinput('zero_lon must be a numeric scalar')
        end
        
        if ~isstruct(Monthly) || (~isfield(Monthly,'Data') && ~isfield(Monthly,'BEHRAMFTrop'))
            E.badinput('Monthly must be a BEHR Data structure, or a structure containing multiple Data structures, as Monthly(1).Data, Monthly(2).Data, etc.')
        elseif ~isstruct(Daily) || (~isfield(Daily,'Data') && ~isfield(Daily,'BEHRAMFTrop'))
            E.badinput('Daily must be a BEHR Data structure, or a structure containing multiple Data structures, as Daily(1).Data, Daily(2).Data, etc.')
        elseif ~isstruct(Hybrid) || (~isfield(Hybrid,'Data') && ~isfield(Hybrid,'BEHRAMFTrop'))
            E.badinput('Hybrid must be a BEHR Data structure, or a structure containing multiple Data structures, as Hybrid(1).Data, Hybrid(2).Data, etc.')
        end
        
        if isfield(Monthly,'BEHRAMFTrop')
            Monthly.Data = Monthly;
        end
        if isfield(Daily,'BEHRAMFTrop')
            Daily.Data = Daily;
        end
        if isfield(Hybrid,'BEHRAMFTrop')
            Hybrid.Data = Hybrid;
        end
        
        n = numel(Daily);
        if numel(Hybrid) ~= n
            E.badinput('The Daily and Hybrid structures must have the same number of entries')
        elseif numel(Monthly) ~= 1 && numel(Monthly) ~= n
            E.badinput('The Monthly structure must have the same number of entries as the Daily and Hybrid structures, or have only one entry');
        end
        
        if numel(Monthly) == 1 && n ~= 1
            Monthly = repmat(Monthly, size(Daily));
        end
        
        if ~exist('plot_mode','var')
            plot_mode = 'allinone';
        else
            plot_mode = lower(plot_mode);
            allowed_modes = {'allinone','colored','separate'};
            if any(~ismember(plot_mode, allowed_modes))
                E.badinput('%s is not a valid value for plot_mode; it must be one of %s',plot_mode,strjoin(allowed_modes,', '));
            end
        end
        
        %%%%% MAIN FUNCTION %%%%%
        perrec = [];
        distance_to_zero = [];
        day_index = [];
        for a=1:n;
            lon_check = any(Hybrid(a).Data.Longitude(:) ~= Monthly(a).Data.Longitude(:)) || any(Daily(a).Data.Longitude(:) ~= Monthly(a).Data.Longitude(:));
            lat_check = any(Hybrid(a).Data.Latitude(:) ~= Monthly(a).Data.Latitude(:)) || any(Daily(a).Data.Latitude(:) ~= Monthly(a).Data.Latitude(:));
            if lon_check || lat_check
                E.badinput('The lat/lon coordinates do not match among index %d of the three structures', a);
            end
                
            this_perrec = (Hybrid(a).Data.BEHRAMFTrop(:) - Monthly(a).Data.BEHRAMFTrop(:)) ./ (Daily(a).Data.BEHRAMFTrop(:) - Monthly(a).Data.BEHRAMFTrop(:)) * 100;
            this_distance = nan(numel(Monthly(a).Data.Longitude),1);
            for b=1:numel(Monthly(a).Data.Longitude)
                this_distance(b) = m_lldist([Monthly(a).Data.Longitude(b), zero_lon], [Monthly(a).Data.Latitude(b), zero_lat]);
            end
            this_day = ones(size(this_perrec)) .* a; % used only for coloring by day
            
            if strcmp(plot_mode, 'separate')
                figure;
                scatter(this_distance, this_perrec);
                set(gca,'fontsize',16)
                xlabel('Distance from city (km)')
                ylabel('Percent difference recovered');
                title(sprintf('Day %d',a));
            else
                perrec = cat(1, perrec, this_perrec);
                distance_to_zero = cat(1, distance_to_zero, this_distance);
                day_index = cat(1, day_index, this_day);
            end
        end
        
        if ~strcmp(plot_mode, 'separate')
            figure;
            if strcmp(plot_mode, 'allinone')
                scatter(distance_to_zero, perrec);
            elseif strcmp(plot_mode, 'colored');
                scatter(distance_to_zero, perrec, 36, day_index);
                cb = colorbar;
                cb.Label.String = 'Day index';
                cb.Label.FontSize = 16;
            end
            set(gca,'fontsize',16)
            xlabel('Distance from city (km)')
            ylabel('Percent difference recovered');
        end

    end


%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% OTHER FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%`

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

    function mag = calc_wind_mag(U,V)
        mag = (U.^2 + V.^2).^0.5;
    end

    function theta = calc_wind_dir(U,V)
        theta = nan(size(U));
        for b=1:numel(U)
            if U(b) >= 0
                theta(b) = atand(V(b)/U(b));
            else
                theta(b) = atand(V(b)/U(b))+180;
            end
        end
    end

end

