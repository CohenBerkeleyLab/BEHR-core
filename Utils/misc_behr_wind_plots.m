function [ varargout ] = misc_behr_wind_plots( plttype, varargin )
%MISC_BEHR_WIND_PLOTS Various plots studying the effect of wind variation on BEHR
%   Like the other misc. plot functions, this collects several related
%   plotting functions in one place.

E = JLLErrors;

C=load('blue_red_cmap');

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
    case 'sectors'
        plot_changes_by_sector(varargin{:});
    case 'dif'
        plot_diff();
    case 'res'
        plot_diff_resolutions();
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
            lon = repmat(lon,1,1,size(U,3));
            lat = repmat(lat,1,1,size(U,3));
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

    function plot_apriori(DataHourly, DataMonthly, DataHybrid, indicies, shape_factor)
        % Function that will plot the a priori profiles from a given
        % satellite pixel to compare the daily and monthly profiles. Takes
        % two Data structures output from BEHR_Main and an n-by-2 matrix of
        % indicies to plot. These are the indicies of the pixels in the
        % matrices in Data. Each row should correspond to a different
        % pixel, and will be plotted in it's own figure.  The profiles will
        % be normalized by their VCDs so that the shape factor is plotted,
        % unless a 0 is given as the optional fourth argument. A 2 as the
        % fourth argument will plot both.
        
        % Input checking
        req_fields = {'BEHRNO2apriori','BEHRPressureLevels','GLOBETerpres'};
        if ~isstruct(DataHourly) || ~isstruct(DataMonthly) ||... %continued next line
                any(~isfield(DataHourly, req_fields)) || any(~isfield(DataMonthly, req_fields))
            E.badinput('DataHourly and DataMonthly must be Data structures output from BEHR_Main (must contain the fields %s)',strjoin(req_fields,', '));
        end
        
        % Handle if the user passes a hybrid profile as well
        if isstruct(DataHybrid)
            if any(~isfield(DataHybrid, req_fields))
                E.badinput('DataHybrid, if given, must be Data structures output from BEHR_Main (must contain the fields %s)',strjoin(req_fields,', '));
            end
            use_hybrid = true;
            
            if nargin < 4
                shape_factor = 1;
            elseif (~isnumeric(shape_factor) && ~islogical(shape_factor)) || ~isscalar(shape_factor)
                E.badinput('shape_factor must be scalar numeric or logical, or be omitted')
            end
        else
            if nargin >= 4
                shape_factor = indicies;
            else
                shape_factor = 1;
            end
            indicies = DataHybrid;
            
            use_hybrid = false;
        end
        
        if ~isnumeric(indicies) || size(indicies, 2) ~= 2
            E.badinput('indicies must be an n-by-2 matrix of subscript indicies')
        end
        
        
        
        % Plotting
        for a=1:size(indicies,1)
            x = indicies(a,1);
            y = indicies(a,2);
            
            apriori_hr = DataHourly.BEHRNO2apriori(:,x,y);
            pres_hr = DataHourly.BEHRPressureLevels(:,x,y);
            surfP_hr = DataHourly.GLOBETerpres(x,y);
            if shape_factor > 0
                vcd_hr = integPr2(apriori_hr, pres_hr, surfP_hr);
            else 
                vcd_hr = 1;
            end
            
            apriori_mn = DataMonthly.BEHRNO2apriori(:,x,y);
            pres_mn = DataMonthly.BEHRPressureLevels(:,x,y);
            surfP_mn = DataMonthly.GLOBETerpres(x,y);
            if shape_factor
                vcd_mn = integPr2(apriori_mn, pres_mn, surfP_mn);
                aesthetic_scale = max(max(apriori_hr/vcd_hr),max(apriori_mn/vcd_mn));
            else
                vcd_mn = 1;
                aesthetic_scale = 1;
            end
            
            if use_hybrid
                apriori_hy = DataHybrid.BEHRNO2apriori(:,x,y);
                pres_hy = DataHybrid.BEHRPressureLevels(:,x,y);
                surfP_hy = DataHybrid.GLOBETerpres(x,y);
                if shape_factor
                    vcd_hy = integPr2(apriori_hy, pres_hy, surfP_hy);
                    aesthetic_scale = max(aesthetic_scale, max(apriori_hy/vcd_hy));
                else
                    vcd_hy = 1;
                end
            end
            
            % This value will be used if only plotting shape factor to
            % scale it such that the x-axis doesn't need a power of 10
            if shape_factor
                aesthetic_scale = 10^(-floor(log10(aesthetic_scale)));
            end
            
            figure; 
            
            if shape_factor > 1
                line(apriori_hr, pres_hr, 'color','r','linewidth',2);
                line(apriori_mn, pres_mn, 'color','b','linewidth',2);
                if use_hybrid
                    line(apriori_hy, pres_hy, 'color',[0 0.5 0],'linewidth',2);
                end
                line(apriori_hr/vcd_hr*1e15, pres_hr, 'color','r','linewidth',2,'linestyle','--');
                line(apriori_mn/vcd_mn*1e15, pres_mn, 'color','b','linewidth',2,'linestyle','--');
                if use_hybrid
                    line(apriori_hy/vcd_hy*1e15, pres_hy, 'color',[0 0.5 0],'linewidth',2,'linestyle','--');
                end
            else
                line(apriori_hr/vcd_hr*aesthetic_scale, pres_hr, 'color','r','linewidth',2);
                line(apriori_mn/vcd_mn*aesthetic_scale, pres_mn, 'color','b','linewidth',2);
                if use_hybrid
                    line(apriori_hy/vcd_hy*aesthetic_scale, pres_hy, 'color',[0 0.5 0],'linewidth',2);
                end
            end
            set(gca,'ydir','reverse');
            set(gca,'fontsize',14);
            if shape_factor == 1
                xlabel('Shape factor')
            elseif shape_factor == 2
                xlabel(sprintf('[NO_2] (mixing ratio) or\nshape factor (unitless, scaled)'));
            else
                xlabel('[NO_2]')
            end
            ylabel('Pressure (hPa)')
            if shape_factor > 1 && ~use_hybrid
                legend('Daily (mixing ratio)','Monthly (mixing ratio)','Daily (shape factor)','Monthly (shape factor)');
            elseif shape_factor > 1 && use_hybrid
                legend('Daily (mixing ratio)','Monthly (mixing ratio)','Hybrid (mixing ratio)','Daily (shape factor)','Monthly (shape factor)','Hybrid (shape factor)');
            elseif shape_factor <= 1 && ~use_hybrid
                legend('Daily','Monthly');
            else
                legend('Daily','Monthly','Hybrid');
            end
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

    function plot_diff()
        % This function uses no input, instead it will ask a series of
        % questions to decide what plots to make.  This is used to make
        % difference pcolor plots of WRF-Chem VCDs or BEHR VCDs or AMFs.
        
        %%%%%%%%%%%%%%%%%
        %%%%% INPUT %%%%%
        %%%%%%%%%%%%%%%%%
        
        % First question: WRF or BEHR. pseudo-BEHR is the one where I used
        % the same set of pixels each day.
        allowed_sources = {'wrf','behr','pseudo-behr'};
        q_str = sprintf('Which source will you be using? %s: ',strjoin(allowed_sources, ', '));
        while true
            source = lower(input(q_str,'s'));
            if ~ismember(source, allowed_sources)
                fprintf('You must select one of the allowed choices. Try again, or press Ctrl+C to cancel\n');
            else
                break
            end
        end
        
        % Follow up if using BEHR: should we compare hybrid or hour-wise
        % data to the monthly average?
        allowed_apriori = {'hourly','hybrid'};
        while true
            q_str = 'Compare retrieval using the hourly or hybrid data to the monthly average? ';
            apriori = lower(input(q_str,'s'));
            if ~ismember(apriori, allowed_apriori)
                fprintf('You must enter one of the allowed choices: %s. Try again, or press Ctrl+C to cancel\n',strjoin(allowed_apriori,', '));
            else
                break
            end
        end
        
        % Second question: which city. Expandable.
        allowed_cities = {'Atlanta'};
        n_cities = numel(allowed_cities);
        while true
            fprintf('Which city should we focus on?\n')
            for a=1:n_cities
                fprintf('\t%d - %s\n', a, allowed_cities{a});
            end
            city_index = str2double(input('Enter the number of your selection: ','s'));
            if isnan(city_index) || city_index < 1 || city_index > n_cities
                fprintf('You must choose a number between 1 and %d, or press Ctrl+C to cancel\n', n_cities);
            else
                break
            end
        end
        
        % Third question: VCDs or AMFs. If WRF, VCDs are the only option
        allowed_quantities = {'amf','vcd'};
        q_str = sprintf('Which source will you be using? %s: ',strjoin(allowed_quantities, ', '));
        if strcmpi(source,'wrf')
            quantity = 'vcd';
        else
            while true
                quantity = lower(input(q_str,'s'));
                if ~ismember(quantity, allowed_quantities)
                    fprintf('You must select one of the allowed choices. Try again, or press Ctrl+C to cancel\n');
                else
                    break
                end
            end
        end
        
        % Follow up only if doing VCDs: which ghost column to use? 
        if ~isempty(regexpi(source,'behr'))
            allowed_ghosts = {'none','new','old'};
            q_str = sprintf('Which ghost column correction do you want to use? %s: ', strjoin(allowed_ghosts, ', '));
            while true 
                ghost = lower(input(q_str, 's'));
                if ~ismember(ghost, allowed_ghosts)
                    fprintf('You must select one of the allowed choices. Try again, or press Ctrl+C to cancel\n');
                else
                    break
                end
            end
        else
            ghost = 'old';
        end
        
        % Fourth, is this a percent or absolute difference?
        while true
            diff_type = lower(input('Do you want absolute or percent differences, or just the daily or monthly values? Type p, a, d, or m: ','s'));
            if ~ismember(diff_type,{'a','p','d','m'});
                fprintf('You must select one of the allowed choices. Try again, or press Ctrl+C to cancel\n');
            else 
                break
            end
        end
        
        % Lastly, we need a date
        while true
            date_in = input('Enter the date to compare, using a format datenum can parse: ', 's');
            try
                datenum_in = datenum(date_in);
                break % only executes if the previous line is successfull
            catch err
                if strcmp(err.identifier, 'MATLAB:datenum:ConvertDateString')
                    fprintf('The format could not be recognized by datenum. Try yyyy-mm-dd, or Ctrl+C to cancel\n')
                else
                    rethrow(err)
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% PARSING INPUT %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        % File paths: need the daily and monthly paths. Will construct a
        % cell array where the first dimension corresponds to the
        % difference sources (e.g. WRF, BEHR, pseudo-BEHR) and the second
        % to cities. Put NaNs in cells that don't have data. The order
        % needs to match the "allowed_sources" and "allowed_cities"
        % variables.
        
        switch ghost
            case 'old'
                daily_path = {  '/Volumes/share2/USERS/LaughnerJ/WRF/SE_US_BEHR/NEI11Emis/hourly';...
                    '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hourly';...
                    '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hourly - No clouds'};
                % this will be concatenated in the third dimensions with daily
                % path, so it needs to be the same size. Fill with NaNs if needed.
                hybrid_path = { NaN;...
                    '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hybrid';...
                    '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hybrid - No clouds'};
                
                monthly_path = {  '/Volumes/share2/USERS/LaughnerJ/WRF/SE_US_BEHR/NEI11Emis/monthly';...
                    '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Monthly';...
                    '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Monthly - No clouds'};
            case 'none' 
                daily_path = {  NaN;...
                                NaN;...
                                '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hourly - No clouds - No ghost'};
                hybrid_path = { NaN;...
                                '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hybrid - No ghost';...
                                '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hybrid - No clouds - No ghost'};
                monthly_path = {NaN;...
                                '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Monthly - No ghost';...
                                '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Monthly - No clouds - No ghost'};
            case 'new'
                daily_path = {  NaN;...
                                NaN;...
                                '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hourly - No clouds - New ghost'};
                hybrid_path = { NaN;...
                                '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hybrid - New ghost';...
                                '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hybrid - No clouds - New ghost'};
                monthly_path = {NaN;...
                                '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Monthly - New ghost';...
                                '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Monthly - No clouds - New ghost'};
        end
        
                    
        daily_path = cat(3, daily_path, hybrid_path);
                    
        source_index = find(strcmp(source, allowed_sources));
        if isempty(regexp(source,'behr','once'))
           apriori_index = 1; 
        else
           apriori_index = find(strcmp(apriori, allowed_apriori));
        end
        % city_index is calculated during the user response
        
        daily_path = daily_path{source_index, city_index, apriori_index}; 
        monthly_path = monthly_path{source_index, city_index};
        
        % File name format: strings should include a %1$04d where the year,
        % %2$02d where the month, and %3$02d where the day get filled in. Only
        % expected to be different for the different sources, again, follow
        % the order of the allowed_sources variable.
        
        daily_file_name_spec = {    'WRF_BEHR_hourly_%1$04d-%2$02d-%3$02d.nc';...
                                    'OMI_BEHR_%1$04d%2$02d%3$02d.mat';...
                                    'OMI_BEHR_%1$04d%2$02d%3$02d.mat'};
        monthly_file_name_spec = {    'WRF_BEHR_monthly_%1$04d-%2$02d-%3$02d.nc';...
                                    'OMI_BEHR_%1$04d%2$02d%3$02d.mat';...
                                    'OMI_BEHR_%1$04d%2$02d%3$02d.mat'};
                                
        daily_file_name = sprintf(daily_file_name_spec{source_index}, year(date_in), month(date_in), day(date_in));
        if strcmp(source,'wrf')
            monthly_file_name = sprintf(monthly_file_name_spec{source_index}, year(date_in), month(date_in), eomday(year(date_in), month(date_in)));
        else
            monthly_file_name = sprintf(monthly_file_name_spec{source_index}, year(date_in), month(date_in), day(date_in));
        end
        
        % Which swath or hour to use (BEHR or WRF respectively). First
        % dimension is by source, second by city.
        
        swath_or_hour = [   2;
                            2;
                            2];
                        
        s = swath_or_hour(source_index, city_index);
        
        % Cutting down the arrays to the right pixels. Again, source down
        % the first dimension, city across the second.
        
        xxes = {    148:185;...
                    761:860;...
                    1:11}; % use all the pixels in the pseudo retrieval
        yyes = {    65:97;...
                    141:210;...
                    1:19}; % use all the pixels in the pseudo retrieval
        
        xx = xxes{source_index, city_index};
        yy = yyes{source_index, city_index};
        
        % Plot x and y limits. Only varies by city.
        city_xlims = {[-87.1 -82]};
        city_ylims = {[31.9 35.5]};
        
        xl = city_xlims{city_index};
        yl = city_ylims{city_index};
        
        % City lon and lat center. Obviously only varies by city.
        all_city_lons = [-84.39];
        all_city_lats = [33.775];
        all_city_names = {'Atlanta'};
        
        city_lon = all_city_lons(city_index);
        city_lat = all_city_lats(city_index);
        city_name = all_city_names{city_index};
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% LOAD DATA, CALCULATE QUANTITIES, AND PLOT %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        daily_file = fullfile(daily_path, daily_file_name);
        monthly_file = fullfile(monthly_path, monthly_file_name);
        
        if strcmp(source,'wrf') && strcmp(quantity,'amf')
            E.notimplemented('%s','WRF output does not contain an AMF value');
        elseif strcmp(source,'wrf') && strcmp(quantity,'vcd')
            
            
            XLONG = ncread(daily_file, 'XLONG');
            xlon = XLONG(xx,yy,1);
            XLAT = ncread(daily_file, 'XLAT');
            xlat = XLAT(xx,yy,1);
            
            daily_no2 = ncread(daily_file, 'no2_ndens'); %[NO2 in number density]
            daily_no2 = daily_no2(xx,yy,:,s); % cut down to the hour we want
            daily_zlev = ncread(daily_file, 'zlev'); % Thickness of each layer in meters
            daily_zlev = daily_zlev(xx,yy,:,s);
            daily_tplev = find_wrf_tropopause(ncinfo(daily_file));
            daily_tplev = daily_tplev(xx,yy,s);
            
            for a=1:size(daily_no2,1)
                for b=1:size(daily_no2,2)
                    tp = daily_tplev(a,b);
                    if tp > 0 % tp is given -1 if the tropopause algorithm cannot find a tropopause
                        daily_no2(a,b,tp:end) = nan;
                    end
                end
            end
            
            daily_no2_columns = nansum2(daily_no2 .* (daily_zlev*100), 3);
            
            monthly_no2 = ncread(monthly_file, 'no2_ndens'); %[NO2 in number density]
            monthly_no2 = monthly_no2(xx,yy,:); % cut down to the hour we want
            monthly_zlev = ncread(monthly_file, 'zlev'); % Thickness of each layer in meters
            monthly_zlev = monthly_zlev(xx,yy,:);
            monthly_tplev = find_wrf_tropopause(ncinfo(monthly_file));
            monthly_tplev = monthly_tplev(xx,yy);
            
            for a=1:size(monthly_no2,1)
                for b=1:size(monthly_no2,2)
                    tp = monthly_tplev(a,b);
                    if tp > 0 % tp is given -1 if the tropopause algorithm cannot find a tropopause
                        monthly_no2(a,b,tp:end) = nan;
                    end
                end
            end
            
            monthly_no2_columns = nansum2(monthly_no2 .* (monthly_zlev*100), 3);
            
            switch diff_type
                case 'a'
                    del = daily_no2_columns - monthly_no2_columns;
                case 'p'
                    del = (daily_no2_columns ./ monthly_no2_columns - 1)*100;
                case 'd'
                    del = daily_no2_columns;
                case 'm'
                    del = monthly_no2_columns;
            end
            
            figure; pcolor(xlon, xlat, del);
            cb=colorbar;
            set(gca,'fontsize',16);
            xlim(xl);
            ylim(yl);
            colormap('jet');
%             if strcmp(diff_type,'a')
%                 cb.Label.String = '\Delta VCD_{NO_2} (molec. cm^{-2})';
%             else
%                 cb.Label.String = '%\Delta VCD_{NO_2}';
%             end
%             l=line(city_lon, city_lat, 'linestyle','none', 'marker','p','markersize',18,'color','k','linewidth',2);
%             legend(l, city_name);
%             if strcmp(diff_type,'a')
%                 cm = max(abs(del(:)));
%                 p10 = (10^floor(log10(cm)));
%                 cm = round(cm/p10)*p10;
%                 caxis([-cm cm]);
%             elseif strcmp(diff_type,'p')
%                 cm = max(abs(del(:)));
%                 cm = round(cm/10)*10;
%                 caxis([-cm cm]);
%             else
%                 cm = max(abs(del(:)));
%                 p10 = (10^floor(log10(cm)));
%                 cm = round(cm/p10)*p10;
%                 caxis([0 cm]);
%             end
            
        elseif ~isempty(regexp(source, 'behr', 'once'))
            if strcmp(source,'behr')
                D = load(daily_file,'OMI');
                M = load(monthly_file,'OMI');
                
                D.OMI = omi_pixel_reject(D.OMI(s),'omi',0.2,'XTrackFlags');
                badpix = D.OMI.Areaweight == 0;
                M.OMI = omi_pixel_reject(M.OMI(s),'omi',0.2,'XTrackFlags');
                
                lon = D.OMI.Longitude(yy,xx);
                lat = D.OMI.Latitude(yy,xx);
                
                if strcmp(quantity, 'amf')
                    D.OMI.BEHRAMFTrop(badpix) = nan;
                    daily_value = D.OMI.BEHRAMFTrop(yy,xx);
                    M.OMI.BEHRAMFTrop(badpix) = nan;
                    monthly_value = M.OMI.BEHRAMFTrop(yy,xx);
                else
                    D.OMI.BEHRColumnAmountNO2Trop(badpix) = nan;
                    daily_value = D.OMI.BEHRColumnAmountNO2Trop(yy,xx);
                    M.OMI.BEHRColumnAmountNO2Trop(badpix) = nan;
                    monthly_value = M.OMI.BEHRColumnAmountNO2Trop(yy,xx);
                end
            elseif strcmp(source,'pseudo-behr')
                D = load(daily_file, 'Data');
                M = load(monthly_file, 'Data');
                
                lon = squeeze(D.Data.Loncorn(1,:,:));
                lat = squeeze(D.Data.Latcorn(1,:,:));
                
                if strcmp(quantity,'amf')
                    daily_value = D.Data.BEHRAMFTrop;
                    monthly_value = M.Data.BEHRAMFTrop;
                else
                    daily_value = D.Data.BEHRColumnAmountNO2Trop;
                    monthly_value = M.Data.BEHRColumnAmountNO2Trop;
                end
            else
                E.notimplemented('The source %s is not understood',source);
            end
            
            switch diff_type
                case 'a'
                    del = daily_value - monthly_value;
                case 'p'
                    del = (daily_value ./ monthly_value - 1) * 100;
                case 'd'
                    del = daily_value;
                case 'm'
                    del = monthly_value;
            end
            
            figure; pcolor(lon, lat, del);
            cb=colorbar;
            set(gca,'fontsize',16);
            xlim(xl);
            ylim(yl);
            if ~isempty(regexpi(source,'pseudo'))
                colormap(C.blue_red_cmap);
            else
                colormap('jet');
                shading flat
            end
            
        end
        
        switch quantity
            case 'vcd'
                label_pt2 = 'VCD_{NO_2}';
                label_unit = '(molec. cm^{-2})';
            case 'amf'
                label_pt2 = 'AMF';
                label_unit = '';
        end
        switch diff_type
            case 'a'
                label_pt1 = '\Delta';
            case 'p'
                label_pt1 = '%\Delta';
                label_unit = ''; % remove the unit if doing a percent difference
            case 'd'
                label_pt1 = 'Daily';
            case 'm'
                label_pt1 = 'Monthly';
        end
        
        cb.Label.String = strjoin({label_pt1, label_pt2, label_unit}, ' ');
        
        l=line(city_lon, city_lat, 'linestyle','none', 'marker','p','markersize',18,'color','k','linewidth',2);
        legend(l,city_name);
        cm = max(abs(del(:)));
        if ~isnan(cm)
            if strcmp(diff_type,'p')
                cm = round(cm/10)*10;
                caxis([-cm cm]);
            else
                p10 = (10^floor(log10(cm)));
                cm = round(cm/p10)*p10;
                if strcmp(diff_type,'a')
                    caxis([-cm cm]);
                else
                    caxis([0 cm]);
                end
            end
        end
    end

    function plot_changes_by_sector(daily_prof_type, ghost_type, colorbydate, d_km)
        % This function will plot changes in AMF and VCD vs. SCD (trop) for
        % 8 wind sectors around Atlanta. The idea is to try to understand
        % if the changes in column density are systematic in terms of a
        % monthly average.
        %
        % All inputs are optional. The first should be 'regular' or
        % 'hybrid', referring to which daily profile to use (defaults to
        % 'regular'). The second is whether to use "new," "old," or "none"
        % ghost column correction. It defaults to "none." The third should
        % be true or false and indicate whether to color the scatter plots
        % by the dates of the observations. It defaults to false. The
        % last represents how far from Atlanta (in kilometers) to include
        % data for. Defaults to 100 km.
        
        homedir = getenv('HOME');
        if ~exist('daily_prof_type','var')
            daily_prof_type = 'regular';
        elseif ~ischar(daily_prof_type)
            E.badinput('The first input must always be one of the strings "regular" or "hybrid"');
        end
        
        if ~exist('ghost_type','var')
            ghost = ' - No ghost';
        else
            if ~strcmpi(ghost_type,'old') && strcmpi(daily_prof_type, 'regular')
                E.badinput('New ghost products not yet available for regular daily profiles, only hybrid')
            end
            switch lower(ghost_type)
                case 'new'
                    ghost = ' - New ghost';
                case 'none'
                    ghost = ' - No ghost';
                case 'old'
                    ghost = '';
                otherwise
                    E.badinput('ghost_type must be one of ''new'', ''none'', or ''old''')
            end
        end
        
        if ~exist('colorbydate','var')
            colorbydate = false;
        else
            if ~isscalar(colorbydate) || ~islogical(colorbydate)
                E.badinput('colorbydate (if given) must be a scalar logical')
            end
        end
        
        if ~exist('d_km','var')
            d_km = 100;
        elseif ~isnumeric(d_km) || ~isscalar(d_km)
            E.badinput('d_km must be a scalar number, if given');
        end
        
        
        
        monthly_path = fullfile(homedir,'Documents','MATLAB','BEHR','Workspaces','Wind speed',sprintf('SE US BEHR Monthly%s',ghost));
        if strcmpi(daily_prof_type, 'regular')
            daily_path = fullfile(homedir,'Documents','MATLAB','BEHR','Workspaces','Wind speed',sprintf('SE US BEHR Hourly%s',ghost));
        elseif strcmpi(daily_prof_type, 'hybrid')
            daily_path = fullfile(homedir,'Documents','MATLAB','BEHR','Workspaces','Wind speed',sprintf('SE US BEHR Hybrid%s',ghost));
        else
            E.badinput('Daily profile type "%s" not recognized. Options are "regular" or "hybrid"', daily_prof_type);
        end
        
        DF = dir(fullfile(daily_path,'OMI_BEHR_*.mat'));
        MF = dir(fullfile(monthly_path,'OMI_BEHR_*.mat'));
        
        directions = {'ne','n','nw','w','sw','s','se','e'};
        
        amfs_m = make_empty_struct_from_cell(directions);
        amfs_d = make_empty_struct_from_cell(directions);
        vcds_m = make_empty_struct_from_cell(directions);
        vcds_d = make_empty_struct_from_cell(directions);
        scds = make_empty_struct_from_cell(directions);
        dnums = make_empty_struct_from_cell(directions);

        
        for a=1:numel(DF)
            fprintf('Loading file %d\n',a);
            D_OMI = load(fullfile(daily_path,DF(a).name),'OMI');
            [ds,de] = regexp(DF(a).name,'\d\d\d\d\d\d\d\d');
            this_datenum = datenum(DF(a).name(ds:de),'yyyymmdd');
            D_OMI = D_OMI.OMI(2); % both reduce the layers of structures and assume Atlanta is in the second swath
            D_OMI = omi_pixel_reject(D_OMI, 'omi', 0.2, 'XTrackFlags');
            D_OMI.BEHRColumnAmountNO2Trop(D_OMI.Areaweight == 0) = nan;
            D_OMI.BEHRAMFTrop(D_OMI.Areaweight == 0) = nan;
            M_OMI = load(fullfile(monthly_path,MF(a).name),'OMI');
            M_OMI = M_OMI.OMI(2);
            M_OMI = omi_pixel_reject(M_OMI, 'omi', 0.2, 'XTrackFlags');
            M_OMI.BEHRColumnAmountNO2Trop(M_OMI.Areaweight == 0) = nan;
            M_OMI.BEHRAMFTrop(M_OMI.Areaweight == 0) = nan;
            if a == 1
                % First time through, set up the logical matrices for which
                % grid cells in the OMI structure to put into each sector.
                y = D_OMI.Latitude - 33.755;
                x = D_OMI.Longitude - -84.39;
                theta = atan2d(y,x);
                
                % Calculate the distance between each point and Atlanta.
                % We'll cut down on the number of calculations by ignoring
                % a box greater than 1.2x of d_km, assuming 100 km in a
                % degree.
                too_far = abs(x) > d_km/100*1.2 | abs(y) > d_km/100*1.2;
                x(too_far) = nan;
                y(too_far) = nan;
                
                r = inf(size(x));
                for b=1:numel(r)
                    if ~isnan(x(b))
                        r(b) = m_lldist([0, x(b)], [0, y(b)]);
                    end
                end
                
                rr = r <= d_km;
                % Now get the 8 sector arrays.
                xx.ne = rr & theta >= 22.5 & theta < 67.5;
                xx.n = rr & theta >= 67.5 & theta < 112.5;
                xx.nw = rr & theta >= 112.5 & theta < 157.5;
                xx.w = rr & (theta >= 157.5 | theta < -157.5); % slightly different handling for the 180/-180 change
                xx.sw = rr & theta >= -157.5 & theta < -112.5;
                xx.s = rr & theta >= -112.5 & theta < -67.5;
                xx.se = rr & theta >= -67.5 & theta < -22.5;
                xx.e = rr & theta >= -22.5 & theta < 22.5;
            end
            
            for b=1:numel(directions)
                sz = size(M_OMI.BEHRAMFTrop(xx.(directions{b})));
                amfs_m.(directions{b}) = cat(1,amfs_m.(directions{b}),M_OMI.BEHRAMFTrop(xx.(directions{b})));
                amfs_d.(directions{b}) = cat(1,amfs_d.(directions{b}),D_OMI.BEHRAMFTrop(xx.(directions{b})));
                vcds_m.(directions{b}) = cat(1,vcds_m.(directions{b}),M_OMI.BEHRColumnAmountNO2Trop(xx.(directions{b})));
                vcds_d.(directions{b}) = cat(1,vcds_d.(directions{b}),D_OMI.BEHRColumnAmountNO2Trop(xx.(directions{b})));
                scds.(directions{b}) = cat(1,scds.(directions{b}),M_OMI.ColumnAmountNO2Trop(xx.(directions{b})) .* M_OMI.AMFTrop(xx.(directions{b})));
                scds.(directions{b})(scds.(directions{b})<0) = nan;
                dnums.(directions{b}) = cat(1, dnums.(directions{b}), repmat(this_datenum, sz));
            end
            
        end
        
        % Make 16 plots (woohoo) - change in AMF and VCD vs. SCD for each
        % sector
        
        for b=1:numel(directions)
            figure; 
            if colorbydate
                scatter(scds.(directions{b}), amfs_d.(directions{b}) - amfs_m.(directions{b}), 16, dnums.(directions{b}));
            else
                scatter(scds.(directions{b}), amfs_d.(directions{b}) - amfs_m.(directions{b}));
            end
            xlabel('Tropospheric slant column density (molec. cm^{-2})');
            ylabel('\Delta AMF (daily - monthly)');
            set(gca,'fontsize',16);
            title(sprintf('%s sector \\Delta AMF using %s daily profiles',upper(directions{b}),daily_prof_type));
            
            figure; 
            if colorbydate
                scatter(scds.(directions{b}), vcds_d.(directions{b}) - vcds_m.(directions{b}), 16, dnums.(directions{b}));
            else
                scatter(scds.(directions{b}), vcds_d.(directions{b}) - vcds_m.(directions{b}));
            end
            xlabel('Tropospheric slant column density (molec. cm^{-2})');
            ylabel('\Delta VCD (daily - monthly)');
            set(gca,'fontsize',16);
            title(sprintf('%s sector \\Delta VCD using %s daily profiles',upper(directions{b}),daily_prof_type));
        end
        
    end

    function plot_diff_resolutions
        % First we must find out what difference to plot. Need to know both
        % what profile type and what resolution, as well as whether or not
        % to use interpolation and what date to do this for.
        try
            allowed_plot_types = {'a','p','v'};
            plot_type = ask_multichoice('Do you want to plot an absolute difference, percent difference, or just the values for one case?', allowed_plot_types);
            
            allowed_prof_types = {'monthly','hourly','hybrid'};
            prof_type_base = ask_multichoice('What profile type to use for the base case?', allowed_prof_types);
            prof_type_base(1) = upper(prof_type_base(1));
            
            allowed_quantities = {'amf','vcd'};
            quantity = ask_multichoice('Which quantity do you want to compare?', allowed_quantities);

            allowed_resolutions = {'12','24','96','288'};
            res_base = ask_multichoice('What resolution to use for the base case? ', allowed_resolutions, 'default', '12');
            
            if ~strcmpi(res_base, '12')
                allowed_interps = {'y','n'};
                interp_base = ask_multichoice('Do you want to use interpolation for the base case?', allowed_interps);
            else
                % The 12-km runs do not have interpolation
                interp_base = 'n';
            end
            
            if ismember(plot_type,{'a','p'})
                prof_type_new = ask_multichoice('What profile type to use for the new case?', allowed_prof_types, 'default', lower(prof_type_base));
                res_new = ask_multichoice('What resolution to use for the new case?', allowed_resolutions);
                if ~strcmpi(res_new, '12')
                    allowed_interps = {'y','n'};
                    interp_new = ask_multichoice('Do you want to use interpolation for the new case?', allowed_interps, 'default', interp_base);
                else
                    % The 12-km runs do not have interpolation
                    interp_new = 'n';
                end
            end
        catch err
            % This lets us exit gracefully if requested without having to
            % check every time that the function returned a 0.
            if strcmpi(err.identifier, 'ask_multichoice:user_cancel')
                return
            else
                rethrow(err)
            end
        end
        
        date_to_comp = input('Enter the date to compare using a format datestr can recognize: ','s');
        try
            date_to_comp = datestr(date_to_comp, 'yyyymmdd');
        catch err
            if strcmpi(err.identifier,'MATLAB:datestr:ConvertToDateNumber')
                E.badinput('Format for date not recognized. Try yyyy-mm-dd.');
            else
                rethrow(err);
            end
        end
        
        
        % Now we start actually comparing. First we need to get the files
        % to load.
        if strcmpi(res_base, '12')
            res_str = '';
        else
            res_str = sprintf(' - %s km resolution', res_base);
        end
        if strcmpi(interp_base,'y')
            base_interp_str = ' - Interpolated';
        else
            base_interp_str = '';
        end
        homedir = getenv('HOME');
        if strcmp(quantity,'amf')
            base_filedir = sprintf('Atlanta BEHR %s - No clouds%s%s', prof_type_base, res_str, base_interp_str);
        else
            base_filedir = sprintf('SE US BEHR %s - No ghost%s', prof_type_base, res_str);
        end
        base_filename = fullfile(homedir,'Documents','MATLAB','BEHR','Workspaces','Wind speed',base_filedir,sprintf('OMI_BEHR_%s.mat', date_to_comp));
        if ~exist(base_filename,'file') && strcmp(quantity,'vcd')
            fprintf('Only a limited subset of options for comparing VCDs are available.\n')
            return
        end
        Base = load(base_filename);
        
        if ismember(plot_type,{'a','p'})
            if strcmpi(res_new, '12')
                res_str = '';
            else
                res_str = sprintf(' - %s km resolution', res_new);
            end
            if strcmpi(interp_new,'y')
                new_interp_str = ' - Interpolated';
            else
                new_interp_str = '';
            end
            homedir = getenv('HOME');
            if strcmp(quantity,'amf')
                new_filedir = sprintf('Atlanta BEHR %s - No clouds%s%s', prof_type_new, res_str, new_interp_str);
            else
                new_filedir = sprintf('SE US BEHR %s - No ghost%s', prof_type_new, res_str);
            end
            new_filename = fullfile(homedir,'Documents','MATLAB','BEHR','Workspaces','Wind speed',new_filedir,sprintf('OMI_BEHR_%s.mat', date_to_comp));
            if ~exist(base_filename,'file') && strcmp(quantity,'vcd')
                fprintf('Only a limited subset of options for comparing VCDs are available.\n')
                return
            end
            New = load(new_filename);
        end
        % Now load and plot. Will need to cut down the VCD array if it is
        % requested
                
        
        if strcmp(quantity,'amf')
            lonxy = squeeze(Base.Data.Loncorn(1,:,:));
            latxy = squeeze(Base.Data.Latcorn(1,:,:));
            if strcmp(plot_type,'p')
                perdiff = (New.Data.BEHRAMFTrop ./ Base.Data.BEHRAMFTrop - 1)*100;
            elseif strcmp(plot_type,'a')
                perdiff = New.Data.BEHRAMFTrop - Base.Data.BEHRAMFTrop;
            else
                perdiff = Base.Data.BEHRAMFTrop;
            end
        else
            xx = 761:860;
            yy = 141:210;
            lonxy = Base.OMI(2).Longitude(yy,xx);
            latxy = Base.OMI(2).Latitude(yy,xx);
            if strcmp(plot_type, 'p')
                perdiff = (New.OMI(2).BEHRColumnAmountNO2Trop(yy,xx) ./ Base.OMI(2).BEHRColumnAmountNO2Trop(yy,xx) - 1)*100;
            elseif strcmp(plot_type, 'a')
                perdiff = New.OMI(2).BEHRColumnAmountNO2Trop(yy,xx) - Base.OMI(2).BEHRColumnAmountNO2Trop(yy,xx);
            else
                perdiff = Base.OMI(2).BEHRColumnAmountNO2Trop(yy,xx);
            end
        end
        
        figure; 
        pcolor(lonxy, latxy, perdiff);
        cb = colorbar;
        xlim([-87.1 -82]);
        ylim([32, 35.6]);
        
        l = line(-84.39, 33.755, 'linestyle','none', 'marker','p','markersize',18,'color','k','linewidth',2);
        legend(l,{'Atlanta'});
        
        set(gca,'fontsize',16);
        
        if ismember(plot_type,{'a','p'})
            title(sprintf('%s: %s km %s vs %s km %s', date_to_comp, res_new, new_interp_str, res_base, base_interp_str));
            % ensure colorbar limits are a multiple of 10 and equal positive
            % and negative if doing a difference
            maxval = max(abs(get(gca,'clim')));
            maxval = ceil(maxval/10)*10;
            caxis([-maxval maxval]);
            
            C = load('blue_red_cmap.mat');
            colormap(C.blue_red_cmap);
            
            if strcmp(quantity,'amf')
                if strcmp(plot_type,'p')
                    cb.Label.String = '%\Delta AMF';
                else
                    cb.Label.String = 'Delta AMF';
                end
            else
                shading flat
                if strcmp(plot_type,'p')
                    cb.Label.String = '%\Delta VCD';
                else
                    cb.Label.String = '\Delta VCD';
                end
            end
        else
            title(sprintf('%s: %s km %s', date_to_comp, res_base, base_interp_str));
            cb.Label.String = 'AMF';
        end
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% OTHER FUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

