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
    case 'plotavgapriori'
        plot_avg_apriori(varargin{:});
    case 'getnearprof'
        varargout{1} = get_nearest_profiles(varargin{:});
    case 'calcwind'
        varargout{1} = calc_wind_mag(varargin{1}, varargin{2});
        varargout{2} = calc_wind_dir(varargin{1}, varargin{2});
    case 'calcavgwind'
        [varargout{1}, varargout{2}] = calc_avg_wind(varargin{:});
    case 'perdiffvstheta'
        plot_delamf_vs_delangle(varargin{:});
    case 'perrec'
        plot_perrec_vs_distance(varargin{:});
    case 'sectors'
        plot_changes_by_sector(varargin{:});
    case 'dif'
        plot_diff(varargin{:});
    case 'dif-ts'
        plot_pseudo_diff_timeser();
    case 'db-avg'
        plot_behr_avg();
    case 'pr-prof'
        [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}] = plot_pseudo_apriori(varargin{:});
    case 'ens-prof'
        varargout{1} = plot_ens_apriori();
    case 'full-prof-diff'
        full_apriori_avg_diff(varargin{:});
    case 'cld'
        plot_cloudfrac(varargin{1});
    case 'res'
        plot_diff_resolutions();
    case 'cat-cld-bins'
        varargout{1} = cat_cloud_frac_bins(varargin{:});
    case 'del-shape-vcd'
        del_shape_vs_vcd();
    case 'mag-delvcd'
        mag_delvcd_stats();
    case 'pix-pos-v-cld'
        pix_pos_vs_cld();
    case 'wind-v-cld'
        wind_cond_vs_cld();
    case 'apriori-surface'
        plot_behr_apriori_surface();
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

    function plot_avg_apriori(indicies, startdate, enddate)
        if ~isnumeric(indicies) || ~ismatrix(indicies) || size(indicies,2) ~= 2
            E.badinput('Indicies must be an n-by-2 matrix')
        end
        sdate=datenum(startdate);
        edate=datenum(enddate);
        
        hourly_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hourly - No clouds - No ghost';
        monthly_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Monthly - No clouds - No ghost';
        
        F_hr = dir(fullfile(hourly_path, 'OMI*.mat'));
        F_mn = dir(fullfile(monthly_path, 'OMI*.mat'));
        
        if numel(F_hr) ~= numel(F_mn)
            E.callError('different_num_files','Different number of hybrid and monthly files');
        end
        
        hr_profs = nan(30,numel(F_hr),size(indicies,1));
        hr_pres = nan(30,size(indicies,1));
        mn_profs = nan(30,numel(F_mn),size(indicies,1));
        mn_pres = nan(30,size(indicies,1));
        
        lon = nan(1,size(indicies,1));
        lat = nan(1,size(indicies,1));
        
        first_time = true;
        
        for a=1:numel(F_hr)
            [s,e] = regexp(F_hr(a).name,'\d\d\d\d\d\d\d\d');
            fdate = datenum(F_hr(a).name(s:e),'yyyymmdd');
            if fdate >= sdate && fdate <= edate
                hr = load(fullfile(hourly_path, F_hr(a).name),'Data');
                mn = load(fullfile(monthly_path, F_mn(a).name),'Data');
                
                for b=1:size(indicies,1)
                    x=indicies(b,1);
                    y=indicies(b,2);
                    hr_profs(:,a,b) = hr.Data.BEHRNO2apriori(:,x,y);
                    mn_profs(:,a,b) = mn.Data.BEHRNO2apriori(:,x,y);
                    if first_time
                        % In the pseudo retrieval, these should stay the
                        % same across all days
                        hr_pres(:,b) = hr.Data.BEHRPressureLevels(:,x,y);
                        mn_pres(:,b) = mn.Data.BEHRPressureLevels(:,x,y);
                        lon(b) = hr.Data.Longitude(x,y);
                        lat(b) = hr.Data.Latitude(x,y);
                        first_time = false;
                    end
                end
            end
        end
        
        hr_avg_prof = squeeze(nanmean(hr_profs,2));
        mn_avg_prof = squeeze(nanmean(mn_profs,2));
        
        for b=1:size(indicies,1)
            figure;
            line(hr_avg_prof(:,b), hr_pres(:,b), 'color','b','linewidth',2);
            line(mn_avg_prof(:,b), mn_pres(:,b), 'color','r','linewidth',2);
            legend('Hourly','Monthly');
            set(gca,'ydir','reverse');
            set(gca,'fontsize',16);
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

    function plot_diff(options)
        % This function can use no input, instead it will ask a series of
        % questions to decide what plots to make.  This is used to make
        % difference pcolor plots of WRF-Chem VCDs or BEHR VCDs or AMFs.
        %
        % If given input in the form of the structure "options" it will
        % only ask the necessary questions to fill in the missing options.
        % Options can have the following fields:
        %   source (wrf, behr, pseudo-behr)
        %   apriori_base (hourly, hybrid, monthly)
        %   apriori_new (hourly, hybrid, monthly)
        %   res (f, c)
        %   timemode (avg, instant)
        %   quantity (vcd, amf, wind)
        %   coarsen (integer > 0)
        %   diff_type (a, p, d, m)
        %   city (Atlanta, Birmingham, Montgomery)
        %   date_in (valid datestring)
        
        %%%%%%%%%%%%%%%%%
        %%%%% INPUT %%%%%
        %%%%%%%%%%%%%%%%%
        
        if exist('options','var') && ~isstruct(options)
            E.badinput('OPTIONS must be a structure (if given)');
        elseif ~exist('options','var')
            options = struct;
        end
        
        % First question: WRF or BEHR. pseudo-BEHR is the one where I used
        % the same set of pixels each day.
        allowed_sources = {'wrf','behr','pseudo-behr'};
        if ~isfield(options,'source')
            source = ask_multichoice('Which source will you be using?', allowed_sources);
        else
            if ~ismember(options.source, allowed_sources)
                E.badinput('options.source must be one of %s',strjoin(allowed_sources,', '))
            end
            source = options.source;
        end
        
        % Follow up if using BEHR: should we compare hybrid or hour-wise
        % data to the monthly average? Allow the user to manually choose
        % two directories to compare if the a priori they want is not
        % listed. This will remove the need for some later questions as
        % well.
        manual_dir = false;
        if ~strcmpi(source,'wrf')
            manual_dir = false;
            if strcmpi(source,'pseudo-behr')
                allowed_apriori = {'hourly','hybrid','hybrid-avg','monthly'};
            else
                allowed_apriori = {'hourly','hybrid','monthly','monthly-converg','monthly-sqrt-converg'};
            end
            if any(~isfield(options,{'apriori_base','apriori_new'}))
                apriori_base = ask_multichoice('Which a priori will be the base case?', [allowed_apriori, {'Choose directory manually'}], 'default', 'monthly','list',true);
                if ~strcmpi(apriori_base,'Choose directory manually') || ~isDisplay
                    if strcmpi(apriori_base,'Choose directory manually') && ~isDisplay
                        fprintf('Not using a display, cannot open UIGETDIR dialogue. Please choose an a priori.\n')
                        apriori_base = ask_multichoice('Which a priori will be the base case?', allowed_apriori, 'default', 'monthly','list',true);
                    end
                    apriori_new = ask_multichoice('Which a priori will be the new case?', allowed_apriori, 'list', true);
                else
                    fprintf('Choose the base apriori first, then the new apriori\n');
                    input('Press ENTER to continue','s');
                    base_dir = uigetdir('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed', 'Choose the base a priori folder');
                    new_dir = uigetdir('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed', 'Choose the new a priori folder');
                    if isnumeric(base_dir) && base_dir == 0 || isnumeric(new_dir) && new_dir == 0
                        E.userCancel;
                    end
                    manual_dir = true;
                    
                    % Check that the directory name includes "Atlanta" if
                    % doing pseudo retrieval or "SE US"/"W US" for full
                    % BEHR. This is just how I named the folders.
                    if strcmpi(source,'pseudo-behr')
                        teststr = 'Atlanta';
                    else
                        teststr = '(SE US|W US)';
                    end
                    btest = isempty(regexp(base_dir, teststr, 'once'));
                    ntest = isempty(regexp(new_dir, teststr, 'once'));
                    if btest || ntest
                        fprintf('One of the directories does not appear to be a %s directory (did not pass the regexp test "%s"\n', source, teststr);
                        if strcmpi(ask_multichoice('Continue?', {'y','n'}, 'default', 'n'),'n')
                            E.userCancel;
                        end
                    end
                end
            else
                if any(~ismember({options.apriori_base, options.apriori_new}, allowed_apriori))
                    E.badinput('options.apriori_base and options.apriori_new must be one of %s',strjoin(allowed_apriori,', '))
                end
                apriori_base = options.apriori_base;
                apriori_new = options.apriori_new;
            end
        else 
            apriori_new = 'hourly';
            apriori_base = 'monthly';
        end
        
        % Do we want to use the fine or coarse WRF simulation?
        allowed_res = {'f','c'};
        if manual_dir
            %do nothing
        elseif ~isfield(options,'res_base') || ~isfield(options,'res_new')
            res_base = ask_multichoice('Do you want the fine (12 km) or coarse (108 km) WRF for the base a priori', allowed_res);
            res_new = ask_multichoice('And for the new a priori?', allowed_res, 'default', res_base);
        else
            if ~ismember(options.res, allowed_res)
                E.badinput('options.res must be one of %s',strjoin(allowed_res,', '))
            end
            res_base = options.res_base;
            res_new = options.res_new;
        end
        
        % Use the hour-average or instantaneous profiles?
        if manual_dir
            %do nothing
        elseif strcmpi(source,'pseudo-behr') && any(ismember({apriori_base, apriori_new},{'hourly','hybrid'}) & strcmpi({res_base, res_new},'f'))
            allowed_timemode = {'avg','instant'};
            if ~isfield(options,'timemode')
                timemode = ask_multichoice('Use profiles averaged over an hour or instantaneous at the top of the hour?', allowed_timemode);
            else
                if ~ismember(options.timemode, allowed_timemode)
                    E.badinput('options.timemode must be one of %s',strjoin(allowed_timemode,', '))
                end
                timemode = options.timemode;
            end
        else
            timemode = 'avg';
        end
        
        % Which quantity. WRF has VCD and winds. Full BEHR
        % has VCDs or AMFs. Pseudo-behr only has AMFs.
        coarsen = 1;
        switch source
            case 'wrf'
                allowed_quantities = {'vcd','wind'};
                if ~isfield(options,'quantity')
                    quantity = ask_multichoice('Which source will you be using?', allowed_quantities);
                else
                    if ~ismember(options.quantity, allowed_quantities)
                        E.badinput('options.quantity must be one of %s',strjoin(allowed_quantities,', '))
                    end
                    quantity = options.quantity;
                end
            case 'behr'
                allowed_quantities = {'amf','vcd'};
                if ~isfield(options,'quantity')
                    quantity = ask_multichoice('Which source will you be using?', allowed_quantities);
                else
                    if ~ismember(options.quantity, allowed_quantities)
                        E.badinput('options.quantity must be one of %s',strjoin(allowed_quantities,', '))
                    end
                    quantity = options.quantity;
                end
            otherwise
                quantity = 'amf';
        end
        
        % Give the option to make winds coarser so the arrows aren't so
        % small
        if strcmp(quantity,'wind')
            coarse_test = @(x) mod(x,1) == 0 && x > 0;
            if ~isfield(options,'coarsen')
                coarsen = ask_number('Do you want to average winds together? Enter an integer of 1 if not, or greater if yes','default',1,'testfxn',coarse_test,'testmsg','Must be an integer >= 1');
            else
                if ~coarse_test(options.coarsen)
                    E.badinput('options.coarsen must be an integer of 1 or greater')
                end
                coarsen = options.coarsen;
            end
            plot_wind_arrow = false;
        else
            plot_wind_arrow = strcmpi(ask_multichoice('Plot an arrow indicating wind direction over the city?',{'y','n'}),'y');
        end
        
        % Fourth, is this a percent or absolute difference?
        if strcmpi(quantity,'wind')
            allowed_diff = {'d','m'};
            ask_str = 'Do you want the daily or monthly winds?';
            
        else
            allowed_diff = {'a','p','d','m'};
            ask_str = 'Do you want absolute or percent differences, or just the daily or monthly values?';
        end
        if ~isfield(options,'diff_type')
            diff_type = ask_multichoice(ask_str, allowed_diff);
        else
            if ~ismember(options.diff_type, allowed_diff)
                E.badinput('options.diff_type must be one of %s',strjoin(allowed_diff,', '))
            end
            diff_type = options.diff_type;
        end
        
        % Which city to focus on?
        allowed_cities = {'Atlanta','Birmingham','Montgomery','SF'};
        if ~isfield(options,'city')
            city = ask_multichoice('Which city to focus on?',allowed_cities,'list',true,'default','Atlanta');
        else
            if ~ismember(options.city, allowed_cities)
                E.badinput('options.city must be one of %s',strjoin(allowed_cities,', '))
            end
            city = options.city;
        end
        
        % Lastly, we need a date
        if ~isfield(options,'date_in')
            date_in = ask_date('Enter the date to compare, using a format datenum can parse');
        else
            try
                datenum(options.date_in)
            catch err
                if strcmpi(err.identifier,'MATLAB:datenum:ConvertDateString')
                    E.badinput('options.date_in must be a valid date string')
                else
                    rethrow(err);
                end
            end
            date_in = options.date_in;
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% PARSING INPUT %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        switch city
            case 'Atlanta'
                %xl = [-86 -82];
                %yl = [32.5 35.5];
                xl = [-87.1 -81.9];
                yl = [31.9 35.5];
                city_lon = -84.39;
                city_lat = 33.775;
                city_swaths = 2;
                windfile = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta-Wind-Conditions-1900UTC-5layers-earthrel.mat';
            case 'Birmingham'
                xl = [-88 -85.5];
                yl = [33 34.5];
                city_lon = -86.80;
                city_lat = 33.52;
                city_swaths = 2;
                windfile = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Birmingham-Wind-Conditions-1900UTC-5layers-earthrel.mat';
            case 'Montgomery'
                xl = [-87 -85];
                yl = [31.5 33];
                city_lon = -86.3;
                city_lat = 32.37;
                city_swaths = 2;
                windfile = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Montgomery-Wind-Conditions-1900UTC-5layers-earthrel.mat';
            case 'SF'
                xl = [-125 -119];
                yl = [32 42];
                city_lon = -122.42;
                city_lat = 37.77;
                city_swaths = [3,4];
                windfile = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SF-Wind-Conditions-1900UTC-5layers-earthrel-earthrel.mat';
            otherwise
                E.notimplemented('city = %s',city);
        end
        city_name = city;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% LOAD DATA, CALCULATE QUANTITIES, AND PLOT %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~manual_dir
            base_file = return_data_file_info(source, city, timemode, res_base, apriori_base, date_in);
            new_file = return_data_file_info(source, city, timemode, res_new, apriori_new, date_in);
        else
            if strcmpi(source,'wrf')
                E.notimplemented('Manual WRF directory - monthly names are different')
            end
            fname = sprintf('OMI_BEHR_%s.mat', datestr(date_in, 'yyyymmdd'));
            base_file = fullfile(base_dir, fname);
            new_file = fullfile(new_dir, fname);
        end
        
        if strcmp(source,'wrf')
            % If we're not doing a difference, point both base and new to
            % the same file. This skirts around an issue with loading the
            % lat/lon
            if strcmpi(diff_type,'d')
                base_file = new_file;
            elseif strcmpi(diff_type,'m')
                new_file = base_file;
            end
            
            XLONG = ncread(new_file, 'XLONG');
            XLAT = ncread(new_file, 'XLAT');
            
            xx = any(XLONG >= min(xl) & XLONG <= max(xl) & XLAT >= min(yl) & XLAT <= max(yl),2);
            yy = any(XLONG >= min(xl) & XLONG <= max(xl) & XLAT >= min(yl) & XLAT <= max(yl),1);
            xx = squeeze(xx(:,:,1));
            yy = squeeze(yy(:,:,1));

            % add a few more to fill out the edges
            xx = find(xx);
            xx = (min(xx)-2):(max(xx)+2);
            yy = find(yy);
            yy = (min(yy)-2):(max(yy)+2);

            
            xlon = XLONG(xx,yy,1);
            xlat = XLAT(xx,yy,1);
            
            utchr = ncread(new_file, 'utchr');
            s = utchr == 19;
            
            if strcmpi(quantity, 'vcd')
                [xlon_pcol, xlat_pcol] = pcolor_pixel_corners(xlon,xlat);
                
                daily_no2 = ncread(new_file, 'no2_ndens'); %[NO2 in number density]
                daily_no2 = nanmean(daily_no2(xx,yy,:,s),4); % cut down to the hour we want. If there are multiple outputs, average them.
                daily_zlev = ncread(new_file, 'zlev'); % Thickness of each layer in meters
                daily_zlev = nanmean(daily_zlev(xx,yy,:,s),4);
                daily_tplev = find_wrf_tropopause(ncinfo(new_file));
                daily_tplev = floor(nanmean(daily_tplev(xx,yy,s),3));
                
                for a=1:size(daily_no2,1)
                    for b=1:size(daily_no2,2)
                        tp = daily_tplev(a,b);
                        if tp > 0 % tp is given -1 if the tropopause algorithm cannot find a tropopause
                            daily_no2(a,b,tp:end) = nan;
                        end
                    end
                end
                
                daily_no2_columns = nansum2(daily_no2 .* (daily_zlev*100), 3);
                
                monthly_no2 = ncread(base_file, 'no2_ndens'); %[NO2 in number density]
                monthly_no2 = monthly_no2(xx,yy,:); % cut down to the hour we want
                monthly_zlev = ncread(base_file, 'zlev'); % Thickness of each layer in meters
                monthly_zlev = monthly_zlev(xx,yy,:);
                monthly_tplev = find_wrf_tropopause(ncinfo(base_file));
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
                        cmap = blue_red_cmap;
                    case 'p'
                        del = (daily_no2_columns ./ monthly_no2_columns - 1)*100;
                        cmap = blue_red_cmap;
                    case 'd'
                        del = daily_no2_columns;
                        cmap = 'jet';
                    case 'm'
                        del = monthly_no2_columns;
                        cmap = 'jet';
                end
            
                f = figure; pcolor(xlon_pcol, xlat_pcol, del);
                cb=colorbar;
                set(gca,'fontsize',20);
                xlim(xl);
                ylim(yl);
                colormap(cmap);
            elseif strcmpi(quantity, 'wind')
                if strcmpi(diff_type,'d')
                    U = ncread(new_file, 'U');
                    U = squeeze(nanmean(U(:,:,1,s),4));
                    V = ncread(new_file, 'V');
                    V = squeeze(nanmean(V(:,:,1,s),4));
                    COSALPHA = ncread(new_file, 'COSALPHA');
                    COSALPHA = COSALPHA(:,:,1); % should not change with time
                    SINALPHA = ncread(new_file, 'SINALPHA');
                    SINALPHA = SINALPHA(:,:,1);
                elseif strcmpi(diff_type,'m')
                    U = ncread(base_file, 'U');
                    U = squeeze(U(:,:,1));
                    V = ncread(base_file, 'V');
                    V = squeeze(V(:,:,1));
                    COSALPHA = ncread(base_file, 'COSALPHA');
                    SINALPHA = ncread(base_file, 'SINALPHA');
                end
                [U_e, V_e] = wrf_winds_transform(U,V,COSALPHA,SINALPHA);
                % want to cut these down only after they are unstaggered.
                U_e = U_e(xx,yy);
                V_e = V_e(xx,yy);
                
                if coarsen > 1
                    xlon_new = nan(floor(size(xlon)/coarsen));
                    xlat_new = nan(floor(size(xlat)/coarsen));
                    U_new = nan(floor(size(U_e)/coarsen));
                    V_new = nan(floor(size(V_e)/coarsen));
                    for a=1:size(U_new,1)
                        i=(a-1)*coarsen + 1;
                        for b=1:size(V_new,2)
                            j=(b-1)*coarsen + 1;
                            xlon_tmp = xlon(i:(i+coarsen-1),j:(j+coarsen-1));
                            xlon_new(a,b) = nanmean(xlon_tmp(:));
                            xlat_tmp = xlat(i:(i+coarsen-1),j:(j+coarsen-1));
                            xlat_new(a,b) = nanmean(xlat_tmp(:));
                            u_tmp = U_e(i:(i+coarsen-1),j:(j+coarsen-1));
                            U_new(a,b) = nanmean(u_tmp(:));
                            v_tmp = V_e(i:(i+coarsen-1),j:(j+coarsen-1));
                            V_new(a,b) = nanmean(v_tmp(:));
                        end
                    end
                    xlon = xlon_new;
                    xlat = xlat_new;
                    U_e = U_new;
                    V_e = V_new;
                end
                
                f = figure; quiver(xlon,xlat,U_e,V_e,'linewidth',2,'color','k');
                set(gca,'fontsize',20);
                xlim(xl);
                ylim(yl);
            end
        elseif ~isempty(regexp(source, 'behr', 'once'))
            if strcmp(source,'behr')
                D = load(new_file,'Data');
                M = load(base_file,'Data');
                lon = cell(1,numel(city_swaths));
                lat = cell(1,numel(city_swaths));
                daily_value = cell(1,numel(city_swaths));
                monthly_value = cell(1,numel(city_swaths));
                for i=1:numel(city_swaths)
                    s = city_swaths(i);
                    minloncorn = squeeze(min(D.Data(s).Loncorn,[],1));
                    maxloncorn = squeeze(max(D.Data(s).Loncorn,[],1));
                    minlatcorn = squeeze(min(D.Data(s).Latcorn,[],1));
                    maxlatcorn = squeeze(max(D.Data(s).Latcorn,[],1));
                    % Get any pixel even slightly in the domain.
                    xx = any(maxloncorn >= min(xl) & minloncorn <= max(xl) & maxlatcorn >= min(yl) & minlatcorn <= max(yl),2);
                    yy = any(maxloncorn >= min(xl) & minloncorn <= max(xl) & maxlatcorn >= min(yl) & minlatcorn <= max(yl),1);
                    
                    D_badpix = find_bad_pixels(D.Data(s));
                    M_badpix = find_bad_pixels(M.Data(s));
                    badpix = D_badpix | M_badpix;
                    
                    lon{i} = squeeze(D.Data(s).Loncorn(1,xx,yy));
                    lat{i} = squeeze(D.Data(s).Latcorn(1,xx,yy));
                    
                    if strcmp(quantity, 'amf')
                        D.Data(s).BEHRAMFTrop(badpix) = nan;
                        daily_value{i} = D.Data(s).BEHRAMFTrop(xx,yy);
                        M.Data(s).BEHRAMFTrop(badpix) = nan;
                        monthly_value{i} = M.Data(s).BEHRAMFTrop(xx,yy);
                    else
                        D.Data(s).BEHRColumnAmountNO2Trop(badpix) = nan;
                        daily_value{i} = D.Data(s).BEHRColumnAmountNO2Trop(xx,yy);
                        M.Data(s).BEHRColumnAmountNO2Trop(badpix) = nan;
                        monthly_value{i} = M.Data(s).BEHRColumnAmountNO2Trop(xx,yy);
                    end
                end
            elseif strcmp(source,'pseudo-behr')
                D = load(new_file, 'Data');
                M = load(base_file, 'Data');
                
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
            
            f = gobjects(1, numel(city_swaths));
            cb = gobjects(1, numel(city_swaths));
            if ~iscell(lon); lon = {lon}; end
            if ~iscell(lat); lat = {lat}; end
            if ~iscell(daily_value); daily_value = {daily_value}; end
            if ~iscell(monthly_value); monthly_value = {monthly_value}; end
            
            for i=1:numel(city_swaths)
                switch diff_type
                    case 'a'
                        del = daily_value{i} - monthly_value{i};
                        cmap = four_color_cmap;
                    case 'p'
                        del = (daily_value{i} ./ monthly_value{i} - 1) * 100;
                        cmap = four_color_cmap;
                    case 'd'
                        del = daily_value{i};
                        cmap = 'jet';
                    case 'm'
                        del = monthly_value{i};
                        cmap = 'jet';
                end

                f(i) = figure; pcolor(lon{i}, lat{i}, del);
                colormap(cmap);
                cb(i)=colorbar;
                cb(i).FontSize = 20;
                set(gca,'fontsize',20);
                xlim(xl);
                ylim(yl);
            end
            
        end
        
        for i=1:numel(f)
            figure(f(i));
            if ~strcmpi(quantity,'wind')
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
                        label_pt1 = '';
                    case 'm'
                        label_pt1 = '';
                end
                
                
                cb(i).Label.String = strjoin({label_pt1, label_pt2, label_unit}, ' ');
                
                cm = max(abs(del(:)));
                if ~isnan(cm)
                    if ismember(diff_type,{'p','a'})
                        caxis([-cm cm]);
                    else
                        caxis([0 cm]);
                    end
                end
            end
            
            % Plot the wind arrow over the city if requested
            if plot_wind_arrow
                W = load(windfile);
                dd = W.dnums == datenum(date_in);
                theta = W.theta(dd);
                alon = [city_lon-0.5*cosd(theta), city_lon+0.5*cosd(theta)];
                alat = [city_lat-0.5*sind(theta), city_lat+0.5*sind(theta)];
                % Move the arrow off the city a little bit.
                alon = alon + 0.25*cosd(theta+90);
                alat = alat + 0.25*sind(theta+90);
                
                draw_arrow(alon, alat, 'linewidth', 3, 'color', 'k','maxheadsize',0.75);
            end
            
            col='k';
            l=line(city_lon, city_lat, 'linestyle','none', 'marker','p','markersize',18,'color',col,'linewidth',2);
            leg=legend(l,city_name);
        end
    end

    function plot_pseudo_diff_timeser()
        allowed_diffs = {'hr-hy','hy-mn','hr-mn','hy-avg','avg-mn','all','Pick directories manually'};
        diff_mode = ask_multichoice('Which difference to consider; hourly vs hybrid or hybrid vs monthly?', allowed_diffs,'list',true);
        if ~strcmpi(diff_mode,'all')
            allowed_modes = {'box','dist','scatter-dist','scatter-dist-wbox','scatter-angle','scatter-angle-wbox','pcolor','pcolor-med','pcolor-apri','pcolor-apri-stdp','combo'};
            plot_mode = ask_multichoice(sprintf('Which type of plot do you want:\n'), allowed_modes);
        else
            plot_mode = 'box-comp';
        end
        start_date = ask_date('Enter the start date');
        end_date = ask_date('Enter the end date');
        city_lon = -84.39;
        city_lat = 33.775;
        
        workdir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/';
        switch lower(diff_mode)
            case 'hr-hy'
                new_dir = fullfile(workdir,'Atlanta BEHR Hourly - No clouds - No ghost - UTC 1800-2200');
                F_new = dir(fullfile(new_dir,'OMI*.mat'));
                old_dir = fullfile(workdir,'Atlanta BEHR Hybrid - No clouds - No ghost - lonweight 13.5 overpass - UTC 1800-2200');
                F_old = dir(fullfile(old_dir,'OMI*.mat'));
            case 'hy-mn'
                new_dir = fullfile(workdir,'Atlanta BEHR Hybrid - No clouds - No ghost - lonweight 13.5 overpass - UTC 1800-2200');
                F_new = dir(fullfile(new_dir,'OMI*.mat'));
                old_dir = fullfile(workdir,'Atlanta BEHR Monthly - No clouds - No ghost - lonweight 13.5 overpass - UTC 1800-2200');
                F_old = dir(fullfile(old_dir,'OMI*.mat'));
            case 'hy-avg'
                new_dir = fullfile(workdir,'Atlanta BEHR Hybrid - No clouds - No ghost - lonweight 13.5 overpass - UTC 1800-2200');
                F_new = dir(fullfile(new_dir,'OMI*.mat'));
                old_dir = fullfile(workdir,'Atlanta BEHR Avg Hybrid - No clouds - No ghost');
                F_old = dir(fullfile(old_dir,'OMI*.mat'));
            case 'hr-mn'
                new_dir = fullfile(workdir,'Atlanta BEHR Hourly - No clouds - No ghost - UTC 1800-2200');
                F_new = dir(fullfile(new_dir,'OMI*.mat'));
                old_dir = fullfile(workdir,'Atlanta BEHR Monthly - No clouds - No ghost - lonweight 13.5 overpass - UTC 1800-2200');
                F_old = dir(fullfile(old_dir,'OMI*.mat'));
            case 'avg-mn'
                new_dir = fullfile(workdir,'Atlanta BEHR Avg Hybrid - No clouds - No ghost');
                F_new = dir(fullfile(new_dir,'OMI*.mat'));
                old_dir = fullfile(workdir,'Atlanta BEHR Monthly - No clouds - No ghost - lonweight 13.5 overpass - UTC 1800-2200');
                F_old = dir(fullfile(old_dir,'OMI*.mat'));
            case 'all'
                hr_dir = fullfile(workdir,'Atlanta BEHR Hourly - No clouds - No ghost - UTC 1800-2200');
                F_hr = dir(fullfile(hr_dir,'OMI*.mat'));
                hy_dir = fullfile(workdir,'Atlanta BEHR Hybrid - No clouds - No ghost - lonweight 13.5 overpass - UTC 1800-2200');
                F_hy = dir(fullfile(hy_dir,'OMI*.mat'));
                mn_dir = fullfile(workdir,'Atlanta BEHR Monthly - No clouds - No ghost - lonweight 13.5 overpass - UTC 1800-2200');
                F_mn = dir(fullfile(mn_dir,'OMI*.mat'));
            case 'pick directories manually'
                if ~isDisplay
                    E.notdisplay('Picking directories manually requires a display.')
                end
                fprintf('Two dialogue boxes will now open. Pick the new a priori first, the base a priori second.\n');
                input('Press ENTER to continue', 's');
                new_dir = uigetdir(workdir, 'Pick the new a priori');
                old_dir = uigetdir(workdir, 'Pick the base a priori');
                
                if isnumeric(new_dir) || isnumeric(old_dir)
                    E.userCancel;
                elseif isempty(regexp(new_dir,'Atlanta','once')) || isempty(regexp(old_dir, 'Atlanta', 'once'))
                    fprintf('One of the directories does not appear to be a pseudo-retrieval directory (does not contain "Atlanta" in the name)\n');
                    if strcmpi(ask_multichoice('Continue?', {'y','n'}, 'default', 'n'),'n')
                        E.userCancel;
                    end
                end
                
                F_new = dir(fullfile(new_dir,'OMI*.mat'));
                F_old = dir(fullfile(old_dir,'OMI*.mat'));
        end
        if ~strcmpi(diff_mode,'all') && numel(F_new) ~= numel(F_old)
            E.callError('unequal_num_files','There are not equal numbers of hourly and hybrid files');
        end
        
        if isDisplay
            wb = waitbar(0,'Loading files');
        end
        
        if strcmpi(diff_mode,'all')
            F_hr = files_in_dates(F_hr,start_date,end_date);
            F_hy = files_in_dates(F_hy,start_date,end_date);
            F_mn = files_in_dates(F_mn,start_date,end_date);
            n = numel(F_hy);
        else
            F_new = files_in_dates(F_new,start_date,end_date);
            F_old = files_in_dates(F_old,start_date,end_date);
            n = numel(F_new);
        end
        
        for a=1:n
            if isDisplay
                waitbar(a/n);
            end
            if strcmpi(diff_mode,'all')
                hr = load(fullfile(hr_dir, F_hr(a).name));
                hy = load(fullfile(hy_dir, F_hy(a).name));
                mn = load(fullfile(mn_dir, F_mn(a).name));
                if a == 1
                    delmat = nan(numel(hr.Data.Longitude), n,3);
                end
                del = (hr.Data.BEHRAMFTrop ./ mn.Data.BEHRAMFTrop - 1)*100;
                delmat(:,a,1) = del(:);
                del = (hy.Data.BEHRAMFTrop ./ mn.Data.BEHRAMFTrop - 1)*100;
                delmat(:,a,2) = del(:);
                del = (hr.Data.BEHRAMFTrop ./ hy.Data.BEHRAMFTrop - 1)*100;
                delmat(:,a,3) = del(:);
            else
                new = load(fullfile(new_dir, F_new(a).name));
                base = load(fullfile(old_dir, F_old(a).name));
                if ~isempty(regexp(plot_mode,'pcolor-apri','once'))
                    new_prof = new.Data.BEHRNO2apriori;
                    base_prof = base.Data.BEHRNO2apriori;
                    if strcmpi(plot_mode,'pcolor-apri-stdp')
                        std_p = BEHR_std_pres;
                        new_pres = new.Data.BEHRPressureLevels;
                        new_prof_tmp = nan(size(new_prof)-[2 0 0]);
                        base_pres = base.Data.BEHRPressureLevels;
                        base_prof_tmp = nan(size(base_prof)-[2 0 0]);
                        for b=1:numel(new.Data.Longitude)
                            xx = ismember(new_pres(:,b),std_p);
                            new_prof_tmp(:,b) = new_prof(xx,b);
                            xx = ismember(base_pres(:,b),std_p);
                            base_prof_tmp(:,b) = base_prof(xx,b);
                        end
                        new_prof = new_prof_tmp;
                        base_prof = base_prof_tmp;
                    end
                    del = new_prof - base_prof;
                    del = squeeze(nansum2(del,1));
                else
                    del = (new.Data.BEHRAMFTrop ./ base.Data.BEHRAMFTrop - 1)*100;
                end
                if a == 1
                    delmat = nan(numel(del), numel(F_new));
                    if strcmp(plot_mode, 'dist')
                        distmat = nan(numel(del), numel(F_new));
                        anglemat = nan(numel(del), numel(F_new));
                        dvec = nan(1, numel(F_new));
                    end
                end
                delmat(:,a) = del(:);
                
                if ismember(plot_mode, {'dist','combo'}) || ~isempty(regexp(plot_mode,'scatter', 'once'));
                    for b=1:numel(del)
                        distmat(b,a) = m_lldist([new.Data.Longitude(b), city_lon], [new.Data.Latitude(b), city_lat]);
                        anglemat(b,a) = atan2d(new.Data.Latitude(b) - city_lat, new.Data.Longitude(b) - city_lon);
                    end
                end
                [s,e] = regexp(F_new(a).name,'\d\d\d\d\d\d\d\d');
                dvec(a) = datenum(F_new(a).name(s:e),'yyyymmdd');
            end
        end
        
        if isDisplay
            close(wb);
        end
        
        
        
        figure;
        if strcmp(plot_mode,'box')
            boxplot(delmat);
            xlab = get(gca,'xtick');
            set(gca,'xtick',[xlab(1), xlab(end)]);
            set(gca,'xticklabels',{datestr(dvec(1),'yyyy-mm-dd'), datestr(dvec(end),'yyyy-mm-dd')});
            ylabel('% difference AMF');
        elseif strcmp(plot_mode, 'dist')
            dvec = repmat(dvec,size(delmat,1),1);
            scatter(dvec(:), distmat(:), 16, delmat(:));
        elseif strcmp(plot_mode, 'combo');
            [distlabels, I] = sort(distmat(:,1));
            boxplot(delmat(I,:)');
            xlab = get(gca,'xtick');
            set(gca,'xtick',[xlab(1), xlab(end)]);
            set(gca,'xticklabels',[distlabels(1), distlabels(end)]);
            xlabel('Pixel dist. from Atlanta (km)');
            ylabel('% difference AMF');
        elseif ~isempty(regexp(plot_mode,'scatter', 'once'))
            if ~isempty(regexp(plot_mode,'dist','once'))
                scatter(distmat(:), delmat(:), 16, 'k');
                xlabel('Distance from Atlanta (km)')
            elseif ~isempty(regexp(plot_mode, 'angle', 'once'))
                scatter(anglemat(:), delmat(:), 16, 'k');
                xlabel('Heading from Atlanta (deg. CCW from east)');
            end
            ylabel('% difference full vs. hybrid');
            if ~isempty(regexp(plot_mode,'wbox', 'once'))
                yl = get(gca,'ylim');
                figure;
                boxplot(delmat(:));
                ylim(yl);
            end
        elseif ~isempty(regexp(plot_mode, 'pcolor', 'once'))
            if ~isempty(regexp(plot_mode,'pcolor-apri','once'))
                cbstr1 = 'column sum abs difference apriori';
                cbstr2 = '1\sigma sum \Delta apriori';
            else
                cbstr1 = 'Mean %\Delta AMF';
                cbstr2 = '1\sigma %\Delta AMF';
            end
            if strcmpi(plot_mode,'pcolor-med')
                delmatmean = reshape(nanmedian(delmat,2), size(new.Data.Longitude));
            else
                delmatmean = reshape(nanmean(delmat,2), size(new.Data.Longitude));
            end
            delmatstd = reshape(nanstd(delmat,0,2), size(new.Data.Longitude));
            pcolor(squeeze(new.Data.Loncorn(1,:,:)),squeeze(new.Data.Latcorn(1,:,:)),delmatmean);
            set(gca,'fontsize',20);
            cb=colorbar;
            cb.Label.String = cbstr1;
            cb.FontSize = 20;
            cmax = max(abs(cb.Limits));
            caxis([-cmax, cmax]);
            colormap(blue_red_cmap);
            
            figure; 
            pcolor(squeeze(new.Data.Loncorn(1,:,:)),squeeze(new.Data.Latcorn(1,:,:)),delmatstd);
            cb=colorbar;
            cb.Label.String = cbstr2;
            cb.FontSize = 20;
            caxis([0 max(cb.Limits)]);
            colormap('jet')
        elseif strcmpi(plot_mode, 'box-comp')
            delmat_final = reshape(delmat,[],3);
            boxplot(delmat_final)
            set(gca,'xticklabels',{'Full daily vs. monthly','Hybrid daily vs. monthly','Full vs. hybrid daily'})
            line(1:3, nanmean(delmat_final,1), 'color', 'k','marker','o','markersize',16,'linestyle','none');
        end
        set(gca,'fontsize',20)
    end

    function plot_behr_avg()
        allowed_apriori = {'hourly','hybrid','monthly'};
        try
            quantity = ask_multichoice('Which quantity do you want to plot?',{'amf','vcd'});
            new_apriori = ask_multichoice(sprintf('Which a priori do you want to use as the "new" apriori?\n'), allowed_apriori, 'default', 'hybrid');
            old_apriori = ask_multichoice(sprintf('Which a priori do you want to use as the base apriori?\n'), [allowed_apriori, {'none'}], 'default', 'monthly');
            clds = ask_multichoice('Filter by clouds?', {'y','n'});
            rowanom = ask_multichoice('Filter for row anomaly (XFlag > 0 or VCD > 1e17)?', {'y','n'});
            vcdqual = ask_multichoice('Filter for bad VCDs (vcdFlag odd or VCD < 0)?', {'y','n'});
            aw = ask_multichoice('Weight by areaweight?', {'y','n'});
            start_date = ask_date('Enter the date for the beginning of the averaging period (yyyy-mm-dd)');
            end_date = ask_date('Enter the end of the averaging period (yyyy-mm-dd)');
        catch err
            if strcmp(err.identifier, 'ask_multichoice:user_cancel')
                return
            else
                rethrow(err)
            end
        end
        
        sdate = datenum(start_date);
        edate = datenum(end_date);
        cldbool = strcmp(clds,'y');
        rowbool = strcmp(rowanom,'y');
        columnbool = strcmp(vcdqual,'y');
        awbool = strcmp(aw,'y');
        
        xl = [-87.1 -81.9];
        yl = [31.9 35.6];
        
        hr_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hourly - No ghost';
        hy_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hybrid - No ghost';
        mn_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Monthly - No ghost';
        switch new_apriori
            case 'hourly'
                np = hr_path;
            case 'hybrid'
                np = hy_path;
            case 'monthly'
                np = mn_path;
        end
        switch old_apriori
            case 'hourly'
                op = hr_path;
            case 'hybrid'
                op = hy_path;
            case 'monthly'
                op = mn_path;
        end
        F_new = dir(fullfile(np,'OMI*.mat'));
        if ~strcmpi(old_apriori,'none')
            F_old = dir(fullfile(op,'OMI*.mat'));
            if numel(F_new) ~= numel(F_old)
                E.callError('wrong_nfiles','There are different numbers of old and new files')
            end
        else 
            F_old = [];
        end
        
        
        n = numel(F_old) + numel(F_new);
        
        datamat_new = [];
        cldmat_new = [];
        awmat_new = [];
        rowmat_new = [];
        vcdmat_new = [];
        first_time = true;
        
        if isDisplay
            wb = waitbar(0,'Concatentaion progress');
        end
        for a=1:numel(F_new)
            if isDisplay; waitbar(a/n); end
            [s,e] = regexp(F_new(a).name,'\d\d\d\d\d\d\d\d');
            fdate = datenum(F_new(a).name(s:e),'yyyymmdd');
            if fdate >= sdate && fdate <= edate
                N=load(fullfile(np,F_new(a).name));
                OMI = N.OMI;
                
                if first_time
                    xx = any(OMI(1).Longitude >= min(xl) & OMI(1).Longitude <= max(xl) & OMI(1).Latitude >= min(yl) & OMI(1).Latitude <= max(yl),2);
                    yy = any(OMI(1).Longitude >= min(xl) & OMI(1).Longitude <= max(xl) & OMI(1).Latitude >= min(yl) & OMI(1).Latitude <= max(yl),1);
                end
                if strcmpi(quantity,'amf')
                    this_data = cat(3,OMI.BEHRAMFTrop);
                    this_data(this_data < 1e-4) = nan;
                elseif strcmpi(quantity,'vcd')
                    this_data = cat(3,OMI.BEHRColumnAmountNO2Trop);
                end
                datamat_new = cat(3, datamat_new, this_data(xx,yy,:));
                if cldbool
                    this_cld = cat(3, OMI.CloudFraction);
                    cldmat_new = cat(3,cldmat_new,this_cld(xx,yy,:));
                end
                if awbool
                    this_aw = cat(3, OMI.Areaweight);
                    awmat_new = cat(3,awmat_new,this_aw(xx,yy,:));
                end
                
                if rowbool || columnbool
                    this_colmat = cat(3, OMI.BEHRColumnAmountNO2Trop);
                end
                if rowbool
                    this_xflags = cat(3, OMI.XTrackQualityFlags);
                    this_rowmat = false(size(this_xflags));
                    for b=1:numel(this_xflags)
                        this_rowmat(b) = any(this_xflags{b}>0) || this_colmat(b) > 1e17;
                    end
                    rowmat_new = cat(3, rowmat_new,this_rowmat(xx,yy,:));
                end
                if columnbool
                    this_vcdflags = cat(3, OMI.vcdQualityFlags);
                    this_vcdmat = false(size(this_vcdflags));
                    for b=1:numel(this_vcdflags)
                        this_vcdmat(b) = any(mod([this_vcdflags{b}],2)~=0) || this_colmat(b) < 0;
                    end
                    vcdmat_new = cat(3, vcdmat_new, this_vcdmat(xx,yy,:));
                end
            end
        end
        
        if strcmpi(old_apriori,'none')
            datamat_old = zeros(size(datamat_new));
            cldmat_old = zeros(size(cldmat_new));
            awmat_old = ones(size(awmat_new)); % prevent divide by zero error
            rowmat_old = zeros(size(rowmat_new));
            vcdmat_old = zeros(size(vcdmat_new));
        else
            datamat_old = [];
            cldmat_old = [];
            awmat_old = [];
            rowmat_old = [];
            vcdmat_old = [];
            for a=1:numel(F_old)
                if isDisplay; waitbar((a+numel(F_new))/n); end
                [s,e] = regexp(F_old(a).name,'\d\d\d\d\d\d\d\d');
                fdate = datenum(F_old(a).name(s:e),'yyyymmdd');
                if fdate >= sdate && fdate <= edate
                    N=load(fullfile(op,F_old(a).name));
                    OMI = N.OMI;
                    if strcmpi(quantity,'amf')
                        this_data = cat(3,OMI.BEHRAMFTrop);
                        this_data(this_data < 1e-4) = nan;
                    elseif strcmpi(quantity,'vcd')
                        this_data = cat(3,OMI.BEHRColumnAmountNO2Trop);
                    end
                    datamat_old = cat(3, datamat_old, this_data(xx,yy,:));
                    if cldbool
                        this_cld = cat(3, OMI.CloudFraction);
                        cldmat_old = cat(3,cldmat_old,this_cld(xx,yy,:));
                    end
                    if awbool
                        this_aw = cat(3, OMI.Areaweight);
                        awmat_old = cat(3,awmat_old,this_aw(xx,yy,:));
                    end
                    if rowbool || columnbool
                        this_colmat = cat(3, OMI.BEHRColumnAmountNO2Trop);
                    end
                    if rowbool
                        this_xflags = cat(3, OMI.XTrackQualityFlags);
                        this_rowmat = false(size(this_xflags));
                        for b=1:numel(this_xflags)
                            this_rowmat(b) = any(this_xflags{b}>0) || this_colmat(b) > 1e17;
                        end
                        rowmat_old = cat(3, rowmat_old,this_rowmat(xx,yy,:));
                    end
                    if columnbool
                        this_vcdflags = cat(3, OMI.vcdQualityFlags);
                        this_vcdmat = false(size(this_vcdflags));
                        for b=1:numel(this_vcdflags)
                            this_vcdmat(b) = any(mod([this_vcdflags{b}],2)~=0) || this_colmat(b) < 0;
                        end
                        vcdmat_old = cat(3, vcdmat_old, this_vcdmat(xx,yy,:));
                    end
                end
            end
        end
        if isDisplay
            close(wb);
        end
        
        if cldbool
            datamat_new(cldmat_new > 0.2) = nan;
            datamat_old(cldmat_old > 0.2) = nan;
        end
        if rowbool
            datamat_new(logical(rowmat_new)) = nan;
            datamat_old(logical(rowmat_old)) = nan;
        end
        if columnbool
            datamat_new(logical(vcdmat_new)) = nan;
            datamat_old(logical(vcdmat_old)) = nan;
        end
        if awbool
            amfmean_new = nansum2(datamat_new .* awmat_new, 3) ./ nansum2(awmat_new, 3);
            amfmean_old = nansum2(datamat_old .* awmat_old, 3) ./ nansum2(awmat_old, 3);
        else
            amfmean_new = nanmean(datamat_new, 3);
            amfmean_old = nanmean(datamat_old, 3);
        end
        
        lon = OMI(1).Longitude(xx,yy);
        lat = OMI(1).Latitude(xx,yy);
        
        figure;
        pcolor(lon, lat, (amfmean_new ./ amfmean_old - 1)*100);
        shading flat
        cb=colorbar;
        cmax = max(abs(cb.Limits));
        caxis([-cmax cmax]);
        colormap('jet');
        title('Percent diff');
        
        figure;
        pcolor(lon, lat, amfmean_new);
        shading flat
        cb=colorbar;
        colormap('jet');
        title('New');
        
        figure;
        pcolor(lon, lat, amfmean_old);
        shading flat
        cb=colorbar;
        colormap('jet');
        title('Old');
    end

    function [new_apriori, base_apriori, new_pres, base_pres, lon, lat] = plot_pseudo_apriori(new_apriori, base_apriori)
        plottype = ask_multichoice('What type of plot to make?',{'slice','profshape'});
        shape_bool = ask_multichoice('Convert to shape factor?',{'y','n'});
        shape_bool = strcmpi(shape_bool,'y');
        if strcmpi(plottype, 'profshape')
            while true
                indstr = input('Enter the indicies to plot the profiles for: ', 's');
                indcell = strsplit(indstr,' ');
                if numel(indcell)==2
                    inds = zeros(size(indcell));
                    for a=1:numel(inds)
                        inds(a) = str2double(indcell{a});
                    end
                    break
                else
                    fprintf('Must be two numbers separated by a space. Try again.\n');
                end
            end
        end
        
        if nargin < 1
            allowed_apriori={'hourly','hybrid','hybrid-avg','monthly'};
            new_case = ask_multichoice('Which apriori is the new case?',allowed_apriori,'default','hybrid');
            base_case = ask_multichoice('Which apriori is the base case?',allowed_apriori,'default','monthly');
            std_p_bool = ask_multichoice('Use only the standard pressures?',{'y','n'});
            std_p_bool = strcmpi(std_p_bool,'y');
            start_date = ask_date('Enter the start date');
            end_date = ask_date('Enter the end date');
        end
        hourly_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hourly - No clouds - No ghost';
        hybrid_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hybrid - No clouds - No ghost';
        avg_hybrid_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Avg Hybrid - No clouds - No ghost';
        monthly_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Monthly - No clouds - No ghost';
        
        switch new_case
            case 'hourly'
                new_path = hourly_path;
            case 'hybrid'
                new_path = hybrid_path;
            case 'hybrid-avg'
                new_path = avg_hybrid_path;
            case 'monthly'
                new_path = monthly_path;
        end
        switch base_case
            case 'hourly'
                base_path = hourly_path;
            case 'hybrid'
                base_path = hybrid_path;
            case 'hybrid-avg'
                base_path = avg_hybrid_path;
            case 'monthly'
                base_path = monthly_path;
        end
        
        if nargin < 1
            [new_apriori, new_pres, loncorn, latcorn, new_terpres] = concat_files(new_path,'OMI*.mat',start_date,end_date,{'Data','BEHRNO2apriori'; 'Data','BEHRPressureLevels'; 'Data','Loncorn'; 'Data','Latcorn'; 'Data','GLOBETerpres'});
            [base_apriori, base_pres, base_terpres] = concat_files(base_path,'OMI*.mat',start_date,end_date,{'Data','BEHRNO2apriori'; 'Data', 'BEHRPressureLevels'; 'Data','GLOBETerpres'});
        end
        
        if std_p_bool
            [new_apriori, new_pres] = remove_interp_pres(new_apriori, new_pres);
            [base_apriori, base_pres] = remove_interp_pres(base_apriori, base_pres);
        end
        
        if shape_bool
            sz = size(new_apriori);
            for a=1:prod(sz(2:end))
                vcd = integPr2(new_apriori(:,a), new_pres(:,a), new_terpres(a)); 
                new_apriori(:,a) = new_apriori(:,a) / vcd;
            end
            for a=1:prod(sz(2:end))
                vcd = integPr2(base_apriori(:,a), base_pres(:,a), base_terpres(a)); 
                base_apriori(:,a) = base_apriori(:,a) / vcd;
            end
        end
        
        % use latcorn to get the pixels in the right place when plotting
        % with pcolor.
        lon = squeeze(loncorn(1,:,:,1));
        lat = squeeze(latcorn(1,:,:,1));
        if strcmpi(plottype,'slice')
            

            mean_new = permute(nanmean(new_apriori,4),[2 3 1]);
            mean_base = permute(nanmean(base_apriori,4),[2 3 1]);

            plot_slice_gui(mean_new - mean_base, lon, lat, 'New - base');
            plot_slice_gui(reldiff(mean_new, mean_base)*100, lon, lat, 'Per. diff. new - base')
        else
            profs_new = squeeze(new_apriori(:,inds(1),inds(2),:));
            pres_new = squeeze(new_pres(:,inds(1),inds(2),:));
            profs_base = squeeze(base_apriori(:,inds(1),inds(2),:));
            pres_base = squeeze(base_pres(:,inds(1),inds(2),:));
            
            figure;
            for a=1:size(profs_new,2)
                l_new=line(profs_new(:,a),pres_new(:,a),'color',[0.5 0.5 0.5]);
            end
            for a=1:size(profs_base,2)
                l_base=line(profs_base(:,a),pres_base(:,a),'color','r','linewidth',1);
            end
            l_new_mean = line(nanmean(profs_new,2), nanmean(pres_new,2), 'color', 'k', 'linewidth', 2, 'linestyle', '--');
            l_base_mean = line(nanmean(profs_base,2), nanmean(pres_base,2), 'color', 'y', 'linewidth', 2, 'linestyle', '--');
            set(gca,'ydir','reverse','fontsize',18);
            legend([l_new; l_new_mean; l_base; l_base_mean],{'New profiles','Mean new','Base profiles','Mean base'})
            title(sprintf('Profiles at %s',mat2str(inds)));
        end
    end

    function plot_cloudfrac(date_in)
        fpath = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hybrid - No ghost';
        fname = sprintf('OMI_BEHR_%s.mat', datestr(date_in,'yyyymmdd'));
        X = load(fullfile(fpath, fname),'OMI');
        xx = 761:860;
        yy = 141:210;
        OMI = X.OMI(2);
        figure;
        subplot(3,1,1);
        pcolor(OMI.Longitude(yy,xx), OMI.Latitude(yy,xx), OMI.CloudFraction(yy,xx));
        colorbar;
        caxis([0 1]);
        colormap('jet');
        shading flat
        title('OMI');
        
        subplot(3,1,2);
        pcolor(OMI.Longitude(yy,xx), OMI.Latitude(yy,xx), OMI.MODISCloud(yy,xx));
        colorbar;
        caxis([0 1]);
        colormap('jet');
        shading flat
        title('MODIS');
        
        OMI = omi_pixel_reject(OMI,'omi',1,'XTrackFlags');
        subplot(3,1,3);
        pcolor(OMI.Longitude(yy,xx), OMI.Latitude(yy,xx), double(OMI.Areaweight(yy,xx) == 0));
        colorbar;
        caxis([0 1]);
        colormap('jet');
        shading flat
        title('Other reject');
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
            
            colormap(blue_red_cmap);
            
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

    function plot_behr_apriori_surface()
        [filename, pathname] = uigetfile();
        D=load(fullfile(pathname,filename),'Data');
        Data = D.Data;
        nswaths = numel(Data);
        swath = ask_number(sprintf('Which swath to plot (1-%d)?', nswaths), 'testfxn', @(x) x > 0 && x <= nswaths, 'testmsg', sprintf('Value must be between 1 and %d', nswaths));
        
        lon = squeeze(Data(swath).Loncorn(1,:,:));
        lat = squeeze(Data(swath).Latcorn(1,:,:));
        apri_surf = behr_apriori_surface(Data(swath));
        
        figure; pcolor(lon,lat,apri_surf);
        cb=colorbar;
        cb.Label.String = '[NO_2] (mixing ratio)';
    end

    function avg_prof = plot_ens_apriori()
        prof_type = ask_multichoice('Plot concentrations or shape factors?',{'conc','shape'},'default','conc');
        sw_bool = ask_multichoice('Convolve with scattering weights?',{'y','n'});
        sw_bool = strcmpi(sw_bool,'y');
        ret_type = ask_multichoice('Which retrieval?',{'full','pseudo'});
        if strcmpi(ret_type,'full')
            allowed_apriori = {'monthly','hybrid','hourly'};
        elseif strcmpi(ret_type,'pseudo')
            allowed_apriori = {'monthly','hybrid','hourly','hybrid-avg'};
        end
        stdP_only_bool = ask_multichoice('Use only the standard pressures?',{'y','n'});
        stdP_only_bool = strcmpi(stdP_only_bool,'y');
        apriori_type = ask_multichoice('Which apriori profiles?',allowed_apriori);
        options.quad_bool = ask_multichoice('Divide into quadrants?',{'y','n'});
        options.dist_limit = ask_number('Only use pixels within x km of Atlanta? 0 if no','default',0,'testfxn',@(x)x>=0,'testmsg','Value must be >= 0');
        options.size_lim_type = ask_multichoice('Limit pixel size by',{'vza','row','none'});
        switch options.size_lim_type
            case 'vza'
                options.lim_crit = ask_number('Enter the maximum VZA to use (0 <= x <= 90)','testfxn',@(x) x>=0 && x<=90);
            case 'row'
                options.lim_crit = ask_number('How many rows to exclude from the edge?','testfxn',@(x) x>0);
            case 'none'
                options.lim_crit = [];
        end
        start_date = ask_date('Enter the start date');
        end_date = ask_date('Enter the ending date');
        
        options.quad_bool = strcmpi(options.quad_bool,'y');
        
        omi_clds_vec = [];
        apriori_mat = [];
        pres_mat = [];
        sw_mat = [];
        terpres_vec = [];
        quadrant_vec = [];
        
        domain.x = [-87.1, -81.9];
        domain.y = [31.9, 35.6];
        
        options.center_lon = -84.39;
        options.center_lat = 33.775;
        
        switch ret_type
            case 'full'
                switch apriori_type
                    case 'monthly'
                        behr_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Monthly - No ghost';
                    case 'hybrid'
                        behr_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hybrid - No ghost';
                    case 'hourly'
                        behr_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hourly - No ghost';
                end
            case 'pseudo'
                switch apriori_type
                    case 'monthly'
                        behr_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Monthly - No clouds - No ghost';
                    case 'hybrid'
                        behr_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hybrid - No clouds - No ghost';
                    case 'hourly'
                        behr_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hourly - No clouds - No ghost';
                    case 'hybrid-avg'
                        behr_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Avg Hybrid - No clouds - No ghost';
                end
        end
        behr_files = dir(fullfile(behr_path,'OMI*.mat'));
        behr_files = files_in_dates(behr_files, start_date, end_date);
        
        % Will describe the filtering conditions
        titlestr = apriori_type;
        switch options.size_lim_type
            case 'vza'
                titlestr = sprintf('%s, VZA < %f',titlestr,options.lim_crit);
            case 'row'
                titlestr = sprintf('%s, excl. %g edge rows', titlestr,options.lim_crit);
        end
        if options.dist_limit > 0
            titlestr = sprintf('%s, w/i %g km', titlestr, options.dist_limit);
        end
        
        if isDisplay
            wb=waitbar(0,'Loading files');
        end
        for f=1:numel(behr_files)
            if isDisplay
                waitbar(f/numel(behr_files));
            end
            D = load(fullfile(behr_path,behr_files(f).name));
            Data = D.Data;
            for s=1:numel(Data)
                [in, quadrant_in] = subset_BEHR_pixels(Data(s), domain.x, domain.y, options);
                omi_clouds_in = Data(s).CloudFraction(in);
                apriori_in = Data(s).BEHRNO2apriori(:,in);
                pres_in = Data(s).BEHRPressureLevels(:,in);
                sw_in = Data(s).BEHRScatteringWeights(:,in);
                terpres_in = Data(s).GLOBETerpres(in);
                
                if stdP_only_bool
                    [apriori_in] = remove_interp_pres(apriori_in, pres_in);
                    [sw_in, pres_in] = remove_interp_pres(sw_in, pres_in);
                end
  
                omi_clds_vec = cat(1, omi_clds_vec, omi_clouds_in(:));
                apriori_mat = cat(2, apriori_mat, apriori_in);
                pres_mat = cat(2, pres_mat, pres_in);
                sw_mat = cat(2, sw_mat, sw_in);
                terpres_vec = cat(1, terpres_vec, terpres_in);
                quadrant_vec = cat(1, quadrant_vec, quadrant_in(:));
            end
        end
        
        if isDisplay
            close(wb)
        end
        
        if strcmpi(prof_type,'shape')
            if isDisplay
                wb=waitbar(0,'Calculating shape factor');
            end
            for a=1:size(apriori_mat,2)
                if isDisplay
                    waitbar(a/size(apriori_mat,2))
                end
                vcd = integPr2(apriori_mat(:,a), pres_mat(:,a), terpres_vec(a));
                apriori_mat(:,a) = apriori_mat(:,a) / vcd;
            end
            if isDisplay
                close(wb)
            end
        end
        
        if sw_bool
            apriori_mat = apriori_mat .* sw_mat;
        end
        
        if options.quad_bool
            q = 4;
            quad_names = {'NE','NW','SW','SE'};
        else
            q = 1;
            quad_names = {'All'};
        end
        
        nbins = 5;
        cldfn=cell(1,nbins);
        cld_ll = nan(1,nbins);
        cld_ul = nan(1,nbins);
        clddiv = 1.0/nbins;
        for a=0:nbins-1;
            for b=1:q
                figure;
                cld_ll(a+1) = a*clddiv;
                cld_ul(a+1) = (a+1)*clddiv;
                xx = omi_clds_vec >= cld_ll(a+1) & omi_clds_vec <= cld_ul(a+1) & quadrant_vec == b; % quadrant_in will be all ones if not sorting by quadrant.
                this_apriori = apriori_mat(:,xx);
                this_pres = pres_mat(:,xx);
                for c=1:size(this_pres,2)
                    line(this_apriori(:,c)*1e9, this_pres(:,c), 'color', [0.5 0.5 0.5]);
                end
                cldfn{a+1} = sprintf('cld%d',a*2);
                avg_prof.(cldfn{a+1}).(quad_names{b}).x = nanmean(this_apriori,2)*1e9;
                avg_prof.(cldfn{a+1}).(quad_names{b}).y = nanmean(this_pres,2);
                line(nanmean(this_apriori,2)*1e9, nanmean(this_pres,2), 'color','r','linewidth',3);
                switch prof_type
                    case 'shape'
                        xlabel('Shape factor');
                    otherwise
                        xlabel('[NO_2] (ppbv)');
                end
                ylabel('Pres (hPa)');
                set(gca,'fontsize',16);
                set(gca,'ydir','reverse');
                if options.quad_bool
                    title(sprintf('%s, between %.1f and %.1f (%s quadrant)',titlestr,cld_ll(a+1),cld_ul(a+1),quad_names{b}));
                else
                    title(sprintf('%s, between %.1f and %.1f',titlestr,cld_ll(a+1),cld_ul(a+1)));
                end
            end
        end
        
        for b=1:q
            figure;
            hold on
            lstr = cell(1,nbins);
            
            overall_avg_x = [];
            overall_avg_y = [];
            
            for a=1:nbins
                plot(avg_prof.(cldfn{a}).(quad_names{b}).x,avg_prof.(cldfn{a}).(quad_names{b}).y,'linewidth',2);
                overall_avg_x = cat(2, overall_avg_x, avg_prof.(cldfn{a}).(quad_names{b}).x);
                overall_avg_y = cat(2, overall_avg_y, avg_prof.(cldfn{a}).(quad_names{b}).y);
                lstr{a} = sprintf('Clds %.1g-%.1g',cld_ll(a),cld_ul(a));
            end
            legend(lstr{:});
            set(gca,'ydir','reverse');
            title(sprintf('%s, avg. prof (%s quadrant)',titlestr,quad_names{b}));
            
            figure;
            plot(overall_avg_x, overall_avg_y, 'k-', 'linewidth', 2);
            set(gca,'ydir','reverse');
            title(sprintf('%s, overall avg. prof (%s quadrant)',titlestr,quad_names{b}));
            avg_prof.Overall.(quad_names{b}).x = overall_avg_x;
            avg_prof.Overall.(quad_names{b}).y = overall_avg_y;
        end
    end
    
    function full_apriori_avg_diff(hycld, mncld, relbool)
        % Takes the structure output from the plot_full_apriori() function
        % and plots the difference in profile shape or concentration at
        % each level.
        if ~exist('relbool','var')
            relbool = false;
        end
        
        cld_fns = fieldnames(hycld);
        quad_fns = fieldnames(hycld.(cld_fns{1}));
        quad_avg_hy_x = cell(1,4);
        quad_avg_hy_y = cell(1,4);
        quad_avg_mn_x = cell(1,4);
        %quad_avg_mn_y = cell(1,4);
        
        for a=1:numel(quad_fns)
            figure;
            for b=1:numel(cld_fns)
                % store the average profiles for each cloud subdivision
                quad_avg_hy_x{a} = cat(2,quad_avg_hy_x{a},hycld.(cld_fns{b}).(quad_fns{a}).x);
                quad_avg_hy_y{a} = cat(2,quad_avg_hy_y{a},hycld.(cld_fns{b}).(quad_fns{a}).y);
                quad_avg_mn_x{a} = cat(2,quad_avg_mn_x{a},mncld.(cld_fns{b}).(quad_fns{a}).x);
                %quad_avg_mn_y{a} = cat(2,quad_avg_mn_y{a},mncld.(cld_fns{b}).(quad_fns{a}).y);
                
                subplot(1,numel(cld_fns),b);
                hy_x = nanmean(hycld.(cld_fns{b}).(quad_fns{a}).x,2);
                mn_x = nanmean(mncld.(cld_fns{b}).(quad_fns{a}).x,2);
                y = hycld.(cld_fns{b}).(quad_fns{a}).y;
                if relbool
                    del = reldiff(hy_x, mn_x)*100;
                else
                    del = hy_x - mn_x;
                end
                plot(del, y, 'k','linewidth',2);
                set(gca,'ydir','reverse')
                title(cld_fns{b});
            end
            suptitle(quad_fns{a});
            
            % plot the overall average, if it's not in the structure
            % already
            if ~ismember('Overall',cld_fns)
                hy_x = nanmean(quad_avg_hy_x{a},2);
                hy_y = nanmean(quad_avg_hy_y{a},2);
                mn_x = nanmean(quad_avg_mn_x{a},2);
                if relbool
                    del = reldiff(hy_x, mn_x);
                    xstr = 'Rel. difference';
                else
                    del = hy_x - mn_x;
                    xstr = 'Abs. difference';
                end
                figure;
                plot(del, hy_y, 'k-','linewidth',2);
                set(gca,'ydir','reverse')
                xlabel(xstr);
                ylabel('Pressure (hPa)');
                title(sprintf('%s - avg over all cldfrac',quad_fns{a}));
            end
        end
    end

    function del_shape_vs_vcd()
        % Get the paths for the desired retrievals
        [new_path, base_path, ret_type, new_ret, base_ret] = return_paths();
        
        % For pseudo retrievals, we can plot the changes for a single
        % pixel, since it's a fixed grid.
        if strcmpi(ret_type,'pseudo')
            inds = input('Enter index values to plot for one pixel, or leave blank to subset: ','s');
        else
            inds = '';
        end
        
        if isempty(inds)
            % If not doing one pixel, give the user the chance to subset the
            % pixels.
            subset_bool = true;
            options = subset_BEHR_pixels();
            quadrant = ask_multichoice('Which quadrant to plot for? (This will override the previous quandrant answer.)',{'NE','NW','SW','SE','All','Don''t divide'},'list',true);
            if strcmpi(quadrant,'Don''t divide')
                options.quad_bool = false;
            else
                options.quad_bool = true;
            end
        else
            % Parse the indicies into a 2-element vector
            subset_bool = false;
            inds = str2double(strsplit(inds));
            if numel(inds) ~= 2
                E.badinput('Input should be two numbers separated by a space')
            end
        end
        
        % Plot a scatter plot of %change in shape factor vs. abs. change in
        % VCD for each level, or a histogram of those changes.
        plot_type = ask_multichoice('Which plot type?',{'scatter','hist','cumdist'});
        if strcmpi(plot_type,'scatter')
            allowed_vars = {'shape','vcd','conc','amf'};
            x_var = ask_multichoice(sprintf('Plot %%delta shape, delta VCD, %%delta [NO_2] or %%delta AMF? on x-axis?'),allowed_vars);
            y_var = ask_multichoice('And on the y-axis?', allowed_vars(~ismember(allowed_vars, x_var)));
        else
            x_var = 'vcd';
            y_var = 'shape';
        end
        start_date = ask_date('Enter the starting date');
        end_date = ask_date('Enter the ending date');
        
        domain_lonlim = [-87.1 -81.9];
        domain_latlim = [31.9 35.6];
        
        F_new = files_in_dates(dir(fullfile(new_path, 'OMI*.mat')),start_date,end_date);
        F_base = files_in_dates(dir(fullfile(base_path, 'OMI*.mat')),start_date,end_date);
        
        if numel(F_new) ~= numel(F_base)
            E.sizeMismatch('F_new','F_base');
        end
        
        apriori_new_mat = [];
        pres_new_mat = [];
        amf_new_vec = [];
        terpres_new_vec = [];
        
        apriori_base_mat = [];
        pres_base_mat = [];
        amf_base_vec = [];
        terpres_base_vec = [];
        
        quad_vec = [];
        
        for a=1:numel(F_new)
            N = load(fullfile(new_path,F_new(a).name),'Data');
            B = load(fullfile(base_path,F_base(a).name),'Data');
            if numel(N.Data) ~= numel(B.Data)
                E.sizeMismatch('N.Data','B.Data');
            end
            for s=1:numel(N.Data)
                apriori_new = N.Data(s).BEHRNO2apriori;
                pres_new = N.Data(s).BEHRPressureLevels;
                amf_new = N.Data(s).BEHRAMFTrop;
                terpres_new = N.Data(s).GLOBETerpres;
        
                [apriori_new, pres_new] = remove_interp_pres(apriori_new, pres_new);
                
                apriori_base = B.Data(s).BEHRNO2apriori;
                pres_base = B.Data(s).BEHRPressureLevels;
                amf_base = B.Data(s).BEHRAMFTrop;
                terpres_base = B.Data(s).GLOBETerpres;
        
                [apriori_base, pres_base] = remove_interp_pres(apriori_base, pres_base);
                
                % Assume that both the new and old Data struct have the
                % same pixels.
                if subset_bool
                    [in, quads] = subset_BEHR_pixels(N.Data, domain_lonlim, domain_latlim, options);
                    
                    apriori_new_mat = cat(2, apriori_new_mat, apriori_new(:,in));
                    pres_new_mat = cat(2, pres_new_mat, pres_new(:,in));
                    amf_new_vec = cat(1, amf_new_vec, amf_new(in));
                    terpres_new_vec = cat(1, terpres_new_vec, terpres_new(in));
                    
                    apriori_base_mat = cat(2, apriori_base_mat, apriori_base(:,in));
                    pres_base_mat = cat(2, pres_base_mat, pres_base(:,in));
                    amf_base_vec = cat(1, amf_base_vec, amf_base(in));
                    terpres_base_vec = cat(1, terpres_base_vec, terpres_base(in));
                    quad_vec = cat(1, quad_vec, quads(:));
                else
                    apriori_new_mat = cat(2, apriori_new_mat, apriori_new(:,inds(1),inds(2)));
                    pres_new_mat = cat(2, pres_new_mat, pres_new(:,inds(1),inds(2)));
                    amf_new_vec = cat(2, amf_new_vec, amf_new(inds(1),inds(2)));
                    terpres_new_vec = cat(2, terpres_new_vec, terpres_new(inds(1),inds(2)));
                    
                    apriori_base_mat = cat(2, apriori_base_mat, apriori_base(:,inds(1),inds(2)));
                    pres_base_mat = cat(2, pres_base_mat, pres_base(:,inds(1),inds(2)));
                    amf_base_vec = cat(2, amf_base_vec, amf_base(inds(1),inds(2)));
                    terpres_base_vec = cat(2, terpres_base_vec, terpres_base(inds(1),inds(2)));
                    quad_vec = cat(1, quad_vec,1); % if plotting one pixel, just make it a default quad number so that all of them are accepted during plotting
                end
                
            end
        end
        % Compute the VCDs, simultaneously converting to shape factors.
        vcd_new_vec = nan(size(terpres_new_vec));
        apriori_new_shape_mat = nan(size(apriori_new_mat));
        for a=1:numel(vcd_new_vec)
            vcd_new_vec(a) = integPr2(apriori_new_mat(:,a), pres_new_mat(:,a), terpres_new_vec(a));
            apriori_new_shape_mat(:,a) = apriori_new_mat(:,a) / vcd_new_vec(a);
        end
        vcd_base_vec = nan(size(terpres_base_vec));
        apriori_base_shape_mat = nan(size(apriori_new_mat));
        for a=1:numel(vcd_base_vec)
            vcd_base_vec(a) = integPr2(apriori_base_mat(:,a), pres_base_mat(:,a), terpres_base_vec(a));
            apriori_base_shape_mat(:,a) = apriori_base_mat(:,a) / vcd_base_vec(a);
        end
        
        % Set which quadrants to make plots for
        quad_names = {'NE','NW','SW','SE'};
        if ~subset_bool
            q1 = 1;
            q2 = 1;
        elseif strcmpi(quadrant,'All')
            q1 = 1;
            q2 = 4;
        elseif strcmpi(quadrant,'Don''t divide')
            % the quadrants vector will just be all 1s if the data wasn't
            % subset into quadrants.
            quad_names = {'All quadrants'};
            q1 = 1;
            q2 = 1;
        else
            q1 = find(strcmpi(quadrant,quad_names));
            q2 = q1;
        end
        
        std_p = BEHR_std_pres;
        for q=q1:q2
            xx = quad_vec == q;
            if any(ismember({'shape','conc'},{x_var,y_var}))
                n = size(apriori_new_mat,1);
            else
                n=1;
            end
            for a=1:n
                delvcd = vcd_new_vec(xx) - vcd_base_vec(xx);
                vcdstr = '\Delta a priori VCD';
                delconc = reldiff(apriori_new_mat(a,xx), apriori_base_mat(a,xx), true)*100;
                concstr = '%\Delta a priori [NO_2]';
                delamf = reldiff(amf_new_vec(xx), amf_base_vec(xx), true)*100;
                amfstr = '%\Delta AMF';
                delshape = reldiff(apriori_new_shape_mat(a,xx), apriori_base_shape_mat(a,xx),true)*100;
                shapestr = '%\Delta a priori shape factor';
                switch x_var
                    case 'vcd'
                        delx = delvcd;
                        xstr = vcdstr;
                    case 'conc'
                        delx = delconc;
                        xstr = concstr;
                    case 'amf'
                        delx = delamf;
                        xstr = amfstr;
                    case 'shape'
                        delx = delshape;
                        xstr = shapestr;
                end
                switch y_var
                    case 'vcd'
                        dely = delvcd;
                        ystr = vcdstr;
                    case 'conc'
                        dely = delconc;
                        ystr = concstr;
                    case 'amf'
                        dely = delamf;
                        ystr = amfstr;
                    case 'shape'
                        dely = delshape;
                        ystr = shapestr;
                end
                if all(isnan(dely)) || all(isnan(delx));
                    continue
                end
                figure;
                switch plot_type
                    case 'scatter'
                        scatter(delx, dely);
                        xlabel(xstr);
                        ylabel(ystr);
                        xl = get(gca,'xlim');
                        l(1) = line(xl, repmat(nanmean(dely),1,2), 'color','k','linewidth',2,'linestyle','--');
                        l(2) = line(xl, repmat(nanmedian(dely),1,2),'color','b','linewidth',2,'linestyle',':');
                        legend(l',{'Mean','Median'});
                    case 'hist'
                        hist(dely,20)
                        xlabel(ystr);
                        yl = get(gca,'ylim');
                        l = line(repmat(nanmean(dely),1,2),yl,'color','k','linewidth',2,'linestyle','--');
                        legend(l,'Mean');
                    case 'cumdist'
                        cumdist(dely,20,'color','b','linewidth',3)
                        xlabel(ystr);
                        yl = get(gca,'ylim');
                        l = line(repmat(nanmean(dely),1,2),yl,'color','k','linewidth',2,'linestyle','--');
                        legend(l,'Mean');
                end
                if n>1
                    title(sprintf('%s vs %s, %d hPa - %s', new_ret, base_ret, std_p(a), quad_names{q}));
                else
                    title(sprintf('%s vs %s - %s', new_ret, base_ret, quad_names{q}));
                end
            end
            if ismember(plot_type,{'hist','cumdist'})
                figure;
                switch plot_type
                    case 'hist'
                        hist(delx,20)
                    case 'cumdist'
                        cumdist(delx,20,'color','b','linewidth',3)
                end
                xlabel('Abs. \Delta VCD new - base');
                l = line(repmat(nanmean(delx),1,2),yl,'color','k','linewidth',2,'linestyle','--');
                legend(l,'Mean');
                title('Absdiff VCD');
            
            end
        end
        
        % Plot the pixel in question so we know where we are
        if ~subset_bool
            lon = squeeze(N.Data.Loncorn(1,:,:));
            lat = squeeze(N.Data.Latcorn(1,:,:));
            location = zeros(size(lon));
            location(inds(1), inds(2)) = 1;
            figure; pcolor(lon,lat,location);
            line(-84.39,33.775,'marker','p','markersize',20,'color','k','linestyle','none')
        end
        
    end

    function mag_delvcd_stats
        % This function will make plots representing how the frequency and
        % magnitude of changes around each of the three cities.
        %
        % Couple modes of operation: first, straight-up plot the
        % distribution of change values as boxplot or histogram. Second,
        % count the number of days where the change exceeds some value.
        % Third, count the days but make a polar plot of the percent of
        % days exceeding that value divided into quadrants or eighths, and
        % not counting days with no valid data in that quadrant.
        op_mode = ask_multichoice('Which plot to make?',{'boxplot','hist','count','polar'});
        [new_path, base_path] = return_paths;
        
        
        city_names = {'Atlanta','Birmingham','Montgomery'};
        city_lons = [   -84.39,...  %Atlanta
                        -86.80,...  %Birmingham
                        -86.30  ];  %Montgomery
        city_lats = [   33.775,...  %Atlanta
                        33.52,...   %Birmingham
                        32.37];     %Montgomery
        city_xl =   [   -87.1, -82;...
                        -88, -85.5;...
                        -87, -85];
        city_yl =   [   31.9, 35.5;...
                        33, 34.5;...
                        31.5, 33];
        
        % Subset BEHR options            
        opts.quad_bool = false;
        opts.dist_limit = 50; %km
        opts.size_lim_type = 'none';
        opts.lim_crit = [];
        % center_lon and _lat will be set for each city.
        
        %%%%% DATA PREP %%%%%
        
        if ismember(op_mode,{'boxplot','hist'})
            % For the boxplot and histogram, since we're looking at pixel
            % over all days, just concatenate the data.
            [lon, lat, loncorn, latcorn, cldfrac, vcdflag, xtrackflag, vcd_new, amf_new] = cat_sat_data(new_path,...
                {'Longitude','Latitude','Loncorn','Latcorn','CloudFraction','vcdQualityFlags','XTrackQualityFlags','BEHRColumnAmountNO2Trop','BEHRAMFTrop'},'prefix','OMI_BEHR','DEBUG_LEVEL','visual');
            [vcd_base, amf_base] = cat_sat_data(base_path,{'BEHRColumnAmountNO2Trop','BEHRAMFTrop'},'prefix','OMI_BEHR','DEBUG_LEVEL','visual');
            
            absdiff = vcd_new - vcd_base;
            DataNew.Longitude = lon;
            DataNew.Latitude = lat;
            DataNew.Loncorn = loncorn;
            DataNew.Latcorn = latcorn;
            DataNew.CloudFraction = cldfrac;
            DataNew.vcdQualityFlags = vcdflag;
            DataNew.XTrackQualityFlags = xtrackflag;
            DataNew.BEHRColumnAmountNO2Trop = vcd_new;
            
            DataBase.Longitude = lon;
            DataBase.Latitude = lat;
            DataBase.Loncorn = loncorn;
            DataBase.Latcorn = latcorn;
            DataBase.CloudFraction = cldfrac;
            DataBase.vcdQualityFlags = vcdflag;
            DataBase.XTrackQualityFlags = xtrackflag;
            DataBase.BEHRColumnAmountNO2Trop = vcd_base;
            
            badpix_new = find_bad_pixels(DataNew);
            badpix_new(amf_new == 1e-6) = true; % remove pixels that have the minimum value AMF
            badpix_base = find_bad_pixels(DataBase);
            badpix_base(amf_base == 1e-6) = true;
            
            absdiff(badpix_new | badpix_base) = nan;
            city_diffs = cell(size(city_names));
            for c=1:numel(city_names)
                opts.center_lon = city_lons(c);
                opts.center_lat = city_lats(c);
                in = subset_BEHR_pixels(DataNew, city_xl(c,:), city_yl(c,:), opts);
                city_diffs{c} = absdiff(in);
            end
        elseif strcmpi(op_mode,'count')
            % For count we need to load each day, reject any bad pixels,
            % subset the data to what we want, and see for how many days
            % are relatively clear and have changes > 1e15.
            dnums=datenum('2013-06-01'):datenum('2013-08-30');
            city_logs = cellmat(1,numel(city_names),1,numel(dnums)); % cell array of matrices filled with 0s
            city_logs_6e14 = cellmat(1,numel(city_names),1,numel(dnums)); % cell array of matrices filled with 0s
            city_logs_20percent = cellmat(1,numel(city_names),1,numel(dnums)); % cell array of matrices filled with 0s
            city_logs_quadsum = cellmat(1,numel(city_names),1,numel(dnums)); % cell array of matrices filled with 0s
            city_logs_quadsum2 = cellmat(1,numel(city_names),1,numel(dnums)); % cell array of matrices filled with 0s
            city_npix = cellmat(1,numel(city_names),1,numel(dnums)); % cell array of matrices filled with 0s
            city_npixclr = cellmat(1,numel(city_names),1,numel(dnums)); % cell array of matrices filled with 0s
            city_ndays_someclear = zeros(size(city_names));
            city_minchange = zeros(size(city_names));
            city_minperchange = zeros(size(city_names));
            city_maxchange = zeros(size(city_names));
            city_maxperchange = zeros(size(city_names));
            for d=1:numel(dnums)
                % concatenate each day's swaths
                prefix = sprintf('OMI_BEHR_%s',datestr(dnums(d),'yyyymmdd'));
                [lon, lat, loncorn, latcorn, cldfrac, vcdflag, xtrackflag, vcd_new, amf_new] = cat_sat_data(new_path,...
                    {'Longitude','Latitude','Loncorn','Latcorn','CloudFraction','vcdQualityFlags','XTrackQualityFlags','BEHRColumnAmountNO2Trop','BEHRAMFTrop'},'prefix',prefix,'DEBUG_LEVEL',0);
                [vcd_base, amf_base] = cat_sat_data(base_path,{'BEHRColumnAmountNO2Trop','BEHRAMFTrop'},'prefix',prefix,'DEBUG_LEVEL',0);
                
                DataNew.Longitude = lon;
                DataNew.Latitude = lat;
                DataNew.Loncorn = loncorn;
                DataNew.Latcorn = latcorn;
                DataNew.CloudFraction = cldfrac;
                DataNew.vcdQualityFlags = vcdflag;
                DataNew.XTrackQualityFlags = xtrackflag;
                DataNew.BEHRColumnAmountNO2Trop = vcd_new;
                
                DataBase.Longitude = lon;
                DataBase.Latitude = lat;
                DataBase.Loncorn = loncorn;
                DataBase.Latcorn = latcorn;
                DataBase.CloudFraction = cldfrac;
                DataBase.vcdQualityFlags = vcdflag;
                DataBase.XTrackQualityFlags = xtrackflag;
                DataBase.BEHRColumnAmountNO2Trop = vcd_base;
                
                badpix_new = find_bad_pixels(DataNew);
                badpix_new(amf_new == 1e-6) = true; % remove pixels that have the minimum value AMF
                badpix_base = find_bad_pixels(DataBase);
                badpix_base(amf_base == 1e-6) = true;
                
                absdiff = vcd_new - vcd_base;
                absdiff(badpix_new | badpix_base) = nan;
                perdiff = (vcd_new - vcd_base)./vcd_base * 100;
                perdiff(badpix_new | badpix_base) = nan;
                
                for c=1:numel(city_names)
                    opts.center_lon = city_lons(c);
                    opts.center_lat = city_lats(c);
                    in = subset_BEHR_pixels(DataNew, city_xl(c,:), city_yl(c,:), opts);
                    
                    if ~isempty(in)
                        % Are any pixels above the single pixel uncertainty for
                        % change in VCD?
                        city_logs{c}(d) = any(abs(absdiff(in)) >= 1e15);
                        city_logs_6e14{c}(d) = any(abs(absdiff(in)) >= 0.6e15); % boersma 04: DOAS uncertainty (0.4e15) plus slant column (0.45e15) added in quadrature
                        city_logs_20percent{c}(d) = any(abs(perdiff(in)) >= 20 & abs(absdiff(in)) >= 0.6e15); % bucsela 13: 20% estimated clear sky AMF error.
                        city_logs_quadsum{c}(d) = any(abs(absdiff(in)) >= sqrt( 0.7e15^2 + 0.2e15^2 + (0.20 * vcd_base(in)).^2 ));
                        city_logs_quadsum2{c}(d) = any(abs(absdiff(in)) >= sqrt( 1e15^2 + (0.25 * vcd_base(in)).^2 ));
                        % If all the pixels were bad, then this day is probably
                        % falling in the row anomaly, but it could be heavily
                        % clouded. In either case, it shouldn't be counted as
                        % part of the statistics b/c we can't see any pixels!
                        city_ndays_someclear(c) = city_ndays_someclear(c) + ~all(isnan(absdiff(in)));
                        % These will give us information about how many pixels
                        % were clear
                        city_npix{c}(d) = numel(in);
                        city_npixclr{c}(d) = sum(~isnan(absdiff(in)));
                        % Update the max and min changes
                        city_maxchange(c) = max([city_maxchange(c), max(absdiff(in))]);
                        m = absdiff == max(absdiff(in));
                        city_maxperchange(c) = max([city_maxperchange(c), perdiff(m)]);
                        city_minchange(c) = min([city_minchange(c), min(absdiff(in))]);
                        m = absdiff == min(absdiff(in));
                        city_minperchange(c) = min([city_minperchange(c), perdiff(m)]);
                    end
                end
            end
        end
        
        %%%%%% PLOTTING %%%%%
        
        switch op_mode
            case 'boxplot'
                % Case 1: make a boxplot for each city
                city_diff_for_box = padcat(2, city_diffs{:});
                figure;
                boxplot(city_diff_for_box, city_names);
                ylabel('\Delta VCD');
                set(gca,'fontsize',16)
            case 'hist'
                figure;
                for c=1:numel(city_names)
                    subplot(numel(city_names),1,c);
                    hist(city_diffs{c},20);
                    xlabel('\Delta VCD');
                    title(city_names{c});
                    set(gca,'fontsize',14);
                end
            case 'count'
                % making a table this time 
                tabcell = cell(numel(city_names), 6);
                varnames = cell(1,6);
                for c=1:numel(city_names)
                    fracclr = city_npixclr{c} ./ city_npix{c};
                    
                    varnames{1,1} = 'PercentDaysDVCDGT1e15';
                    tabcell{c,1} = sum(city_logs{c})/city_ndays_someclear(c)*100;
                    
                    varnames{1,2} = 'PercentDaysDVCDGT6e14';
                    tabcell{c,2} = sum(city_logs_6e14{c})/city_ndays_someclear(c)*100;
                    
                    varnames{1,3} = 'PercentDaysDVCDGT20percentOR1e15';
                    tabcell{c,3} = sum(city_logs_20percent{c})/city_ndays_someclear(c)*100;
                    
                    varnames{1,4} = 'PercentDaysDVCDGT20percentPLUS6e14';
                    tabcell{c,4} = sum(city_logs_quadsum{c})/city_ndays_someclear(c)*100;
                    
                    varnames{1,5} = 'PercentDaysDVCDGT25percentPLUS1e15';
                    tabcell{c,5} = sum(city_logs_quadsum2{c})/city_ndays_someclear(c)*100;
                    
                    varnames{1,6} = 'MinChange';
                    tabcell{c,6} = city_minchange(c);
                    varnames{1,7} = 'MinPercentChange';
                    tabcell{c,7} = city_minperchange(c);
                    varnames{1,8} = 'MaxChange';
                    tabcell{c,8} = city_maxchange(c);
                    varnames{1,9} = 'MaxPercentChange';
                    tabcell{c,9} = city_maxperchange(c);
                    
                    varnames{1,10} = 'NumDaysWithGT0percentClearPix';
                    tabcell{c,10} = city_ndays_someclear(c);
                    varnames{1,11} = 'FracDaysGT50percentClearChange';
                    tabcell{c,11} = sum(fracclr > 0.5 & city_logs{c})/sum(fracclr > 0.5);
                    varnames{1,12} = 'FracDaysGT80percentClearChange';
                    tabcell{c,12} = sum(fracclr > 0.8 & city_logs{c})/sum(fracclr > 0.8);
                    
                    dVCDStats.dVCD_logical.(city_names{c}) = city_logs{c};
                    dVCDStats.npix.(city_names{c}) = city_npix{c};
                    dVCDStats.npixclr.(city_names{c}) = city_npixclr{c};
                end
                
                dVCDStats.table = cell2table(tabcell,'VariableNames',varnames,'RowNames',city_names);
                fprintf('Placing output variable "dVCDStats" in the base workspace\n')
                putvar(dVCDStats);
        end
        
    end
    
    function pix_pos_vs_cld()
        % This function will plot the pixel position as a function of
        % distance or direction from Atlanta vs. cloud fraction and color
        % by % change in AMF or VCD.
        diff_type = ask_multichoice('Which diff type to use?',{'amf','vcd'},'default','vcd','list',true);
        max_dist = ask_number('How far in km from Atlanta should we consider pixels for?','testfxn',@(x) x>0);
        hybrid_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hourly - No ghost'; 
        monthly_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Monthly - No ghost'; 
        
        city_lon = -84.39;
        city_lat = 33.775;
        
        switch lower(diff_type)
            case 'vcd'
                behrfield = 'BEHRColumnAmountNO2Trop';
                cblabel = '%\Delta VCD';
            case 'amf'
                behrfield = 'BEHRAMFTrop';
                cblabel = '%\Delta AMF';
            otherwise
                E.notimplemented('No corresponding field for the diff type %s is available. Did you add that option recently and forget to update the fields?', diff_type);
        end
        
        pixdist = [];
        pixangle = [];
        clds = [];
        hyquant = [];
        mnquant = [];
        
        opts.dist_limit = max_dist;
        opts.quad_bool = 0;
        opts.size_lim_type = 'none';
        opts.lim_crit = [];
        opts.center_lon = city_lon;
        opts.center_lat = city_lat;
        
        if isDisplay
            wb = waitbar(0,'Loading clouds');
        end
        
        F = dir(fullfile(hybrid_path, 'OMI_BEHR*.mat'));
        
        for d=1:numel(F)
            if isDisplay
                waitbar(d/numel(F));
            end
                
            H = load(fullfile(hybrid_path,F(d).name),'Data');
            Data = H.Data;
            M = load(fullfile(monthly_path,F(d).name),'Data');
            
            for s=1:numel(Data)
                badpix = find_bad_pixels(Data(s),1);
                Data(s).CloudFraction(badpix) = nan;
                
                in = subset_BEHR_pixels(Data(s),[-90 -80],[30 40],opts);
                
                dlon = Data(s).Longitude(in) - city_lon;
                dlat = Data(s).Latitude(in) - city_lat;
                this_dist = sqrt(dlon.^2 + dlat.^2);
                this_angle = atan2d(dlat, dlon);
                
                pixdist = cat(1, pixdist, this_dist);
                pixangle = cat(1, pixangle, this_angle);
                clds = cat(1, clds, Data(s).CloudFraction(in));
                hyquant = cat(1, hyquant, Data(s).(behrfield)(in));
                mnquant = cat(1, mnquant, M.Data(s).(behrfield)(in));
            end
        end
        
        if isDisplay
            close(wb)
        end
        
        % I'm guessing these figures are going to be slow to draw b/c
        % they'll have a lot of points, so be ready for that.
        del = reldiff(hyquant, mnquant)*100;
        figure; 
        scatter(pixdist, clds, 32, del, 'filled');
        xlabel('Distance (degrees)')
        ylabel('Cloud fraction')
        cb = colorbar;
        cb.Label.String = cblabel;
        set(gca,'fontsize',16);
        
        figure; 
        scatter(pixangle, clds, 32, del, 'filled');
        xlabel('Angle (degrees CCW from east)')
        ylabel('Cloud fraction')
        cb = colorbar;
        cb.Label.String = cblabel;
        set(gca,'fontsize',16);
    end

    function wind_cond_vs_cld()
        % This function will plot wind speed and direction based off of the
        % wind conditions file for atlanta vs. cloud fraction.
        diff_type = ask_multichoice('Which diff type to use?',{'rel','abs'},'default','rel','list',true);
        hybrid_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hourly - No ghost'; 
        monthly_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Monthly - No ghost'; 
        W = load('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta-Wind-Conditions-1900UTC-5layers.mat');
        
        clds50 = nan(size(W.dnums));
        dAMF50 = nan(size(W.dnums));
        dVCD50 = nan(size(W.dnums));
        clds150 = nan(size(W.dnums));
        dAMF150 = nan(size(W.dnums));
        dVCD150 = nan(size(W.dnums));
        
        opts.quad_bool = 0;
        opts.size_lim_type = 'none';
        opts.lim_crit = [];
        opts.center_lon = -84.39;
        opts.center_lat = 33.775;
        
        if isDisplay
            wb = waitbar(0,'Loading clouds');
        end
        
        for d=1:numel(W.dnums)
            if isDisplay
                waitbar(d/numel(W.dnums));
            end
                
            behr_file = sprintf('OMI_BEHR_%s.mat',datestr(W.dnums(d),'yyyymmdd'));
            H = load(fullfile(hybrid_path,behr_file),'Data');
            Data = H.Data;
            M = load(fullfile(monthly_path,behr_file),'Data');
            tmp_clds50 = [];
            tmp_hyamf50 = [];
            tmp_mnamf50 = [];
            tmp_hyvcd50 = [];
            tmp_mnvcd50 = [];
            tmp_clds150 = [];
            tmp_hyamf150 = [];
            tmp_mnamf150 = [];
            tmp_hyvcd150 = [];
            tmp_mnvcd150 = [];
            for s=1:numel(Data)
                badpix = find_bad_pixels(Data(s),1);
                Data(s).CloudFraction(badpix) = nan;
                
                opts.dist_limit = 50;
                in50 = subset_BEHR_pixels(Data(s),[-90 -80],[30 40],opts);
                tmp_clds50 = cat(1, tmp_clds50, Data(s).CloudFraction(in50));
                tmp_hyamf50 = cat(1, tmp_hyamf50, Data(s).BEHRAMFTrop(in50));
                tmp_mnamf50 = cat(1, tmp_mnamf50, M.Data(s).BEHRAMFTrop(in50));
                tmp_hyvcd50 = cat(1, tmp_hyvcd50, Data(s).BEHRColumnAmountNO2Trop(in50));
                tmp_mnvcd50 = cat(1, tmp_mnvcd50, M.Data(s).BEHRColumnAmountNO2Trop(in50));
                
                
                opts.dist_limit = 150;
                in150 = subset_BEHR_pixels(Data(s),[-90, -80],[30 40],opts);
                tmp_clds150 = cat(1, tmp_clds150, Data(s).CloudFraction(in150));
                tmp_hyamf150 = cat(1, tmp_hyamf150, Data(s).BEHRAMFTrop(in150));
                tmp_mnamf150 = cat(1, tmp_mnamf150, M.Data(s).BEHRAMFTrop(in150));
                tmp_hyvcd150 = cat(1, tmp_hyvcd150, Data(s).BEHRColumnAmountNO2Trop(in150));
                tmp_mnvcd150 = cat(1, tmp_mnvcd150, M.Data(s).BEHRColumnAmountNO2Trop(in150));
            end
            
            clds50(d) = nanmean(tmp_clds50);
            clds150(d) = nanmean(tmp_clds150);
            if strcmpi(diff_type,'abs')
                dAMF50(d) = nanmean(tmp_hyamf50 - tmp_mnamf50);
                dAMF150(d) = nanmean(tmp_hyamf150 - tmp_mnamf150);
                dVCD50(d) = nanmean(tmp_hyvcd50 - tmp_mnvcd50);
                dVCD150(d) = nanmean(tmp_hyvcd150 - tmp_mnvcd150);
            elseif strcmpi(diff_type,'rel')
                dAMF50(d) = nanmean(reldiff(tmp_hyamf50, tmp_mnamf50))*100;
                dAMF150(d) = nanmean(reldiff(tmp_hyamf150, tmp_mnamf150))*100;
                dVCD50(d) = nanmean(reldiff(tmp_hyvcd50, tmp_mnvcd50))*100;
                dVCD150(d) = nanmean(reldiff(tmp_hyvcd150, tmp_mnvcd150))*100;
            else
                E.notimplemented('diff_type = %s',diff_type)
            end
        end
        
        if isDisplay
            close(wb)
        end
        % Make figures. For now we'll do two for each distance, until I
        % upgrade to 2016a and have access to polar plots. The two will be
        % cloud fraction vs. wind speed and direction.
        figure;
        scatter(W.windvel, clds50, 32, dAMF50, 'filled');
        xlabel('Wind speed (m/s)')
        ylabel('Cloud fraction')
        title('Within 50 km')
        cb=colorbar;
        cb.Label.String = '%\Delta AMF';
        colormap('jet')
        caxis([-50 50])
        
        figure;
        scatter(W.theta, clds50, 32, dAMF50, 'filled');
        xlabel('Wind direction (CCW from E)')
        ylabel('Cloud fraction')
        title('Within 50 km');
        cb=colorbar;
        cb.Label.String = '%\Delta AMF';
        colormap('jet')
        caxis([-50 50])
        
        figure;
        scatter(W.windvel, clds150, 32, dAMF50, 'filled');
        xlabel('Wind speed (m/s)')
        ylabel('Cloud fraction')
        title('Within 150 km')
        cb=colorbar;
        cb.Label.String = '%\Delta AMF';
        colormap('jet')
        caxis([-50 50])
        
        figure;
        scatter(W.theta, clds150, 32, dAMF50, 'filled');
        xlabel('Wind direction (CCW from E)')
        ylabel('Cloud fraction')
        title('Within 150 km');
        cb=colorbar;
        cb.Label.String = '%\Delta AMF';
        colormap('jet')
        caxis([-50 50])
        
        figure;
        scatter(W.windvel, clds50, 32, dVCD50, 'filled');
        xlabel('Wind speed (m/s)')
        ylabel('Cloud fraction')
        title('Within 50 km')
        cb=colorbar;
        cb.Label.String = '%\Delta VCD';
        colormap('jet')
        caxis([-50 50])
        
        figure;
        scatter(W.theta, clds50, 32, dVCD50, 'filled');
        xlabel('Wind direction (CCW from E)')
        ylabel('Cloud fraction')
        title('Within 50 km');
        cb=colorbar;
        cb.Label.String = '%\Delta VCD';
        colormap('jet')
        caxis([-50 50])
        
        figure;
        scatter(W.windvel, clds150, 32, dVCD150, 'filled');
        xlabel('Wind speed (m/s)')
        ylabel('Cloud fraction')
        title('Within 150 km')
        cb=colorbar;
        cb.Label.String = '%\Delta VCD';
        colormap('jet')
        caxis([-50 50])
        
        figure;
        scatter(W.theta, clds150, 32, dVCD150, 'filled');
        xlabel('Wind direction (CCW from E)')
        ylabel('Cloud fraction')
        title('Within 150 km');
        cb=colorbar;
        cb.Label.String = '%\Delta VCD';
        colormap('jet')
        caxis([-50 50])
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

    function [windspd, winddir] = calc_avg_wind(xlon, xlat, U, V, clon, clat)
        [xx1,yy1] = find_square_around(xlon, xlat, clon, clat, 1);
        windspd = nan(1, size(U,3));
        winddir = nan(1, size(U,3));
        for a=1:size(U,3)
            Ubar = nanmean(reshape(U(xx1, yy1, a),[],1));
            Vbar = nanmean(reshape(V(xx1, yy1, a),[],1));
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

    
    function avg_prof = cat_cloud_frac_bins(avg_prof)
        cldfn = fieldnames(avg_prof);
        quad_names = fieldnames(avg_prof.(cldfn{1}));
        for b=1:numel(quad_names)
            overall_avg_x = [];
            overall_avg_y = [];
            for a=1:numel(cldfn)
                overall_avg_x = cat(2, overall_avg_x, avg_prof.(cldfn{a}).(quad_names{b}).x);
                overall_avg_y = cat(2, overall_avg_y, avg_prof.(cldfn{a}).(quad_names{b}).y);
            end
            avg_prof.Overall.(quad_names{b}).x = overall_avg_x;
            avg_prof.Overall.(quad_names{b}).y = overall_avg_y;
        end
    end
end

function files = files_in_dates(files, start_date, end_date)
% Will remove files outside of the given date range. Dates can be as
% datenumbers or strings. Right now assumes the date in the filename is
% yyyymmdd or yyyy-mm-dd
E = JLLErrors;

file_dates = nan(size(files));
for a=1:numel(files)
    [s,e] = regexp(files(a).name,'\d\d\d\d\d\d\d\d');
    if ~isempty(s)
        file_dates(a) = datenum(files(a).name(s:e),'yyyymmdd');
    else
        [s,e] = regexp(files(a).name,'\d\d\d\d-\d\d-\d\d');
        if isempty(s)
            E.callError('unknown_date_format','Cannot find the file''s date in the file name')
        end
        file_dates(a) = datenum(files(a).name(s:e),'yyyy-mm-dd');
    end
end

sdate = datenum(start_date);
edate = datenum(end_date);

xx = file_dates >= sdate & file_dates <= edate;
files = files(xx);

end

function varargout = concat_files(fpath,fpattern,start_date,end_date,fields,fieldinds)
% Will handle file concatenation. fpath should be the path to the files and
% fpattern the pattern of the files to concatenate. Start and end date are
% the first and last date to concatenate. fields should be a cell array of
% the necessary fields as strings to get at the data in the files, and
% fieldinds the scalar indicies for any intervening fields. So if in each
% file, you want Data.TraceGases.NO2, fields should be
% {'Data','TraceGases','NO2'} and if you need to specify
% Data(2).TraceGases(1).NO2, fieldinds should be [2 1]. If you want to load
% multiple fields, specified each as a row of fields. If some rows have
% more elements, fill the shorter rows with empty spaces, but always make
% the final variable the last entry.
F = dir(fullfile(fpath,fpattern));
n = numel(F);

sdate = datenum(start_date);
edate = datenum(end_date);

E=JLLErrors;

if ~iscellstr(fields)
    E.badinput('fields must be a cell array of strings')
end

if ~exist('fieldinds','var')
    fieldinds = ones(size(fields)-[0 1]);
elseif numel(fieldinds) ~= (numel(fields)-1)
    E.numelMismatch('fields','fieldinds');
end

if isDisplay
    wb = waitbar(0,'Concatenation progress');
end

varargout = cell(1,nargout);
for a=1:n
    if isDisplay
        waitbar(a/n)
    end
    [s,e] = regexp(F(a).name,'\d\d\d\d\d\d\d\d');
    fdate = datenum(F(a).name(s:e),'yyyymmdd');
    if fdate >= sdate && fdate <= edate
        D = load(fullfile(fpath,F(a).name));
        for b=1:nargout
            var = D;
            for f=1:size(fields,2)
                if isempty(fields{b,f})
                    continue
                elseif f < size(fields,2)
                    var = var.(fields{b,f})(fieldinds(b,f));
                else
                    var = var.(fields{b,f});
                end
            end
            varargout{b} = cat(ndims(var)+1, varargout{b}, var);
        end
    end
end

if isDisplay
    close(wb)
end
end

function utchrs = parse_apriori_mode(ap_mode)
% Returns a matrix of the UTC hours defined in the BEHRaprioriMode field
% Pass the string from that field as the input or one element of the Data
% structure

E = JLLErrors;

if isstruct(ap_mode)
    if ~isscalar(ap_mode)
        E.badinput('Only pass one element of the Data structure if given as input');
    elseif ~isfield(ap_mode,'BEHRaprioriMode')
        E.badinput('Input structure does not have the field BEHRaprioriMode');
    else
        ap_mode = ap_mode.BEHRaprioriMode;
    end
elseif ~ischar(ap_mode)
    E.badinput('ap_mode must either be a scalar structure containing the field BEHRaprioriMode or the string contained in the field');
end

s = regexp(ap_mode,'[');
e = regexp(ap_mode,']');
utchrs = eval(ap_mode(s:e));

end

function [apriori_cut, pres_cut] = remove_interp_pres(apriori, pres)
% Removes non-standard pressure levels (added by interpolation during
% integPr2)
E=JLLErrors;
sz = size(apriori);
if any(sz ~= size(pres))
    E.sizeMismatch('apriori','pres')
end
std_P = BEHR_std_pres;
apriori_cut = nan([numel(std_P),sz(2:end)]);
pres_cut = nan([numel(std_P),sz(2:end)]);
for a=1:prod(sz(2:end))
    xx = ismember(pres(:,a),std_P);
    apriori_cut(:,a) = apriori(xx,a);
    pres_cut(:,a) = pres(xx,a);
end
end

function [new_path, base_path, ret_type, new_apriori, base_apriori] = return_paths()
ret_type = ask_multichoice('Which retrieval?',{'full','pseudo'});
if strcmpi(ret_type,'full')
    allowed_apriori = {'monthly','hybrid','hourly'};
elseif strcmpi(ret_type,'pseudo')
    allowed_apriori = {'monthly','hybrid','hourly','hybrid-avg'};
end

% so that there's something to return in any case
new_path = '';
base_path = '';
new_apriori = '';
base_apriori = '';

for a=1:min(nargout,2) % if no second output is requested, don't ask for a "new" and "base" case, only ask for one case.
    if a==1 && nargout == 1
        casestr = '?';
    elseif a==1
        casestr = ' for the new case?';
    elseif a==2
        casestr = ' for the base case?';
    end
    apriori_type = ask_multichoice(sprintf('Which apriori to use%s',casestr),allowed_apriori);
    switch ret_type
        case 'full'
            switch apriori_type
                case 'monthly'
                    this_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Monthly - No ghost - lw 13.5 overpass - 18-22 UTC';
                    %this_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Monthly - No ghost';
                case 'hybrid'
                    this_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hybrid - No ghost';
                    warning('The hybrid case has not been updated to include UTC 1800 profiles');
                case 'hourly'
                    this_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hourly - No ghost - lw 13.5 overpass - 18-22 UTC';
                    %this_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hourly - No ghost';
            end
        case 'pseudo'
            switch apriori_type
                case 'monthly'
                    this_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Monthly - No clouds - No ghost - lonweight 13.5 overpass - UTC 1800-2200';
                case 'hybrid'
                    %this_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hybrid - No clouds - No ghost';
                    this_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hybrid - No clouds - No ghost - lonweight 13.5 overpass - UTC 1800-2200';
                case 'hourly'
                    this_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Hourly - No clouds - No ghost - UTC 1800-2200';
                case 'hybrid-avg'
                    this_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta BEHR Avg Hybrid - No clouds - No ghost';
                    warning('The hybrid-avg case has not been updated to include UTC 1800 profiles');
            end
    end
    if a==1
        new_path = this_path;
        new_apriori = apriori_type;
    else
        base_path = this_path;
        base_apriori = apriori_type;
    end
end
end

function badpix = find_bad_pixels(Data, cldfrac)
if ~exist('cldfrac','var')
    cldfrac = 0.2;
end
badpix = Data.CloudFraction > cldfrac | Data.XTrackQualityFlags ~= 0 | mod(Data.vcdQualityFlags,2) ~= 0 | Data.BEHRColumnAmountNO2Trop < 0 | Data.BEHRColumnAmountNO2Trop > 1e17;
end

function [interp_val, interp_pres] = interp_to_std_p(vals, P)
behr_stdp = BEHR_std_pres;
vals_tmp = nan(numel(behr_stdp), size(vals,2));
for a=1:size(vals,2)
    vals_tmp(:,a) = interp1(P(:,a), vals(:,a), behr_stdp);
end
P = repmat(behr_stdp',1,size(vals,2));
vals = vals_tmp;
end

function [file_out] = return_data_file_info(source, city, timemode, res, apriori, date_in)
% File paths: need the daily and monthly paths. Will construct a
% cell array where the first dimension corresponds to the
% difference sources (e.g. WRF, BEHR, pseudo-BEHR) and the second
% to cities. Put NaNs for cases that don't exist.
%
% File name format: strings should include a %1$04d where the year,
% %2$02d where the month, and %3$02d where the day get filled in.
% Only expected to be different for the different sources.
E = JLLErrors;
sharedir = '/Volumes/share2/USERS/LaughnerJ/';
workdir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/';
if strcmpi(city,'SF')
    wrf_coast = 'W';
    behr_coast = 'W';
else
    wrf_coast = 'E';
    behr_coast = 'SE';
end
switch source
    case 'wrf'
        daily_file_name_spec = 'WRF_BEHR_hourly_%1$04d-%2$02d-%3$02d.nc';
        monthly_file_name_spec = 'WRF_BEHR_monthly_%1$04d-%2$02d-%3$02d.nc';
        switch res
            case 'f'
                daily_path = fullfile(sharedir,'WRF',sprintf('%s_US_BEHR',wrf_coast),'hourly-14.0-lonwt-1822UTC');
                hybrid_path = NaN;
                monthly_path = fullfile(sharedir,'WRF',sprintf('%s_US_BEHR',wrf_coast),'monthly-13.5-lonwt-1822UTC');
            case 'c'
                daily_path = fullfile(sharedir,'WRF',sprintf('%s_US_BEHR_COARSE_13lonwt_1822UTC',wrf_coast),'hourly');
                hybrid_path = NaN;
                monthly_path = fullfile(sharedir,'WRF',sprintf('%s_US_BEHR_COARSE_13lonwt_1822UTC',wrf_coast),'monthly');
                if strcmpi(apriori, 'hourly')
                    warning('The coarse hourly apriori using lonweight assuming OMI overpass at 1330 LST or profiles from 1800-2200 UTC has not been downloaded');
                end
        end
    case 'behr'
        daily_file_name_spec = 'OMI_BEHR_%1$04d%2$02d%3$02d.mat';
        monthly_file_name_spec = 'OMI_BEHR_%1$04d%2$02d%3$02d.mat';
        s = 2;
        switch res
            case 'f'
                daily_path = fullfile(workdir,sprintf('%s US BEHR Hourly - No ghost - lw 13.5 overpass - 18-22 UTC',behr_coast));
                hybrid_path = fullfile(workdir,sprintf('%s US BEHR Hybrid - No ghost',behr_coast));
                monthly_path = fullfile(workdir,sprintf('%s US BEHR Monthly - No ghost - lw 13.5 overpass - 18-22 UTC',behr_coast));
                monthly_converg_path = fullfile(workdir,sprintf('%s US BEHR Monthly - Convergence',behr_coast));
                monthly_sqrt_converg_path = fullfile(workdir,sprintf('%s US BEHR Monthly - Sqrt Convergence',behr_coast));
                if any(strcmpi(apriori, {'hybrid','monthly-converg','monthly-sqrt-converg'}))
                    warning('The %s apriori has not been updated to use lonweight assuming OMI overpass at 1330 LST or profiles from 1800-2200 UTC', apriori);
                end
            case 'c'
                daily_path = fullfile(workdir,sprintf('%s US BEHR Hourly - No ghost - Coarse WRF',behr_coast));
                hybrid_path = fullfile(workdir,sprintf('%s US BEHR Hybrid - No ghost - Coarse WRF',behr_coast));
                monthly_path = fullfile(workdir,sprintf('%s US BEHR Monthly - No ghost - Coarse WRF - lw 13.5 overpass - 18-22 UTC',behr_coast));
                if any(strcmpi(apriori, {'hybrid','hourly'}))
                    warning('The coarse %s apriori has not been updated to use lonweight assuming OMI overpass at 1330 LST or profiles from 1800-2200 UTC', apriori);
                end
        end
    case 'pseudo-behr'
        if strcmpi(city,'SF')
            E.notimplemented('pseudo-behr for SF');
        end
        daily_file_name_spec = 'OMI_BEHR_%1$04d%2$02d%3$02d.mat';
        monthly_file_name_spec = 'OMI_BEHR_%1$04d%2$02d%3$02d.mat';
        xx = 1:11;
        yy = 1:19;
        s = 1;
        switch res
            case 'f'
                if strcmpi(timemode,'avg')
                    daily_path = fullfile(workdir, 'Atlanta BEHR Hourly - No clouds - No ghost - UTC 1800-2200');
                    hybrid_path = fullfile(workdir, 'Atlanta BEHR Hybrid - No clouds - No ghost - lonweight 13.5 overpass - UTC 1800-2200');
                    hybrid_avg_path = fullfile(workdir, 'Atlanta BEHR Avg Hybrid - No clouds - No ghost');
                    if any(strcmpi(apriori,{'hybrid-avg'}))
                        warning('The hour-averaged %s pseudo-retrieval apriori has not been updated to use lonweight assuming OMI overpass at 1330 LST or profiles from 1800-2200 UTC', apriori);
                    end
                else
                    daily_path = fullfile(workdir, 'Atlanta BEHR Hourly - No clouds - No ghost - Instantaneous');
                    hybrid_path = fullfile(workdir, 'Atlanta BEHR Hybrid - No clouds - No ghost - Instantaneous');
                    warning('The instantaneous pseudo-retrieval apriori have not been updated to use lonweight assuming OMI overpass at 1330 LST or profiles from 1800-2200 UTC');
                end
                monthly_path = fullfile(workdir, 'Atlanta BEHR Monthly - No clouds - No ghost - lonweight 13.5 overpass - UTC 1800-2200');
            case 'c'
                daily_path = fullfile(workdir, 'Atlanta BEHR Hourly - No clouds - No ghost - Coarse WRF');
                hybrid_path = fullfile(workdir, 'Atlanta BEHR Hybrid - No clouds - No ghost - Coarse WRF');
                hybrid_avg_path = fullfile(workdir, 'Atlanta BEHR Avg Hybrid - No clouds - No ghost - Coarse WRF');
                monthly_path = fullfile(workdir, 'Atlanta BEHR Monthly - No clouds - No ghost - Coarse WRF');
                warning('The coarse pseudo-retrieval apriori have not been updated to use lonweight assuming OMI overpass at 1330 LST or profiles from 1800-2200 UTC');
        end
end

daily_file_name = sprintf(daily_file_name_spec, year(date_in), month(date_in), day(date_in));
if strcmp(source,'wrf')
    monthly_file_name = sprintf(monthly_file_name_spec, year(date_in), month(date_in), eomday(year(date_in), month(date_in)));
else
    monthly_file_name = sprintf(monthly_file_name_spec, year(date_in), month(date_in), day(date_in));
end

switch apriori
    case 'monthly'
        file_out = fullfile(monthly_path, monthly_file_name);
    case 'hybrid'
        file_out = fullfile(hybrid_path, daily_file_name);
    case 'hybrid-avg'
        file_out = fullfile(hybrid_avg_path, daily_file_name);
    case 'hourly'
        file_out = fullfile(daily_path, daily_file_name);
    case 'monthly-lonwt13.5'
        file_out = fullfile(monthly_path_lonwt13, monthly_file_name);
    case 'monthly-converg'
        file_out = fullfile(monthly_converg_path, monthly_file_name);
    case 'monthly-sqrt-converg'
        file_out = fullfile(monthly_sqrt_converg_path, monthly_file_name);
end
end
