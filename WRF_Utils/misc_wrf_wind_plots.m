function [ varargout ] = misc_wrf_wind_plots( plttype, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;
wrf_path = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US';

if ~exist('plttype','var')
    allowed_plots = {'wrf_wind-cld','wrf-sat-cld','wrf-var-corr','3dvar-v-cld'};
    plttype = ask_multichoice('Select a plot type',allowed_plots,'list',true);
end

switch plttype
    case 'wrf_wind-cld'
        wrf_wind_cld_corr();
    case 'wrf-sat-cld'
        wrf_cloud_corr();
    case 'wrf-var-corr'
        wrf_var_corr();
    case '3dvar-v-cld'
        var3D_vs_cld(varargin{:});
    case 'match2behr'
        [varargout{1}, varargout{2}, varargout{3}] = wrf_match_to_behr(varargin{:});
end

    function wrf_wind_cld_corr
        sdate = datenum(ask_date('Enter the starting date'));
        edate = datenum(ask_date('Enter the ending date'));
        wrf_path = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US';
        
        city_lon = -84.39;
        city_lat = 33.775;
        
        winddir = [];
        cldfra = [];
        
        for d=sdate:edate
            fprintf('Loading files for %s\n',datestr(d));
            wrf_pat = sprintf('wrfout_d01_%s*',datestr(d,'yyyy-mm-dd'));
            F = dir(fullfile(wrf_path,wrf_pat));
            [this_U, this_V, this_cosa, this_sina, this_cldfra, this_xlon, this_xlat] = read_wrf_vars(wrf_path, F, {'U','V','COSALPHA','SINALPHA','CLDFRA','XLONG','XLAT'},false,0);
            city_inds = find_nearest_gridpoint(city_lon, city_lat, this_xlon(:,:,1), this_xlat(:,:,1));
            cx = city_inds(1);
            cy = city_inds(2);
            [Ue, Ve] = wrf_winds_transform(this_U, this_V, this_cosa, this_sina);
            Ue_bar = nanmean(reshape(Ue(cx-1:cx+1,cy-1:cy+1,1:3,1,:),[],1));
            Ve_bar = nanmean(reshape(Ve(cx-1:cx+1,cy-1:cy+1,1:3,1,:),[],1));
            
            theta = atan2d(Ve_bar, Ue_bar);
            winddir = cat(1, winddir, theta);
            
            max_overlap_cldfra = squeeze(max(this_cldfra,[],3));
            cldfra = cat(1,cldfra,nanmean(max_overlap_cldfra(:)));
        end
        
        figure;
        scatter(winddir, cldfra);
        xlabel('Winddir (CCW from E)')
        ylabel('Average cloud fraction');
    end

    function wrf_cloud_corr()
        start_date = ask_date('Enter the start date');
        end_date = ask_date('Enter the ending date');
        cld_calc = ask_multichoice('Use the random or maximum overlap cloud total?',{'random','max'});
        
        domain_x = [-87.1, -81.9];
        domain_y = [31.9, 35.6];
        
        % Load each BEHR files and WRF file in the date range. Find the
        % BEHR pixels within the domain of interest. Average the WRF
        % clouds to them and add those values to the respective vectors.
        omi_clds_vec = [];
        modis_clds_vec = [];
        wrf_clds_vec = [];
        
        behr_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hybrid - No ghost';
        behr_files = dir(fullfile(behr_path,'OMI*.mat'));
        behr_dates = nan(size(behr_files));
        for a=1:numel(behr_files)
            [s,e] = regexp(behr_files(a).name,'\d\d\d\d\d\d\d\d');
            behr_dates(a) = datenum(behr_files(a).name(s:e),'yyyymmdd');
        end
        wrf_path = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US_BEHR/hourly-new';
        wrf_files = dir(fullfile(wrf_path,'WRF*.nc'));
        wrf_dates = nan(size(wrf_files));
        for a=1:numel(wrf_dates)
            [s,e] = regexp(wrf_files(a).name,'\d\d\d\d-\d\d-\d\d');
            wrf_dates(a) = datenum(wrf_files(a).name(s:e), 'yyyy-mm-dd');
        end
        
        for d=datenum(start_date):datenum(end_date)
            fprintf('Working on %s\n',datestr(d));
            bb = behr_dates == d;
            if sum(bb) < 1
                fprintf('\tNo BEHR file for %s; skipping\n',datestr(d));
                continue
            elseif sum(bb) > 1
                E.toomanyfiles(sprintf('BEHR file for %s',datestr(d)));
            end
            ww = wrf_dates == d;
            if sum(ww) < 1
                fprintf('\tNo WRF files for %s; skipping\n',datestr(d));
            elseif sum(ww) > 1
                E.toomanyfiles(sprintf('WRF file for %s',datestr(d)));
            end
            
            % Read in the necessary WRF variables
            wi = ncinfo(fullfile(wrf_path,wrf_files(a).name));
            wrf_lon = double(ncread(wi.Filename, 'XLONG'));
            wrf_lon = wrf_lon(:,:,1); % lon doesn't change with time.
            wrf_lat = double(ncread(wi.Filename, 'XLAT'));
            wrf_lat = wrf_lat(:,:,1); % neither does lat.
            wrf_hr = double(ncread(wi.Filename, 'utchr'));
            if strcmpi(cld_calc,'random')
                wrf_cld = wrf_total_cldfra(wi.Filename,'r',0);
            elseif strcmpi(cld_calc,'max')
                tmp_cld = double(ncread(wi.Filename,'CLDFRA'));
                wrf_cld = squeeze(max(tmp_cld,[],3));
            end
            
            % Find the BEHR pixels for each swath that lie entirely in the
            % domain
            D = load(fullfile(behr_path, behr_files(bb).name),'Data');
            Data = D.Data;
            for s=1:numel(Data)
                minloncorn = squeeze(min(Data(s).Loncorn,[],1));
                maxloncorn = squeeze(max(Data(s).Loncorn,[],1));
                minlatcorn = squeeze(min(Data(s).Latcorn,[],1));
                maxlatcorn = squeeze(max(Data(s).Latcorn,[],1));
                pp = minloncorn >= domain_x(1) & maxloncorn <= domain_x(2) & minlatcorn >= domain_y(1) & maxlatcorn <= domain_y(2);
                if sum(pp) < 1;
                    continue
                end
                
                loncorn_in = Data(s).Loncorn(:,pp);
                latcorn_in = Data(s).Latcorn(:,pp);
                omi_clouds_in = Data(s).CloudFraction(pp);
                modis_clouds_in = Data(s).MODISCloud(pp);
                minloncorn_in = minloncorn(pp);
                maxloncorn_in = maxloncorn(pp);
                minlatcorn_in = minlatcorn(pp);
                maxlatcorn_in = maxlatcorn(pp);
                
                % Find the WRF grid cells that match each pixel using our
                % usual method of eliminating grid cells through logical
                % operations first. Also find the time from the
                % BEHRapriorimode field.
                apr_str = Data(s).BEHRaprioriMode;
                st = regexp(apr_str,'[')+1;
                ed = regexp(apr_str,']')-1;
                apr_hrs_cell = strsplit(apr_str(st:ed),';');
                apr_hrs = nan(size(apr_hrs_cell));
                for a=1:numel(apr_hrs)
                    apr_hrs(a) = str2double(apr_hrs_cell{a});
                end
                
                for p=1:numel(omi_clouds_in)
                    omi_clds_vec(end+1) = omi_clouds_in(p);
                    modis_clds_vec(end+1) = modis_clouds_in(p);
                    
                    tt = ismember(wrf_hr,apr_hrs);
                    wrf_cld_slice = nanmean(wrf_cld(:,:,tt),3);
                    xx = wrf_lon >= minloncorn_in(p) & wrf_lon <= maxloncorn_in(p) & wrf_lat >= minlatcorn_in(p) & wrf_lat <= maxlatcorn_in(p);
                    in = inpolygon(wrf_lon(xx), wrf_lat(xx), loncorn_in(:,p), latcorn_in(:,p));
                    wrf_cld_pix = wrf_cld_slice(xx);
                    wrf_clds_vec(end+1) = nanmean(wrf_cld_pix(in));
                end
            end
        end
        
        figure;
        scatter(omi_clds_vec, wrf_clds_vec);
        set(gca,'fontsize',20);
        xlabel('OMI clouds');
        ylabel('WRF clouds');
        
        figure;
        scatter(modis_clds_vec, wrf_clds_vec);
        set(gca,'fontsize',20);
        xlabel('MODIS clouds');
        ylabel('WRF clouds');
    end

    function wrf_var_corr()
        % See if there is any correlation between given WRF variables and
        % total cloudfraction.  Will average over days to account for
        % average cloudiness, since WRF output is instantaneous.
        %
        % A better way might be to use the TEMPO output to look at average
        % cloud fraction preceeding the measurement.
        yvar = ask_multichoice('Which WRF quantity to plot again total cldfrac?',{'no2:no','loss-nox'});
        
        wrf_path = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US';
        F = dir(fullfile(wrf_path,'wrfout*'));
        
        xdata = [];
        ydata = [];
        
        for a=1:numel(F)
            fprintf('Loading file %d of %d\n',a,numel(F));
            wi = ncinfo(fullfile(wrf_path,F(a).name));
            min_cldfra = ncread(wi.Filename, 'CLDFRA');
            min_cldfra = max(min_cldfra,[],3);
            xdata = cat(1, xdata, min_cldfra(:));
            
            switch yvar
                case 'no2:no'
                    no = double(ncread(wi.Filename,'no'));
                    no = no(:,:,1);
                    no2 = double(ncread(wi.Filename,'no2'));
                    no2 = no2(:,:,1);
                    ydata = cat(1, ydata, no2(:)./no(:));
                    ystr = 'NO2/NO';
                case 'loss-nox'
                    loss_nox = nanmean(double(ncread(wi.Filename,'LNOX')),4);
                    ydata = cat(1,ydata,loss_nox(:));
                    ystr = 'L(NO_x)';
            end
        end
        
        figure;
        scatter(xdata,ydata)
        set(gca,'fontsize',20);
        xlabel('Total cloudfrac')
        ylabel(ystr);
    end

    function var3D_vs_cld(omi_cld, quadrant_vec, wrf_var, P, PB)
        % Makes quadrant-wise plot of average vertical wind at each level
        % sorted by OMI cloud fraction. Can either pass the OMI cloud,
        % quadrant, W, P, and PB variables in directly or it will
        % automatically call wrf_match_to_behr to get these.
        
        
        plot_mode = ask_multichoice('Plot average, median w/quartiles, or avg. abs val?',{'avg','med','abs'});
        
        if nargin < 5
            wrf_var_name = input('Enter the WRF variable to plot: ','s');
            if strcmpi(wrf_var_name,'q')
                fprintf('Aborting plot\n');
                return
            end
            vars = [wrf_var_name, ' P PB']; % get the pressure as well (base and perturbation) as the vertical coordinate.
            
            opts.quad_bool = true;
            opts.dist_limit = 150;
            opts.size_lim_type = 'row';
            opts.lim_crit = 5;
            opts.center_lon = -84.39;
            opts.center_lat = 33.775;
            [behr_vars, wrf_vars, quadrant_vec] = wrf_match_to_behr(opts, '2013-06-01','2013-08-30',vars);
            omi_cld = behr_vars{1};
            wrf_var = wrf_vars{1};
            P = wrf_vars{2};
            PB = wrf_vars{3};
        else
            xstr = input('How should the x-axis be labeled? ','s');
            
            n(1) = numel(omi_cld);
            n(2) = numel(quadrant_vec);
            n(3) = size(wrf_var,2);
            n(4) = size(P,2);
            n(5) = size(PB,2);
            if any(n ~= n(1))
                E.badinput('omi_cld and quadrant_vec must have the same number of elements as the length of the 2nd dimension of W, P, and PB')
            end
            s(1) = size(wrf_var,1)-1;
            s(2) = size(P,1);
            s(3) = size(PB,1);
            if any(s ~= s(1))
                E.badinput('W is expected to have one more element in the first dimension than P and PB');
            end
        end
        
        % unstagger the variable in the vertical direction if needed and
        % combine P and PB to get box center pressure, converted from Pa to
        % hPa
        if size(wrf_var,1) == size(P,1)+1
            wrf_var = unstagger(wrf_var,1);
        end
        if any(size(wrf_var) ~= size(P))
            E.sizeMismatch('wrf_var','P')
        end
        box_P = (P+PB)/100;
        
        % make four plots (one for each quadrant), each with five subplots
        % (one for each cloud fraction bin).
        nbins = 5;
        cld_ll = (0:4)/nbins;
        cld_ul = (1:5)/nbins;
        quad_names = {'NE','NW','SW','SE'};
        
        % get the x label with the variable name and units if the
        % variable was read directly. Otherwise we got the label back in
        % the beginning.
        if exist('wrf_var_name','var')
            W = dir(fullfile(wrf_path,'wrfout*'));
            wi = ncinfo(fullfile(wrf_path,W(1).name));
            xx = strcmp(wrf_var_name, {wi.Variables.Name});
            uu = strcmpi('units',{wi.Variables(xx).Attributes.Name});
            units = wi.Variables(xx).Attributes(uu).Value;
            xstr = sprintf('%s (%s)',wrf_var_name,units);
        end
        
        
        for a=1:4
            figure;
            qq = quadrant_vec == a;
            for b=1:nbins
                cc = omi_cld >= cld_ll(b) & omi_cld <= cld_ul(b);
                xx = qq & cc;
                
                subplot(1,nbins,b);
                switch plot_mode
                    case 'avg'
                        this_w = nanmean(wrf_var(:,xx),2);
                        this_w_std = nanstd(wrf_var(:,xx),0,2);
                    case 'med'
                        this_w = nanmedian(wrf_var(:,xx),2);
                        this_w_quant = quantile(wrf_var(:,xx),[0.25 0.75],2);
                    case 'abs'
                        this_w = nanmean(abs(wrf_var(:,xx)),2);
                        this_w_std = nanstd(abs(wrf_var(:,xx)),0,2);
                end
                this_p = nanmean(box_P(:,xx),2);
                
                line(this_w, this_p, 'color','k','linewidth',2);
                if strcmpi(plot_mode,'med')
                    scatter_errorbars(this_w, this_p, this_w - this_w_quant(:,1), this_w_quant(:,2) - this_w, 'direction', 'x');
                else
                    scatter_errorbars(this_w, this_p, this_w_std, 'direction', 'x');
                end
                
                title(sprintf('Cld. frac. %.1f - %.1f',cld_ll(b),cld_ul(b)));
                
                
                xlabel(xstr);
                if b==1
                    ylabel('Pressure (hPa)')
                end
                ylim([200 1000]);
                set(gca,'ydir','reverse');
                set(gca,'xgrid','on');
            end
            if exist('wrf_var_name','var')
                suptitle(sprintf('%s - %s quadrant - %s',wrf_var_name, quad_names{a},plot_mode));
            else
                suptitle(sprintf('%s quadrant - %s',quad_names{a},plot_mode));
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% DATA FUNCTIONS %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [behr_vars, wrf_vars, quadrant_vec] = wrf_match_to_behr(options, start_date, end_date, wrf_var_names)
        % Returns vectors or matrices of BEHR and WRF data matched so that
        % every WRF data point has a corresponding BEHR data point.
        %
        % If no input given, will ask interactively for all options.
        % Give wrf_var_names as a space separated string of variable names.
        if ~exist('start_date','var')
            start_date = ask_date('Enter the starting date');
        end
        if ~exist('end_date','var')
            end_date = ask_date('Enter the ending date');
        end
        % could be turned on later if needed
        %behr_var_names = input('Enter the BEHR variables to save, separated by spaces: ', 's');
        %behr_var_names = strsplit(behr_var_names, ' ');
        behr_var_names = {'CloudFraction'};
        if ~exist('wrf_var_names','var')
            wrf_var_names = input('Enter the WRF variables to match to BEHR pixels, separated by spaces: ','s');
        end
        wrf_var_names = strsplit(wrf_var_names,' ');
        
        domain_lonlim = [-87.1 -81.9];
        domain_latlim = [31.9 35.6];
        
        behr_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hybrid - No ghost';
        wrf_path = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US';
        
        % Loop over all data. Load the BEHR file first to determine which
        % WRF file we need. The first time through it will trigger further
        % subsetting options unless these have been passed in as input.
        quadrant_vec = [];
        behr_vars = cell(size(behr_var_names));
        wrf_vars = cell(size(wrf_var_names));
        
        if isDisplay
            wb=waitbar(0, 'Concatenating data');
        end
        
        for d=datenum(start_date):datenum(end_date)
            if isDisplay
                wbfrac = (d-datenum(start_date))/(datenum(end_date)-datenum(start_date));
                waitbar(wbfrac);
            end
            B = load(fullfile(behr_path, sprintf('OMI_BEHR_%s.mat',datestr(d,'yyyymmdd'))),'Data');
            Data = B.Data;
            for s=1:numel(Data)
                if exist('options','var')
                    [in, quadrants] = subset_BEHR_pixels(Data(s), domain_lonlim, domain_latlim, options);
                else
                    [in, quadrants, options] = subset_BEHR_pixels(Data(s), domain_lonlim, domain_latlim);
                end
                
                if isempty(in)
                    continue
                end
                
                loncorn_in = Data(s).Loncorn(:,in);
                latcorn_in = Data(s).Latcorn(:,in);
                
                % Identify the WRF files we want from the UTC hour given in
                % the BEHRaprioriMode field
                st = regexp(Data(s).BEHRaprioriMode,'[');
                ed = regexp(Data(s).BEHRaprioriMode,']');
                utchrs = unique(eval(Data(s).BEHRaprioriMode(st:ed)));
                for h=1:numel(utchrs)
                    wrf_pattern = sprintf('wrfout_d01_%s_%02d*',datestr(d,'yyyy-mm-dd'),utchrs(h));
                    F = dir(fullfile(wrf_path,wrf_pattern));
                    for f=1:numel(F)
                        xlon = ncread(fullfile(wrf_path,F(f).name),'XLONG');
                        xlat = ncread(fullfile(wrf_path,F(f).name),'XLAT');
                        % Read in each WRF variable so we don't read them
                        % from the NETCDF files for every pixel
                        wrf_var_tmp = cell(1,numel(wrf_var_names));
                        for v=1:numel(wrf_var_names)
                            this_wrf_var = ncread(fullfile(wrf_path,F(f).name),wrf_var_names{v});
                            % unstagger if neccessary
                            if size(this_wrf_var,1) == size(xlon,1)+1
                                this_wrf_var = unstagger(this_wrf_var,1);
                            elseif size(this_wrf_var,2) == size(xlon,2)+1
                                this_wrf_var = unstagger(this_wrf_var,2);
                            end
                            wrf_var_tmp{v} = this_wrf_var;
                        end
                        
                        % Now loop over every remaining pixel. Find the WRF
                        % profiles that match it. Replicate the BEHR
                        % variables to match the number of WRF profiles and
                        % concatenate both as vectors or matrices depending
                        % if they are 2- or 3-D.
                        for p=1:size(loncorn_in,2)
                            xx = is_in_pixel(xlon,xlat,loncorn_in(:,p),latcorn_in(:,p));
                            quadrant_vec = cat(1, quadrant_vec, repmat(quadrants(p),sum(xx(:)),1));
                            
                            for v=1:numel(behr_var_names)
                                if ~ismatrix(Data(s).(behr_var_names{v}))
                                    % 3D BEHR variables should have the
                                    % "profile" dimension first
                                    this_var = Data(s).(behr_var_names{v})(:,p);
                                    behr_vars{v} = cat(2, behr_vars{v}, repmat(this_var, 1, sum(xx(:))));
                                else
                                    this_var = Data(s).(behr_var_names{v})(p);
                                    behr_vars{v} = cat(1, behr_vars{v}, repmat(this_var, sum(xx(:)), 1));
                                end
                            end
                            
                            for v=1:numel(wrf_var_names)
                                this_wrf_var = wrf_var_tmp{v};
                                nd = ndims(this_wrf_var);
                                if nd == 2
                                    this_wrf_var = this_wrf_var(xx);
                                    wrf_vars{v} = cat(1,wrf_vars{v},this_wrf_var);
                                elseif nd == 3
                                    this_wrf_var = permute(this_wrf_var, [3 1 2]); % put the vertical dimension first
                                    this_wrf_var = this_wrf_var(:,xx);
                                    wrf_vars{v} = cat(2,wrf_vars{v},this_wrf_var);
                                else
                                    E.notimplemented(sprintf('nd = %d',nd));
                                end
                            end
                        end
                    end
                end
            end
        end
        if isDisplay
            close(wb);
        end
    end
end

