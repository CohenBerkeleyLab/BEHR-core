function [ varargout ] = misc_wrf_wind_plots( plttype, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;
wrf_path = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US';

if ~exist('plttype','var')
    allowed_plots = {'wrf_wind-cld','wrf-sat-cld','wrf-var-corr','2dvar-v-cld','3dvar-v-cld','nonlin-nox','wrf_stats','pick-rxns','mean-react','daily-react','sum-emis'};
    plttype = ask_multichoice('Select a plot type',allowed_plots,'list',true);
end

switch plttype
    case 'wrf_wind-cld'
        wrf_wind_cld_corr();
    case 'wrf-sat-cld'
        wrf_cloud_corr(varargin{:});
    case 'wrf-var-corr'
        wrf_var_corr();
    case '2dvar-v-cld'
        var2D_vs_cld();
    case '3dvar-v-cld'
        var3D_vs_cld(varargin{:});
    case 'nonlin-nox'
        nonlin_nox();
    case 'match2behr'
        [varargout{1}, varargout{2}, varargout{3}] = wrf_match_to_behr(varargin{:});
    case 'wrf_stats'
        varargout{1} = wrf_cell_stats(varargin{:});
    case 'pick-rxns'
        pick_wrf_rxns();
    case 'mean-react'
        [varargout{1}, varargout{2}] = calc_mean_total_react();
    case 'daily-react'
        [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = calc_daily_total_reactivity();
    case 'sum-emis'
        varargout{1} = compute_wrf_emis();
    case 'lifetime'
        compute_avg_lifetime;
end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% PLOTTING FUNCTION %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

    function wrf_cloud_corr(opts)
        start_date = ask_date('Enter the start date');
        end_date = ask_date('Enter the ending date');
        cld_calc = ask_multichoice('Use the random or maximum overlap cloud total?',{'random','max'});
        
        domain_x = [-87.1, -81.9];
        domain_y = [31.9, 35.6];
        
        quad_avg_bool = strcmpi(ask_multichoice('Do you want to compare daily quadrant averages rather than pixel-by-pixel?',{'y','n'}),'y');
        if quad_avg_bool && ~exist('opts','var')
            opts = subset_BEHR_pixels;
            fprintf('*** Center lat/lon will be set to 33.775, -84.39 ***\n')
        end
        % Overwrite the center lat/lon since this is only set up for
        % Atlanta.
        opts.center_lon = -84.39;
        opts.center_lat = 33.775;
        
        % Load each BEHR files and WRF file in the date range. Find the
        % BEHR pixels within the domain of interest. Average the WRF
        % clouds to them and add those values to the respective vectors.
        
        behr_path = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/SE US BEHR Hybrid - No ghost';
        behr_files = dir(fullfile(behr_path,'OMI*.mat'));
        behr_dates = nan(size(behr_files));
        for a=1:numel(behr_files)
            [s,e] = regexp(behr_files(a).name,'\d\d\d\d\d\d\d\d');
            behr_dates(a) = datenum(behr_files(a).name(s:e),'yyyymmdd');
        end
        wrf_path = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US_BEHR/hourly';
        wrf_files = dir(fullfile(wrf_path,'WRF*.nc'));
        wrf_dates = nan(size(wrf_files));
        for a=1:numel(wrf_dates)
            [s,e] = regexp(wrf_files(a).name,'\d\d\d\d-\d\d-\d\d');
            wrf_dates(a) = datenum(wrf_files(a).name(s:e), 'yyyy-mm-dd');
        end
        
        datevec = datenum(start_date):datenum(end_date);
        
        if ~quad_avg_bool
            omi_clds_vec = [];
            modis_clds_vec = [];
            wrf_clds_vec = [];
        else
            omi_clds_vec = nan(4,numel(datevec));
            modis_clds_vec = nan(4,numel(datevec));
            wrf_clds_vec = nan(4,numel(datevec));
        end
        nd=0;
        for d=datevec;
            nd=nd+1; % lazy and not refactoring code to loop over number of elements of datevec
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
            if ~quad_avg_bool
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
                    apr_hrs = get_apri_hr(Data(s));
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
            else
                for s=1:numel(Data)
                    apr_hr = unique(get_apri_hr(Data(s)));
                    if ~isscalar(apr_hr)
                        dum=1;
                    end
                    Data(s).Apriori_hour = repmat(apr_hr(:),[1,size(Data(s).Longitude)]);
                end
                AllData = struct('Longitude',[],'Latitude',[],'Loncorn',[],'Latcorn',[],'CloudFraction',[],'MODISCloud',[],'Row',[],'Apriori_hour',[]);
                fns = fieldnames(AllData);
                for f=1:numel(fns)
                    % Loncorn and latcorn have corners along the first
                    % dimension so we need to concatenated them along the
                    % second
                    AllData.(fns{f}) = padcat(ndims(Data(1).(fns{f}))-1, Data.(fns{f}));
                end
                [in, behr_quadrants] = subset_BEHR_pixels(AllData, domain_x, domain_y, opts);
                omi_subset = AllData.CloudFraction(in);
                modis_subset = AllData.MODISCloud(in);
                apriori_hours = unique(AllData.Apriori_hour(:,in));
                apriori_hours(isnan(apriori_hours)) = [];
                
                [in, wrf_quadrants] = subset_WRF_grid_cells(wrf_lon, wrf_lat, domain_x, domain_y, opts);
                tt = ismember(wrf_hr,apriori_hours);
                wrf_cld_slice = nanmean(wrf_cld(:,:,tt),3);
                wrf_subset = wrf_cld_slice(in);
                
                for q=1:4
                    qqb = behr_quadrants == q;
                    qqw = wrf_quadrants == q;
                    omi_clds_vec(q,nd) = nanmean(omi_subset(qqb));
                    modis_clds_vec(q,nd) = nanmean(modis_subset(qqb));
                    wrf_clds_vec(q,nd) = nanmean(wrf_subset(qqw));
                end
            end
        end
        
        if quad_avg_bool
            nq = 4;
            titlestr = 'Daily average in %s quadrant';
        else
            nq = 1;
            titlestr = 'WRF matched to BEHR pixels';
        end
        quad_names = {'NE','NW','SW','SE'};
        
        for q=1:nq
            figure;
            scatter(omi_clds_vec(q,:), wrf_clds_vec(q,:));
            set(gca,'fontsize',20);
            xlabel('OMI clouds');
            ylabel('WRF clouds');
            title(sprintf(titlestr, quad_names{q}));
            
            figure;
            scatter(modis_clds_vec(q,:), wrf_clds_vec(q,:));
            set(gca,'fontsize',20);
            xlabel('MODIS clouds');
            ylabel('WRF clouds');
            title(sprintf(titlestr, quad_names{q}));
        end
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

    function var2D_vs_cld(omi_cld, quadrant_vec, wrf_var)
        % Makes quadrant-wise scatter plot of any 2D in WRF against OMI
        % cloud fraction
        if nargin < 3
            wrf_var_names = input('Enter the WRF variables to plot, separated by spaces: ', 's');
            if strcmpi(wrf_var_names,'q')
                fprintf('Aborting plot\n');
                return
            end
            opts.quad_bool = true;
            opts.dist_limit = 150;
            opts.size_lim_type = 'row';
            opts.lim_crit = 5;
            opts.center_lon = -84.39;
            opts.center_lat = 33.775;
            [behr_vars, wrf_vars, quadrant_vec] = wrf_match_to_behr(opts, '2013-06-01','2013-08-30','CloudFraction',wrf_var_names);
            omi_cld = behr_vars{1};
        else
            ystr = input('How should the x-axis be labeled? ','s');
            
            n(1) = numel(omi_cld);
            n(2) = numel(quadrant_vec);
            for a=1:numel(wrf_var)
                n(a+2) = numel(wrf_var{a});
            end
            if any(n ~= n(1))
                E.badinput('omi_cld, quadrant_vec, and all WRF vars must have the same number of elements')
            end
        end
        
        % get the x label with the variable name and units if the
        % variable was read directly. Otherwise we got the label back in
        % the beginning, but convert to a cell array
        wrf_var_names = strsplit(wrf_var_names,' ');
        if exist('wrf_var_names','var')
            ystr = cell(numel(wrf_var_names));
            W = dir(fullfile(wrf_path,'wrfout*'));
            wi = ncinfo(fullfile(wrf_path,W(1).name));
            for a=1:numel(wrf_var_names)
                xx = strcmp(wrf_var_names{a}, {wi.Variables.Name});
                uu = strcmpi('units',{wi.Variables(xx).Attributes.Name});
                units = wi.Variables(xx).Attributes(uu).Value;
                ystr{a} = sprintf('%s (%s)',wrf_var_names{a},units);
            end
        else
            ystr = {ystr};
        end
        
        % make four plots (one for each quadrant) times the number of WRF
        % variables.
        quad_names = {'NE','NW','SW','SE'};
        for a=1:numel(wrf_vars)
            for q=1:4
                figure;
                qq = quadrant_vec == q;
                scatter(omi_cld(qq), wrf_vars{a}(qq))
                xlabel('OMI cloud fraction')
                ylabel(ystr{a})
                title(quad_names{a})
            end
        end
    end

    function var3D_vs_cld(omi_cld, quadrant_vec, wrf_var, P, PB)
        % Makes quadrant-wise plot of average vertical wind at each level
        % sorted by OMI cloud fraction. Can either pass the OMI cloud,
        % quadrant, W, P, and PB variables in directly or it will
        % automatically call wrf_match_to_behr to get these.
        
        prof_bool = strcmpi(ask_multichoice('Plot as profiles or scatter?',{'prof','scatter'}),'prof');
        if prof_bool
            plot_mode = ask_multichoice('Plot average, median w/quartiles, or avg. abs val?',{'avg','med','abs'});
            layer_mode = ask_multichoice('Average over layer or interp to std. pres vec first?',{'layer','pres'});
            interp_to_std_p = strcmpi(layer_mode,'pres');
        else
            plot_mode = 'scatter';
            interp_to_std_p = false;
        end
            
        
        if nargin < 5
            wrf_var_name = input('Enter the WRF variable to plot (windmag will get horizontal wind speed): ','s');
            if strcmpi(wrf_var_name,'q')
                fprintf('Aborting plot\n');
                return
            elseif strcmpi(wrf_var_name,'windmag')
                vars = 'U V P PB';
            else
                vars = [wrf_var_name, ' P PB']; % get the pressure as well (base and perturbation) as the vertical coordinate.
            end
            
            opts.quad_bool = true;
            opts.dist_limit = 150;
            opts.size_lim_type = 'row';
            opts.lim_crit = 5;
            opts.center_lon = -84.39;
            opts.center_lat = 33.775;
            [behr_vars, wrf_vars, quadrant_vec] = wrf_match_to_behr(opts, '2013-06-01','2013-08-30','CloudFraction',vars);
            omi_cld = behr_vars{1};
            if strcmpi(wrf_var_name,'windmag')
                wrf_var = sqrt(wrf_vars{1}.^2 + wrf_vars{2}.^2);
                P = wrf_vars{3};
                PB = wrf_vars{4};
            else
                wrf_var = wrf_vars{1};
                P = wrf_vars{2};
                PB = wrf_vars{3};
            end
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
        
        % get the x label with the variable name and units if the
        % variable was read directly. Otherwise we got the label back in
        % the beginning.
        if exist('wrf_var_name','var')
            if ~strcmpi(wrf_var_name,'windmag')
                W = dir(fullfile(wrf_path,'wrfout*'));
                wi = ncinfo(fullfile(wrf_path,W(1).name));
                xx = strcmp(wrf_var_name, {wi.Variables.Name});
                uu = strcmpi('units',{wi.Variables(xx).Attributes.Name});
                units = wi.Variables(xx).Attributes(uu).Value;
                xstr = sprintf('%s (%s)',wrf_var_name,units);
            else
                xstr = 'UV wind speed (m/s)';
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
        
        if interp_to_std_p
            [wrf_var, box_P] = interp_to_std_p(wrf_var, box_P);
        end
        
        % make four plots (one for each quadrant), each with five subplots
        % (one for each cloud fraction bin).
        nbins = 5;
        cld_ll = (0:4)/nbins;
        cld_ul = (1:5)/nbins;
        quad_names = {'NE','NW','SW','SE'};
        
        
        for a=1:4
            figure;
            qq = quadrant_vec == a;
            if prof_bool
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
            else
                omi_cld_plot = permute(omi_cld,[3 1 2]);
                omi_cld_plot = repmat(omi_cld_plot, size(wrf_var,1), 1, 1);
                omi_cld_plot = reshape(omi_cld_plot(:,qq),[],1);
                wrf_var_plot = reshape(wrf_var(:,qq),[],1);
                p_plot = reshape(box_P(:,qq),[],1);
                scatter3(omi_cld_plot, wrf_var_plot, p_plot)
                xlabel('Cloud fraction');
                ylabel(xstr);
                zlabel('Pressure (hPa)');
                set(gca,'zdir','reverse')
            end
            if exist('wrf_var_name','var')
                suptitle(sprintf('%s - %s quadrant - %s',wrf_var_name, quad_names{a},plot_mode));
            else
                suptitle(sprintf('%s quadrant - %s',quad_names{a},plot_mode));
            end
        end
    end

    function nonlin_nox
        % Will plot [OH] vs. [NOx] for all three cities in the SE region
        % I'm considering for the wind effects paper. This will hopefully
        % tell me something about whether Montgomery and Birmingham are
        % more sensitive to NOx concentrations in WRF than Atlanta is.
        radius = ask_number('Radius in WRF grid cells to consider (1 will use a 3x3 subset)', 'default', 0, 'testfxn', @(x) x>=0);
        model_layers = ask_number('How many model layers to consider?', 'default', 1, 'testfxn', @(x) x>=1);
        unit = ask_multichoice('Plot in mixing ratio or number density?',{'mr','nd'},'default','mr');
        
        city_names = {'Atlanta','Birmingham','Montgomery','Background'};
        city_lons = [   -84.39,...  %Atlanta
                        -86.80,...  %Birmingham
                        -86.30,...  %Montgomery
                        -83.00];    %Background
        city_lats = [   33.775,...  %Atlanta
                        33.52,...   %Birmingham
                        32.37,...   %Montgomery
                        32.25];     %Background
        city_xx = cell(size(city_lons));
        city_yy = cell(size(city_lons));
        
        F = dir(fullfile(wrf_path,'wrfout*'));
        first_time = true;
        
        city_oh = cell(size(city_lons));
        city_nox = cell(size(city_lons));
        for c=1:numel(city_lons)
            city_oh{c} = nan(size(F));
            city_nox{c} = nan(size(F));
        end
        
        if isDisplay
            wb = waitbar(0, 'Loading files');
        end
        
        for a=1:numel(F)
            if isDisplay
                waitbar(a/numel(F));
            end
            wrfname = fullfile(wrf_path,F(a).name);
            % restrict to after 1 June (treat end of May as spinup)
            [s,e] = regexp(F(a).name,'\d\d\d\d-\d\d-\d\d');
            if datenum(F(a).name(s:e),'yyyy-mm-dd') < datenum('2013-06-01')
                continue
            end
            
            if first_time
                % Assume that all the files have the same lat/lon grid as
                % each other so that we can just find the cities once.
                first_time = false;
                xlon = ncread(wrfname,'XLONG');
                xlat = ncread(wrfname,'XLAT');
                for c=1:numel(city_lons)
                    inds = find_nearest_gridpoint(city_lons(c), city_lats(c), xlon, xlat);
                    city_xx{c} = inds(1)-radius:inds(1)+radius;
                    city_yy{c} = inds(2)-radius:inds(2)+radius;
                end
            end
            
            OH = ncread(wrfname,'ho')*1e3; % ppm -> ppb
            NO = ncread(wrfname,'no')*1e3; % ppm -> ppb
            NO2 = ncread(wrfname,'no2')*1e3; %ppm -> ppb
            % Calculate [OH] and [NOx] in number density.
            if strcmpi(unit, 'nd')
                ndens_air = calculate_wrf_air_ndens(wrfname);
                OH = OH .* 1e-9 .* ndens_air; % ppb -> number density
                NO = NO .* 1e-9 .* ndens_air; % ppb -> number density
                NO2 = NO2 .* 1e-9 .* ndens_air; %ppb -> number density
            end
            % Subset for each city
            for c=1:numel(city_lons)
                oh_subset = OH(city_xx{c}, city_yy{c}, 1:model_layers);
                city_oh{c}(a) = nanmean(oh_subset(:));
                
                no_subset = NO(city_xx{c}, city_yy{c}, 1:model_layers);
                no2_subset = NO2(city_xx{c}, city_yy{c}, 1:model_layers);
                city_nox{c}(a) = nanmean(no_subset(:) + no2_subset(:));
            end
            
        end
        
        if isDisplay
            close(wb)
        end
        
        % Make the plots
        for c=1:numel(city_lons)
            figure;
            scatter(city_nox{c}, city_oh{c}, 32, 'k');
            title(sprintf('%1$s (%2$d x %2$d x %3$d)',city_names{c},radius*2+1,model_layers));
            xlabel('[NO_x] (molec. cm^{-3})');
            ylabel('[OH] (molec. cm^{-3})');
            set(gca,'fontsize',16);
        end
    end

    function compute_avg_lifetime
        % Will calculate the average lifetimes across WRF grid cells vs.
        % loss to HNO3, RONO2, or total.
        life_type = ask_multichoice('Which lifetime to calculate?',{'HNO3','RONO2','total'},'list',true);
        
        % Find all WRF files for the month of June at 1900 UTC
        F = dir(fullfile(wrf_path,'wrfout_d01_2013-06*_19-*'));
        
        wrf_varnames = {'ho'}; % any calculation needs OH, either for OH + NO2 -> HNO3 or RO2 steady-state calculation
        coeffs = [];
        if ismember(lower(life_type),{'hno3','total'});
            wrf_varnames = [wrf_varnames, {'ho'}];
            coeffs = [coeffs, 1];
            k_hno3 = @(TEMP, C_M) wrf_rate_expr('TROE', 1.49e-30 , 1.8 , 2.58e-11 , 0.0 , TEMP, C_M) ;
        end
        if ismember(lower(life_type),{'rono2','total'})
            [prods,prod_coeff,reacts,~,ro2_rates] = pick_wrf_rxns('LNOX-ANs',false);
            ro2_ss = cell(size(reacts));
            ro2_alpha = zeros(size(reacts));
            for a=1:numel(reacts)
                xx = ~strcmpi('NO',reacts{a});
                if numel(xx)>2
                    E.notimplemented('3 body reaction')
                elseif sum(~xx) ~= 1
                    continue
                else
                    ro2_ss{a} = compute_ss_ro2(reacts{a}{xx}, F, -1);
                    xx_no2 = strcmpi(prods{a},'NO2');
                    if sum(xx_no2) > 0
                        ro2_alpha(a) = 1 - prod_coeff{a}(xx_no2);
                    end
                end
            end
            xx = iscellcontents(ro2_ss,'isempty') | iscellcontents(ro2_ss, @(x) isscalar(x) && x==0);
            ro2_ss = ro2_ss(~xx);
            ro2_rates = ro2_rates(~xx);
        end
        
        
        
        % We'll also need temperature and pressure to calculate number
        % density of air and thus rate constants. Plus lat/lon for plotting
        wrf_varnames = [wrf_varnames, {'XLONG','XLAT','T','P','PB'}];
        
        % Load all WRF files for the month of June at 1900 UTC (or summer
        % if I ever download those). Get the date of each file to match to
        % wind speed later
        
        F_dnums = nan(size(F));
        for a=1:numel(F)
            [s,e] = regexp(F(a).name,'\d\d\d\d-\d\d-\d\d');
            F_dnums(a) = datenum(F(a).name(s:e));
        end
        wrfvars = cell(size(wrf_varnames));
        [wrfvars{:}] = read_wrf_vars(wrf_path, F, wrf_varnames);
        wrflon = wrfvars{end-4};
        wrflat = wrfvars{end-3};
        wrfT = wrfvars{end-2};
        wrfP = wrfvars{end-1};
        wrfPB = wrfvars{end};
        wrfvars(end-4:end) = [];
        
        wrf_temp = convert_wrf_temperature(wrfT, wrfP, wrfPB);
        wrf_ndens = calculate_wrf_air_ndens(wrfT, wrfP, wrfPB);
        
        opts.quad_bool = 0;
        opts.dist_limit = 50;%25;
        
        
        bins_lb = [-Inf, 3, 4, 5]; %m/s
        
        if ismember(lower(life_type),{'hno3','total'})
            k_hno3 = k_hno3(wrf_temp, wrf_ndens);
            tau_hno3 = 1 ./ (k_hno3 .* wrfvars{1} .* 1e-6 .* wrf_ndens); % concentration given in ppmv, convert to number density
            tau_hno3 = squeeze(nanmean(tau_hno3(:,:,1:5,:),3)) ./ 3600; % convert from seconds to hours
        end
        if ismember(lower(life_type),{'rono2','total'})
            tau_rono2 = zeros(size(ro2_ss{1}));
            % Lifetimes add in reciprocal (1/tau = 1/tau_1 + 1/tau_2 + ...)
            % but are themselves the reciprocal of alpha*k*[RO2], so add up
            % all the alpha*k*[RO2] first, then invert.
            for a=1:numel(ro2_ss)
                k_ro2 = ro2_rates{a}(wrf_temp, wrf_ndens);
                tau_rono2 = tau_rono2 + (ro2_alpha(a) .* ro2_ss{a} .* k_ro2);
            end
            tau_rono2 = 1 ./ tau_rono2;
            tau_rono2 = squeeze(nanmean(tau_rono2(:,:,1:5,:),3)) ./ 3600;
        end
        
        switch lower(life_type)
            case 'hno3'
                tau = tau_hno3;
                cb_label1 = '\tau vs HNO3';
                cb_label2 = '\sigma(\tau_{HNO3})';
            case 'rono2'
                tau = tau_rono2;
                cb_label1 = '\tau vs RONO2';
                cb_label2 = '\sigma(\tau_{RONO2})';
            case 'total'
                tau = 1 ./ (1 ./ tau_hno3 + 1./ tau_rono2);
                cb_label1 = '\tau vs HNO3 & RONO2';
                cb_label2 = '\sigma(\tau_{HNO3+RONO2})';
        end
        taubar = nanmean(tau, ndims(tau));
        tausd = nanstd(tau, 0, ndims(tau));
        
        figure;
        pcolor(wrflon(:,:,1), wrflat(:,:,1), taubar);
        cb=colorbar; cb.Label.String = cb_label1;
        city_names = {'Atlanta','Birmingham','Montgomery'};
        for a=1:numel(city_names)
            [city_lon, city_lat] = return_city_info(city_names{a});
            line(city_lon, city_lat, 'linestyle','none','marker','x','linewidth',2,'color','r')
        end
        figure;
        pcolor(wrflon(:,:,1), wrflat(:,:,1), tausd);
        cb=colorbar; cb.Label.String = cb_label2;
        city_names = {'Atlanta','Birmingham','Montgomery'};
        for a=1:numel(city_names)
            [city_lon, city_lat, city_windvel, ~, city_dnums] = return_city_info(city_names{a});
            line(city_lon, city_lat, 'linestyle','none','marker','x','linewidth',2,'color','r')
            opts.center_lon = city_lon;
            opts.center_lat = city_lat;
            in = subset_WRF_grid_cells(wrflon(:,:,1),wrflat(:,:,1),[-90,-80],[30 40],opts);
            sz = size(tau);
            tau_vec = reshape(tau,[sz(1)*sz(2), sz(3:end)]);
            fprintf('Average tau around %s:\n',city_names{a});
            
            for b=1:numel(bins_lb)
                xx = city_windvel >= bins_lb(b);
                tt = ismember(F_dnums, city_dnums(xx));
                city_taubar = nanmean(reshape(tau_vec(in,tt),[],1));
                city_tausd = nanstd(reshape(tau_vec(in,tt),[],1));
                fprintf('\t Wind >= %.2f m/s: %.2f +/- %.2f\n',bins_lb(b),city_taubar,city_tausd);
            end
        end
                
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% DATA FUNCTIONS %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [behr_vars, wrf_vars, quadrant_vec] = wrf_match_to_behr(options, start_date, end_date, behr_var_names, wrf_var_names)
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
        if ~exist('behr_var_names','var')
            behr_var_names = input('Enter the BEHR variables to save, separated by spaces: ', 's');
        end
        behr_var_names = strsplit(behr_var_names, ' ');
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

    function wrfstats = wrf_cell_stats(options)
        % This function will calculate statistics on all requested WRF
        % variables for the desired date range and UTC time.
        if ~exist('options','var') || ~isfield(options,'subset_lon') || ~isfield(options,'subset_lat')
            fprintf('Pixels can be subset by the nearest to a lat/lon or within a domain.\n')
            fprintf('For the next two questions, entering one number will find the grid cell nearest that point\n')
            fprintf('and two numbers will treat that as a domain.\n')
            sub_lon = ask_number('Enter the longitude criteria','testfxn',@(x) numel(x)<=2 && all(x >= -180) && all(x <= 180),'testmsg','Must be one or two values between -180 and +180');
            sub_lat = ask_number('Enter the longitude criteria','testfxn',@(x) numel(x)<=2 && all(x >= -180) && all(x <= 180),'testmsg','Must be one or two values between -180 and +180');
        else
            sub_lon = options.subset_lon;
            sub_lat = options.subset_lat;
        end
        
        if ~exist('options','var') || ~isfield(options,'model_layers')
            model_layers = ask_number('What model layers to consider?', 'testfxn', @(x) numel(x)==2 && x(1) <= x(2) && all(mod(x,1) == 0), 'testmsg', 'You must enter a 2-element vector of integers with the smaller value first');
        else
            model_layers = options.model_layers;
        end
        model_layers = min(model_layers):max(model_layers);
        
        if ~exist('options','var') || ~isfield(options,'wrf_vars');
            wrf_vars_names = strsplit(input('Enter the WRF variables to compute statistics for, separated by spaces: ', 's'));
        else
            if iscellstr(options.wrf_vars)
                wrf_vars_names = options.wrf_vars;
            else
                wrf_vars_names = strsplit(options.wrf_vars);
            end
        end
        
        if ~exist('options','var') || ~isfield(options,'utc_hr')
            utc_hr = ask_number('Enter the UTC hour to consider','testfxn',@(x) isscalar(x) && x>=0 && x<=23);
        else
            utc_hr = options.utc_hr;
        end
        
        if ~exist('options','var') || ~isfield(options, 'start_date') || ~isfield(options, 'start_date')
            start_date = ask_date('Enter the starting date');
            end_date = ask_date('Enter the ending_date');
        else
            start_date = options.start_date;
            end_date = options.end_date;
        end
        
        %%%% GET FILES AND READ %%%%%
        wrf_pattern = sprintf('wrfout_*_%02d-00-00',utc_hr);
        F = dir(fullfile(wrf_path,wrf_pattern));
        F = files_in_dates(F, start_date, end_date);
        
        wrf_vars = cell(size(wrf_vars_names));
        [wrf_vars{:}] = read_wrf_vars(wrf_path,F,wrf_vars_names);
        [xlon, xlat] = read_wrf_vars(wrf_path,F(1),{'XLONG','XLAT'});
        
        % Subset by finding one pixel or all grid cells in the domain
        % (later if needed)
        if isscalar(sub_lon) && isscalar(sub_lat)
            [~,mi] = min(sqrt((xlon(:) - sub_lon).^2) + sqrt((xlat(:) - sub_lat).^2));
            [mix, miy] = ind2sub(size(xlon),mi);
        elseif xor(isscalar(sub_lon), isscalar(sub_lat))
            E.badinput('sub_lon and sub_lat must be both scalar or not')
        else 
            E.notimplemented('subsetting by domain')
        end
        
        % Prep the output structure.
        substruct = struct('mean',0,'stddev',0,'median',0,'quartiles',[0 0],'percentile_1090',[0 0],'minmax',[0 0],'units','');
        wrfstats = make_empty_struct_from_cell(wrf_vars_names, substruct);
        
        % Go through each variable; slice as necessary, then compute the
        % stats. Also read the unit
        for a=1:numel(wrf_vars)
            if size(wrf_vars{a},3) > 1
                this_slice = reshape(wrf_vars{a}(mix,miy,model_layers,:),[],1);
            else
                this_slice = reshape(wrf_vars{a}(mix,miy,1,:),[],1);
            end
            wrfstats.(wrf_vars_names{a}).mean = nanmean(this_slice);
            wrfstats.(wrf_vars_names{a}).stddev = nanstd(this_slice);
            wrfstats.(wrf_vars_names{a}).median = nanmedian(this_slice);
            wrfstats.(wrf_vars_names{a}).quartiles = quantile(this_slice,[0.25 0.75]);
            wrfstats.(wrf_vars_names{a}).percentile_1090 = quantile(this_slice, [0.1 0.9]);
            wrfstats.(wrf_vars_names{a}).minmax = [min(this_slice), max(this_slice)];
            wrfstats.(wrf_vars_names{a}).units = ncreadatt(fullfile(wrf_path,(F(1).name)),wrf_vars_names{a},'units');
        end
    end

    function varargout = pick_wrf_rxns(subset_type, eval_k)
        % This will choose a subset of the reactions defined in the
        % r2smh.eqn file. If you need to use it automatically, pass the
        % subset desired as the first argument. If you wish to use the 'One
        % species' subset, pass it as a cell array with 'One species' as
        % the first element, the species name as the second and the side
        % ('products','reactants', or 'either') as the third.
        %
        % eval_k is optional; defaults to true meaning that the rate
        % constants will be evaluated at 298 K and 2x10^19 molec/cm^3
        % number density of air.
        %
        allowed_subsets = {'NO+RO2','VOC+OH','LNOX-ANs','One species'};
        allowed_sides = {'Products','Reactants','Either'};
        if ~exist('subset_type','var')
            subset_type = ask_multichoice('Which subset of rxns would you like?',allowed_subsets,'list',true);
            if strcmpi(subset_type,'One species')
                spec_name = input('    Enter the species name to look for: ', 's');
                spec_side = ask_multichoice('    Which side of the reaction should it be on?', allowed_sides);
            end
        else
            if iscell(subset_type)
                subset_opts = subset_type;
                subset_type = subset_opts{1};
            end
            if ~ismember(lower(subset_type), lower(allowed_subsets));
                E.badinput('If given, SUBSET_TYPE must be one of : %s',strjoin(allowed_subsets,', '));
            end
            if strcmpi(subset_type,'One species')
                spec_name = subset_opts{2};
                spec_side = subset_opts{3};
                if ~ismember(lower(spec_side), lower(allowed_sides))
                    E.badinput('If given SUBSET_TYPE{3} must be one of : %s',strjoin(allowed_sides,', ')); 
                end
            end
        end
        % If eval_k does not exist, assume that we want to go ahead and
        % evalute the rate constants for surface conditions.
        if ~exist('eval_k','var')
            eval_k = true;
        elseif ~isscalar(eval_k) || ~islogical(eval_k)
            E.badinput('eval_k must be a scalar logical')
        end
        
        % If output requested, we will be saving the products, reactants,
        % and rate constant functions for each reaction.
        if nargout > 0
            save_bool = true;
            prod_out = {};
            prod_c_out = {};
            react_out = {};
            react_c_out = {};
            k_out = {};
        else
            save_bool = 0;
        end
        
        homedir = getenv('HOME');
        eqn_file = fullfile(homedir,'Documents','MATLAB','WRF_Parsing','R2SMH','r2smh.eqn');
        fid = fopen(eqn_file);
        tline = fgetl(fid);
        while ischar(tline)
            if ismember('#',tline) || ~isempty(regexp(tline,'JUNK','once'))
            else
                [products, prod_coeff, reactants, react_coeff, rate_fxn] = read_wrf_mech_line(tline);
                % Get reactions of organics with NO to produce NO2. This
                % will let us calculate organic reactivity down the road.
                switch lower(subset_type)
                    case 'no+ro2'
                        eqn_bool = ismember('NO',reactants) && ismember('NO2',products) && ~any(ismember({'O3','O3P','HO2','NO3','M'},reactants));
                    case 'voc+oh'
                        eqn_bool = ismember('HO',reactants) && ~any(ismember({'O3','M','HO2','H2O2','NO','HONO','NO2','HNO3','NO3','HNO4','SO2'},reactants));
                    case 'lnox-ans'
                        eqn_bool = ismember('LNOXA',products);
                    case 'one species'
                        switch lower(spec_side)
                            case 'products'
                                eqn_bool = ismember(spec_name, products);
                            case 'reactants'
                                eqn_bool = ismember(spec_name, reactants);
                            case 'either'
                                eqn_bool = ismember(spec_name, products) || ismember(spec_name, reactants);
                        end
                end
                if eqn_bool
                    if save_bool
                        prod_out{end+1} = products;
                        prod_c_out{end+1} = prod_coeff;
                        react_out{end+1} = reactants;
                        react_c_out{end+1} = react_coeff;
                        if eval_k
                            k_out{end+1} = rate_fxn(298, 2e19); % evaluate for typical surface conditions.
                        else
                            k_out{end+1} = rate_fxn;
                        end
                    else
                        fprintf('%s\n',tline);
                    end
                end
            end
            tline = fgetl(fid);
        end
        
        if save_bool
            varargout{1} = prod_out;
            varargout{2} = prod_c_out;
            varargout{3} = react_out;
            varargout{4} = react_c_out;
            varargout{5} = k_out;
        end
    end

    function [reactivity, react_uncert] = calc_mean_total_react()
        % This will calculate the total organic reactivity for the
        % requested city.
        city_name = ask_multichoice('Which city to calculate for?',{'Atlanta','Birmingham','Montgomery'},'list',true);
        [city_lon, city_lat] = return_city_info(city_name);
        r_type = ask_multichoice('Which type of reactivity to calculate for?',{'NO+RO2','VOC+OH'},'list',true);
        
        
        [~, ~, all_react, ~, all_rates] = pick_wrf_rxns('NO+RO2');
        
        % Go through all the reactants and extract the names of the organic
        % compounds.
        org_cmpds = {};
        switch lower(r_type)
            case 'no+ro2'
                rej_react = 'NO';
            case 'voc+oh'
                rej_react = 'HO';
        end
        for a=1:numel(all_react)
            xx = ~strcmpi(rej_react,all_react{a});
            org_cmpds = cat(2, org_cmpds, all_react{a}(xx));
        end
        
        % Get the statistics on the concentration of each organic compound
        opts.subset_lon = city_lon;
        opts.subset_lat = city_lat;
        opts.model_layers = [1 1];
        opts.wrf_vars = lower(org_cmpds); % in the output, the chemical compounds are usually in lower case
        opts.utc_hr = 19;
        opts.start_date = '2013-06-01';
        opts.end_date = '2013-06-30';
        
        wrfstats = wrf_cell_stats(opts);
        
        % Now go through and match up concentrations, converting them from
        % ppm or ppmv to molec./cm^3.
        org_conc = nan(size(org_cmpds));
        org_std = nan(size(org_cmpds));
        for a=1:numel(org_cmpds)
            org_name = lower(org_cmpds{a});
            if ~ismember(wrfstats.(org_name).units, 'ppm','ppmv')
                E.notimplemented(sprintf('Conversion from unit %s to molec./cm^3',wrfstats.(org_name).units));
            end
            org_conc(a) = wrfstats.(org_name).mean * 1e-6 * 2e19;
            org_std(a) = wrfstats.(org_name).stddev * 1e-6 * 2e19;
        end
        
        % Sum up organic reactivity and calculate its uncertainty in
        % quadrature based on the spread of the data.
        reactivity = nansum2(org_conc .* all_rates);
        react_uncert = sqrt(nansum2((org_std .* all_rates).^2));
    end

    function [reactivity, reactivity_uncer, speciated_reactivity, species_rates] = calc_daily_total_reactivity()
        % This will calculate the total organic reactivity for the
        % requested city.
        city_name = ask_multichoice('Which city to calculate for?',{'Atlanta','Birmingham','Montgomery'},'list',true);
        [city_lon, city_lat] = return_city_info(city_name);
        r_type = ask_multichoice('Which type of reactivity to calculate for?',{'NO+RO2','VOC+OH'},'list',true);
        start_date = ask_date('Enter the starting date');
        end_date = ask_date('Enter the ending date');
        
        
        [~, ~, all_react, ~, all_rates] = pick_wrf_rxns(r_type);
        all_rates = cell2mat(all_rates);
        % Go through all the reactants and extract the names of the organic
        % compounds. Also remove those not available in the output file.
        org_cmpds = {};
        W = dir(fullfile(wrf_path,'wrfout*'));
        wi = ncinfo(fullfile(wrf_path, W(1).name));
        wrf_vars = {wi.Variables.Name};
        
        switch lower(r_type)
            case 'no+ro2'
                rej_react = 'NO';
            case 'voc+oh'
                rej_react = 'HO';
        end
        for a=1:numel(all_react)
            xx = ~strcmpi(rej_react,all_react{a});
            if sum(xx) > 1
                E.notimplemented('3 body reactions')
            end
            this_voc = lower(all_react{a}(xx));
            if ismember(upper(this_voc), wrf_vars)
                if ~ismember(this_voc, wrf_vars)
                    this_voc = upper(this_voc);
                else
                    E.callError('multiple species','The species %s can be found as both upper and lower case in a wrfout file');
                end
            elseif ~ismember(this_voc, wrf_vars)
                fprintf('\t %s is not available in wrfout\n', this_voc)
                continue
            end
            org_cmpds = cat(2, org_cmpds, this_voc);
        end
        u_org_cmpds = unique(org_cmpds);
        
        % Get the statistics on the concentration of each organic compound
        opts.subset_lon = city_lon;
        opts.subset_lat = city_lat;
        opts.model_layers = [1 1];
        opts.wrf_vars = u_org_cmpds;
        opts.utc_hr = 19;
        
        dvec = datenum(start_date):datenum(end_date);
        reactivity = nan(size(dvec));
        reactivity_uncer = nan(size(dvec));
        
        speciated_reactivity = make_empty_struct_from_cell(u_org_cmpds, zeros(size(dvec)));
        species_rates = make_empty_struct_from_cell(u_org_cmpds, 0);
        
        for d=1:numel(dvec)
            opts.start_date = datestr(dvec(d));
            opts.end_date = datestr(dvec(d));
            wrfstats = wrf_cell_stats(opts);
            
            % Now go through and match up concentrations, converting them from
            % ppm or ppmv to molec./cm^3.
            org_conc = nan(size(org_cmpds));
            org_std = nan(size(org_cmpds));
            for a=1:numel(org_cmpds)
                org_name = org_cmpds{a};
                if ~ismember(wrfstats.(org_name).units, {'ppm','ppmv'})
                    if isempty(wrfstats.(org_name).units)
                        if d==1
                            warning('No units for %s, assuming to be in PPM',org_name)
                        end
                    else
                        E.notimplemented(sprintf('Conversion from unit %s to molec./cm^3',wrfstats.(org_name).units));
                    end
                end
                org_conc(a) = wrfstats.(org_name).mean * 1e-6 * 2e19;
                org_std(a) = wrfstats.(org_name).stddev * 1e-6 * 2e19;
                speciated_reactivity.(org_name)(d) = speciated_reactivity.(org_name)(d) + org_conc(a) * all_rates(a);
                if d == 1
                    species_rates.(org_name) = species_rates.(org_name) + all_rates(a);
                end
            end
            
            reactivity(d) = nansum2(org_conc .* all_rates);
            reactivity_uncer(d) = sqrt(nansum2( (org_std .* all_rates).^2 ));
        end
    end

    function emis = compute_wrf_emis
    % This will add up the emissions from the two large enough cities in
    % the domain: both by assuming a simple 50 km radius and using
    % find_plume to identify emissions above a given threshold.
    
    city_names = {'Atlanta','Birmingham'};
    substruct = struct('within50km',0,'floodfill',0,'units','mol NO h^{-1}');
    emis = make_empty_struct_from_cell(city_names, substruct);
    % Emissions in WRF do not change by day, so just load one file at 1900
    % UTC.
    F = dir(fullfile(wrf_path,'wrfout*_19-00*'));
    [xlon, xlat, e_no] = read_wrf_vars(wrf_path, F(1), {'XLONG','XLAT','E_NO'});
    e_no = squeeze(nansum2(e_no,3));
    
    % Options for subsetting for the 50 km approach
    options.quad_bool = false;
    options.dist_limit = 50;
    
    in_50 = [];
    in_flood = false(size(xlon));
    
    [xloncorn, xlatcorn] = wrf_grid_corners(xlon,xlat);
    earth_ellip = referenceEllipsoid('wgs84','kilometers');
    for c=1:numel(city_names)
        [city_lon, city_lat] = return_city_info(city_names{c});
        options.center_lon = city_lon;
        options.center_lat = city_lat;
        in = subset_WRF_grid_cells(xlon,xlat,[-88 -80],[30 36], options);
        in_50 = [in_50; in(:)];
        city_loncorn = xloncorn(:,in);
        city_latcorn = xlatcorn(:,in);
        city_areas = nan(1,numel(in));
        for a=1:numel(in)
            city_areas(a) = areaint(city_latcorn(:,a), city_loncorn(:,a), earth_ellip);
        end
        
        
        emis.(city_names{c}).within50km = nansum2(e_no(in) .* city_areas');
        
        in = find_plume(e_no, xlon, xlat, 50, city_lon, city_lat);
        in_flood = in_flood | in;
        city_loncorn = xloncorn(:,in);
        city_latcorn = xlatcorn(:,in);
        city_areas = nan(1,sum(in(:)));
        for a=1:numel(city_areas)
            city_areas(a) = areaint(city_latcorn(:,a), city_loncorn(:,a), earth_ellip);
        end
        emis.(city_names{c}).floodfill = nansum2(e_no(in) .* city_areas');
    end
    
    % Plot which emissions were used in each case
    figure;
    in_log = false(size(e_no));
    in_log(in_50) = true;
    e_no_50 = e_no;
    e_no_50(~in_log) = nan;
    pcolor(xlon, xlat, e_no_50);
    title('50 km')
    
    figure;
    e_no_flood = e_no;
    e_no_flood(~in_flood) = nan;
    pcolor(xlon,xlat,e_no_flood);
    title('Floodfill')
    
    
    end

    function ro2_conc = compute_ss_ro2(ro2_name, wrf_files, voc_fatal)
        % This function will compute steady-state RO2 concentration for the
        % given RO2 species, assuming:
        %
        %   d[RO2]/dt = 0 = \sum_i k_{RH_i+OH}*[RH_i]*[OH] - k_{RO2+NO}*[RO2]*[NO]
        %
        % where RH_i represents all the possible VOCs that can produce the
        % desired RO2 via reaction with OH. Inputs are the RO2 name from
        % pick_wrf_rxns and a structure output from dir() with a list of
        % WRF files to consider. voc_fatal is optional, it defaults to 1
        % which means that if the required VOC cannot be found. 0 disables
        % any mention of the failure to find that VOC. -1 will print out a
        % message without erroring.
        
        % Method: find all reactions that produce the specified RO2 via
        % reaction with OH. Make sure that the VOC reacting with OH exists
        % in the WRF files (assuming the first one is representative).
        wi = ncinfo(fullfile(wrf_path,wrf_files(1).name));
        wrf_varnames = {wi.Variables.Name};
        vars_to_read = {'T','P','PB','ho','no'};
        [prods,prod_coeffs,reacts,~,rates] = pick_wrf_rxns({'One species',ro2_name,'Products'},false);
        voc_names = cell(size(reacts));
        ro2_coeffs = zeros(size(reacts));
        for a=1:numel(reacts)
            if numel(reacts{a}) > 2
                rates{a} = 0;
                continue % shouldn't be any 3 body reactions.
            end
            xx = strcmpi(reacts{a},'OH') | strcmpi(reacts{a},'HO');
            if sum(xx) == 0
                rates{a} = 0;
                continue % only want reactions with OH to produce this
            end
            voc_a = reacts{a}{~xx};
            if any(strcmp(wrf_varnames,voc_a))
                voc_names{a} = voc_a;
            elseif any(strcmp(wrf_varnames,lower(voc_a))) %#ok<STCI>
                voc_names{a} = lower(voc_a);
            elseif any(strcmp(wrf_varnames,upper(voc_a))) %#ok<STCI>
                voc_names{a} = upper(voc_a);
            else
                rates{a} = 0;
                msg = sprintf('Species "%s" from %s -> %s cannot be found in any form in the first WRF file',voc_a,strjoin(reacts,' + '),strjoin(prods,' + '));
                if voc_fatal > 0
                    E.callError('voc_not_found',msg);
                elseif voc_fatal < 0
                    fprintf(msg);
                end
            end
            
            rr = strcmpi(prods{a}, ro2_name);
            if sum(rr) ~= 1
                E.callError('ro2_count','%s not defined exactly once in the reaction %s -> %s',ro2_name,strjoin(reacts,' + '),strjoin(prods,' + '));
            end
            ro2_coeffs(a) = prod_coeffs{a}(rr);
        end
        xx1 = iscellcontents(voc_names,'isempty');
        xx2 = ro2_coeffs == 0;
        xx3 = iscellcontents(rates,@(x) ~isa(x,'function_handle'));
        if any(xor(xx1,xx2)) || any(xor(xx1,xx3))
            % If (somehow) we mismatched the VOC names and their respective
            % rate constants, this should catch it by finding any locations
            % where we can remove an empty VOC name or 0 rate constant.
            E.callError('react_rate_mismatch','There is disagreement between which VOC names, RO2 coefficiets, and rate constants to remove')
        end
        
        if sum(~xx1) == 0
            ro2_conc = 0;
            return
        end
        
        voc_names = voc_names(~xx1);
        ro2_coeffs = ro2_coeffs(~xx1);
        rates = rates(~xx3);
        
        if sum(~xx1) == 0
            ro2_conc = 0;
            return
        end
        
        vars_to_read = [vars_to_read, voc_names];
        wrf_vars = cell(size(vars_to_read));
        [wrf_vars{:}] = read_wrf_vars(wrf_path, wrf_files, vars_to_read, 0, 'visual');
        T = wrf_vars{1};
        P = wrf_vars{2};
        PB = wrf_vars{3};
        oh_conc = wrf_vars{4};
        no_conc = wrf_vars{5};
        voc_conc = wrf_vars(6:end);
        wrf_temp = convert_wrf_temperature(T,P,PB);
        wrf_ndens = calculate_wrf_air_ndens(T,P,PB);
        for a=1:numel(voc_conc)
            voc_unit = ncreadatt(fullfile(wrf_path,wrf_files(1).name),vars_to_read{a+5},'units');
            if ~isempty(voc_unit)
                conversion = convert_units(1, voc_unit, 'ppp'); % convert from usually ppm too unscaled mixing ratio
            else
                conversion = 1e-6; % if no unit given, assume ppm
            end
            voc_conc{a} = conversion * voc_conc{a} .* wrf_ndens;
        end
        
        % Calculate the rate of reaction of VOC + OH --> RO2
        sum_rate_RH_OH = zeros(size(voc_conc{1}));
        for a=1:numel(voc_conc)
            k = rates{a}(wrf_temp, wrf_ndens);
            sum_rate_RH_OH = sum_rate_RH_OH + ro2_coeffs(a) .* k .* voc_conc{a} .* oh_conc;
        end
        
        % Now find the rate for RO2 + NO 
        [~, ~, reacts, ~, rates] = pick_wrf_rxns({'One species',ro2_name,'Reactants'},false);
        for a=1:numel(reacts)
            xx = strcmpi(reacts{a},'NO');
            if sum(xx) == 1 && numel(reacts{a}) == 2
                k_ro2_no = rates{a}(wrf_temp, wrf_ndens);
            end
        end
        rate_RO2_NO = k_ro2_no .* no_conc;
        
        % Steady state concentration is the ratio of the production and
        % loss.
        ro2_conc = sum_rate_RH_OH ./ rate_RO2_NO;
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

function [city_lon, city_lat, windvel, theta, wind_dnums] = return_city_info(city_name)
switch lower(city_name)
    case 'atlanta'
        city_lon = -84.39;
        city_lat = 33.775;
        F=load('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta-Wind-Conditions-1900UTC-5layers.mat');
    case 'birmingham'
        city_lon = -86.80;
        city_lat = 33.52;
        F=load('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Birmingham-Wind-Conditions-1900UTC-5layers.mat');
    case 'montgomery'
        city_lon = -86.30;
        city_lat = 32.37;
        F=load('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Montgomery-Wind-Conditions-1900UTC-5layers.mat');
    otherwise
        E=JLLErrors;
        E.badinput('%s not recognized as a city',city_name);
end
windvel = F.windvel;
theta = F.theta;
wind_dnums = F.dnums;
end

function apr_hrs = get_apri_hr(Data)
apr_str = Data.BEHRaprioriMode;
st = regexp(apr_str,'[')+1;
ed = regexp(apr_str,']')-1;
apr_hrs_cell = strsplit(apr_str(st:ed),';');
apr_hrs = nan(size(apr_hrs_cell));
for a=1:numel(apr_hrs)
    apr_hrs(a) = str2double(apr_hrs_cell{a});
end
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