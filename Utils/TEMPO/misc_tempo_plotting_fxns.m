function [ varargout ] = misc_tempo_plotting_fxns( plottype, varargin )
%[ VARARGOUT ] = MISC_TEMPO_PLOTTING_FXNS( PLOTTYPE, ... )
%   This function contains a collection of miscellaneous plotting functions
%   for TEMPO analysis. By collection them here, this helps keep me from
%   forgetting what I called each function. The plot type is specified by a
%   string in the first argument.
%
%   MISC_TEMPO_PLOTTING_FXNS('dif') will allow you to create difference
%   plots of various quantities, such as AMFs or WRF columns. This will
%   produce a grid of subplots, one for each hour.
%       ('dif', clim) allows you to specify the color limits of all the
%       subplots using a 2-element vector for clim as you would for the
%       input to CAXIS.
%
%       ('dif', clim, pub) allows you to produce the figures individually
%       for publication (with an eye towards compositing them in another
%       program like Adobe Illustrator or GNU Inkscape). pub should be true
%       if you wish the individual plots in this way. If you wish the color
%       limits to still be determined automatically, pass an empty matrix
%       for clim.
%
%   [ AVG_VAL ] = MISC_TEMPO_PLOTTING_FXNS('avg') will produce monthly
%   averages for a given field in the Data structures for TEMPO. 
%       ('avg', apriori, start_date, end_date, fieldname) allows you to
%       specify the answers to each question ahead of time, useful if using
%       this function in some sort of automated capacity.
%
%   MISC_TEMPO_PLOTTING_FXNS('plotavg') will let you produce a difference
%   plot of the average values for certain fields over a time range. It
%   will ask for the parameters interactively.


E = JLLErrors;
switch lower(plottype)
    case 'dif'
        make_difference(varargin{:})
    case 'plotavg'
        plot_avg_difference();
    case 'avg'
        varargout{1} = make_average();
    otherwise
        fprintf('Plottype not recognized\n');
        return
end


    function make_difference(clim, pub)
        % Similar to plot_diff() in misc_behr_wind_plots, this will ask the
        % user to specify what kind of difference plot they want, then will
        % produce it.  Since TEMPO hasn't actually launched, there'll be
        % fewer options than in plot_diff. Can pass a two-element CLIM
        % vector as the first argument to set the color bar limits on all
        % plots (otherwise it will choose bounds based on the extremes and
        % set all plots to those limits). Also, send the string 'pub' as
        % the second argument to make each plot its own figure with only
        % the lat/lon limits as ticks. If you wish to plot this way with
        % automatic color limits, pass an empty matrix as the first
        % argument.
        
        if exist('clim','var') && (~isnumeric(clim) || ~isvector(clim) || numel(clim) ~= 2 || clim(1) > clim(2))
            E.badinput('If given, clim must be a 2-element numeric vector with the first element less than the second to be a valid definition of the color range')
        elseif ~exist('clim','var')
            clim = [];
        end
        
        if exist('pub','var') && strcmpi(pub,'pub')
            pub_bool = true;
        else
            pub_bool = false;
        end
        
        % There should be 10 hours per day selected for TEMPO. If not,
        % you'll need to rearrange the subplots
        utchrs = 13:22;
            
        
        city_lon = -84.39;
        city_lat = 33.755;
        
        % Ask the questions!
        try
            allowed_data = {'wrf','wind','tempo'};
            data_type = ask_multichoice('Which data do you wish to plot? WRF columns, WRF wind fields, or TEMPO AMFs?', allowed_data);
            
            if ~strcmpi(data_type,'wind')
                allowed_diffs = {'a','p','d','m'};
                diff_type = ask_multichoice('Do you want to see an absolute difference, percent difference, just the daily values, or just the monthly values?',allowed_diffs);
            else
                allowed_diffs = {'d','m'};
                diff_type = ask_multichoice('Do you want to see the daily winds or the monthly average winds?',allowed_diffs);
            end
            
            if strcmpi(data_type,'tempo')
                allowed_apriori = {'hourly','hybrid'};
                apriori_type = ask_multichoice('Select which apriori type to wish to use',allowed_apriori);
            end
        catch err
            if strcmp(err.identifier,'ask_multichoice:user_cancel')
                fprintf('Aborting plot.\n')
                return
            else
                rethrow(err)
            end
        end
        
        while true
            date_in = input('Enter the date to compare in a format Matlab can recognize: ','s');
            try 
                if strcmpi(data_type,'wrf') || strcmpi(data_type,'wind')
                    filedate = datestr(date_in, 'yyyy-mm-dd');
                else
                    filedate = datestr(date_in, 'yyyymmdd');
                end
                break
            catch err
                if strcmpi(err.identifier, 'MATLAB:datestr:ConvertToDateNumber')
                    fprintf('Date format not recognized. Try yyyy-mm-dd.\n')
                else
                    rethrow(err)
                end
            end
        end
        
        % Now load the proper data.
        if strcmpi(data_type,'wrf')
            filepath = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US_TEMPO';
            filename = 'WRF_TEMPO_%s_%s.nc';
            D = ncinfo(fullfile(filepath,'hourly',sprintf(filename,'hourly',filedate)));
            mfiledate = datestr(filedate,'yyyy-mm');
            M = ncinfo(fullfile(filepath,'monthly',sprintf(filename,'monthly',mfiledate)));
           
            daily_lon = ncread(D.Filename, 'XLONG');
            daily_lat = ncread(D.Filename, 'XLAT');
            daily_no2 = squeeze(ncread(D.Filename,'no2_ndens'));
            daily_zlev = squeeze(ncread(D.Filename,'zlev'));
            daily_tplev = find_wrf_tropopause(D);
            
            monthly_lon = ncread(M.Filename, 'XLONG');
            monthly_lat = ncread(M.Filename, 'XLAT');
            monthly_no2 = ncread(M.Filename, 'no2_ndens');
            monthly_zlev = ncread(M.Filename, 'zlev');
            monthly_tplev = find_wrf_tropopause(M);
            
            xx = 16:60;
            yy = 19:59;
            
            daily_lon = daily_lon(xx,yy);
            daily_lat = daily_lat(xx,yy);
            daily_no2 = daily_no2(xx,yy,:,:);
            daily_zlev = daily_zlev(xx,yy,:,:);
            daily_tplev = daily_tplev(xx,yy,:,:);
            for a=1:size(daily_no2,1)
                for b=1:size(daily_no2,2)
                    for h=1:size(daily_no2,4)
                        tp = daily_tplev(a,b,h);
                        if tp>0
                            daily_no2(a,b,tp:end,h) = nan;
                        end
                    end
                end
            end
            daily_data = nansum2(daily_no2 .* daily_zlev*100,3);
            
            monthly_no2 = monthly_no2(xx,yy,:,:);
            monthly_zlev = monthly_zlev(xx,yy,:,:);
            monthly_tplev = monthly_tplev(xx,yy,:,:);
            for a=1:size(monthly_no2,1)
                for b=1:size(monthly_no2,2)
                    for h=1:size(monthly_no2,4)
                        tp = monthly_tplev(a,b,h);
                        if tp>0
                            monthly_no2(a,b,tp:end,h) = nan;
                        end
                    end
                end
            end
            monthly_data = nansum2(monthly_no2 .* monthly_zlev*100,3);
            
            % Used for labeling the color bar
            data_description = 'VCD (molec. cm^{-2})';
            
        elseif strcmpi(data_type,'wind')
            filepath = '/Volumes/share2/USERS/LaughnerJ/WRF/E_US_TEMPO';
            filename = 'WRF_TEMPO_%s_%s.nc';
            if strcmpi(diff_type,'d')
                F = ncinfo(fullfile(filepath,'hourly',sprintf(filename,'hourly',filedate)));
            elseif strcmpi(diff_type,'m')
                mfiledate = datestr(filedate,'yyyy-mm');
                F = ncinfo(fullfile(filepath,'monthly',sprintf(filename,'monthly',mfiledate)));
            else
                E.notimplemented('diff type = %s', diff_type)
            end
            lon = ncread(F.Filename, 'XLONG');
            lat = ncread(F.Filename, 'XLAT');
            U = squeeze(ncread(F.Filename,'U'));
            V = squeeze(ncread(F.Filename,'V'));
            COSALPHA = squeeze(ncread(F.Filename, 'COSALPHA'));
            SINALPHA = squeeze(ncread(F.Filename, 'SINALPHA'));
            [U,V] = wrf_winds_transform(U,V,COSALPHA,SINALPHA);
            
            xx = 16:60;
            yy = 19:59;
            
            % calculate surface wind speed around atlanta
            [~,atlanta] = min((lon(:) - -84.39).^2 + (lat(:) - 33.775).^2);
            [at_x, at_y] = ind2sub(size(lon),atlanta);
            avg_wind = nan(size(utchrs));
            for h=1:numel(utchrs)
                U_at = U((at_x-1):(at_x+1), (at_y-1):(at_y)+1, 1, h);
                V_at = V((at_x-1):(at_x+1), (at_y-1):(at_y)+1, 1, h);
                Ubar = nanmean(U_at(:));
                Vbar = nanmean(V_at(:));
                avg_wind(h) = nanmean(sqrt(Ubar.^2 + Vbar.^2));
            end
            
            % Not used for winds
            data_description = '';
        elseif strcmpi(data_type,'tempo')
            cap_apriori = apriori_type;
            cap_apriori(1) = upper(cap_apriori(1));
            filepath = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/Atlanta Tempo %s - No clouds/';
            daily_path = sprintf(filepath, cap_apriori);
            filename = sprintf('TEMPO_BEHR-SIM_ATLANTA_%s.mat',filedate);
            D = load(fullfile(daily_path, filename));
            
            monthly_path = sprintf(filepath, 'Monthly');
            M = load(fullfile(monthly_path, filename));
            
            % Append along 4th dim to be like WRF
            daily_lon = D.Data(1).Longitude;
            daily_lat = D.Data(1).Latitude;
            daily_data = D.Data(1).BEHRAMFTrop;
            daily_data = repmat(daily_data,1,1,1,numel(D.Data));
            monthly_data = M.Data(1).BEHRAMFTrop;
            monthly_data = repmat(monthly_data,1,1,1,numel(D.Data));
            
            for a=1:numel(D.Data)
                daily_data(:,:,1,a) = D.Data(a).BEHRAMFTrop;
                monthly_data(:,:,1,a) = M.Data(a).BEHRAMFTrop;
            end
            data_description = 'AMF';
        end
        
        
        maxval = 0;
        minval = 0;
        if ~pub_bool
            figure;
        else
            f = gobjects(size(utchrs));
        end
        for h=1:numel(utchrs)
            if ~pub_bool
                subplot(5,2,h)
            else
                f(h) = figure;
            end
            if ~strcmpi(data_type,'wind')
                if strcmpi(diff_type,'p')
                    datadiff = (daily_data(:,:,:,h) - monthly_data(:,:,:,h))./monthly_data(:,:,:,h) * 100;
                    datadiff = squeeze(datadiff);
                    cblabel = sprintf('%%\\Delta %s',data_description);
                elseif strcmpi(diff_type,'a')
                    datadiff = daily_data(:,:,:,h) - monthly_data(:,:,:,h);
                    datadiff = squeeze(datadiff);
                    cblabel = sprintf('\\Delta %s',data_description);
                elseif strcmpi(diff_type,'d')
                    datadiff = daily_data(:,:,:,h);
                    datadiff = squeeze(datadiff);
                    cblabel = sprintf('Daily %s', data_description);
                elseif strcmpi(diff_type,'m')
                    datadiff = monthly_data(:,:,:,h);
                    datadiff = squeeze(datadiff);
                    cblabel = sprintf('Monthly %s', data_description);
                end
                maxval = max(maxval, max(datadiff(:)));
                minval = min(minval, min(datadiff(:)));
                pcolor(daily_lon, daily_lat, datadiff);
                shading flat;
                colormap('jet')
            else
                quiver(lon(xx,yy),lat(xx,yy),U(xx,yy,1,h),V(xx,yy,1,h));
                xl = get(gca,'xlim');
                yl = get(gca,'ylim');
                text(xl(1) + abs(diff(xl)*0.03), yl(1) + abs(diff(yl)* 0.03), sprintf('Avg. wind = %.1f m/s', avg_wind(h)));
            end
            if ~pub_bool
                msize = 9;
                lwidth = 1;
            else
                msize = 16;
                lwidth = 2.5;
            end
            line(city_lon, city_lat, 'marker','p','color','k','linewidth',lwidth,'markersize',msize);
            if pub_bool
                set(gca,'xtick',[-87 -81]);
                set(gca,'ytick',[32.5 35.5]);
                set(gca,'fontsize',18);
            elseif ~strcmpi(data_type,'wind')
                cb=colorbar;
                cb.Label.String = cblabel;
            end
            title(sprintf('UTC Hour %d',utchrs(h)));
        end
        
        if ismember(diff_type,{'a','p'});
            truemax = max(abs(maxval), abs(minval));
            maxval = truemax;
            minval = -truemax;
        end
        
        for h=1:numel(utchrs)
            if ~pub_bool
                subplot(5,2,h);
            else
                figure(f(h));
            end
            if isempty(clim)
                caxis([minval, maxval]);
            else
                caxis(clim);
            end
        end
        
        if pub_bool && ~strcmpi(data_type,'wind')
            % Create figure with just the colorbar
            figure;
            cb = colorbar('southoutside');
            axis off
            cb.Label.String = cblabel;
            cb.Position = [0.1 0.5 0.8 0.4]; % fill most of the figure with room for ticks and label
            set(gcf,'units','centimeters');
            pos = get(gcf,'position');
            set(gcf,'position',[pos(1:2), 10, 2]);
            caxis(clim);
            colormap('jet')
            set(gca,'fontsize',14)
            cb.AxisLocation = 'in';
            fprintf('To get the proper size .png, it seems to work best using the menu saveas option\n');
        end
    end

    function plot_avg_difference()
        % Get as input which a priori case we are averaging over and start
        % and end dates. These can also be passed as arguments to allow
        % this function to be used by others within this file.
        allowed_apriori = {'hourly','hybrid','monthly'};
        try
            new_apriori = ask_multichoice(sprintf('Which a priori do you want to use as the "new" apriori?\n'), allowed_apriori, 'default', 'hybrid');
            old_apriori = ask_multichoice(sprintf('Which a priori do you want to use as the base apriori?\n'), allowed_apriori, 'default', 'monthly');
            start_date = ask_date('Enter the date for the beginning of the averaging period (yyyy-mm-dd)');
            end_date = ask_date('Enter the end of the averaging period (yyyy-mm-dd)');
        catch err
            if strcmp(err.identifier, 'ask_multichoice:user_cancel')
                return
            else
                rethrow(err)
            end
        end
        % Assuming AMFs for now, but could expand this to other fields
        % later
        allowed_fields = {'amf'};
        field = ask_multichoice('What field do you want to compare averages for?', allowed_fields, 'default', 'amf');
        switch field
            case 'amf'
                fieldname = 'BEHRAMFTrop';
        end
        
        % Percent or absolute difference?
        allowed_difftypes = {'p','a'};
        difftype = ask_multichoice('Do you want a percent or absolute difference?', allowed_difftypes);
        
        % Use this to define the colorbar labels. Expand if more diff types
        % or fields added (fields down dim 1, diff types across dim 2)
        cblabelcell = {'%\Delta AMF', '\Delta AMF'};
        xx_label = strcmp(field, allowed_fields);
        yy_label = strcmp(difftype, allowed_difftypes);
        
        [new_avg, lon, lat] = make_average(new_apriori, start_date, end_date, fieldname);
        old_avg = make_average(old_apriori, start_date, end_date, fieldname);
        
        if strcmp(difftype,'a')
            dif = new_avg - old_avg;
        elseif strcmp(difftype,'p')
            dif = (new_avg - old_avg)./old_avg * 100;
        end
        
        figure;
        cmap = 'jet';%flipud(colormap('jet'));
        mlim = 0;
        utchr = 13:22;
        for a=1:size(dif,3)
            subplot(5,2,a)
            pcolor(lon, lat, dif(:,:,a));
            title(sprintf('UTC %04d',utchr(a)*100));
            shading flat
            cb = colorbar;
            colormap(cmap)
            cb.Label.String = cblabelcell{xx_label,yy_label};
            set(gca,'fontsize',16)
            mlim = max(max(abs(cb.Limits)),mlim);
        end
        for a=1:size(dif,3)
            subplot(5,2,a)
            caxis([-mlim, mlim])
        end
    end

    function [avg_mat, lon, lat] = make_average(apriori, start_date, end_date, fieldname)
        % initialize output in case user cancels and we need to return
        % SOMETHING
        avg_mat = []; lon = []; lat = [];
        
        % Get as input which a priori case we are averaging over and start
        % and end dates, as well as what field to average. These can also
        % be passed as arguments to allow this function to be used by
        % others within this file.
        try
            if ~exist('apriori','var')
                allowed_apriori = {'hourly','hybrid','monthly'};
                apriori = ask_multichoice('Which a priori do you want to use?', allowed_apriori);
            end

            if ~exist('start_date','var')
                start_date = ask_date('Enter the date for the beginning of the averaging period (yyyy-mm-dd)');
            end
            if ~exist('end_date','var')
                end_date = ask_date('Enter the end of the averaging period (yyyy-mm-dd)');
            end
        catch err
            if strcmp(err.identifier, 'ask_multichoice:user_cancel')
                return
            else
                rethrow(err)
            end
        end
        
        % (we need to load at least one file to get the fields available)
        switch apriori
            case 'hourly'
                p='Atlanta TEMPO Hourly - No clouds';
            case 'hybrid'
                p='Atlanta TEMPO Hybrid - No clouds';
            case 'monthly'
                p='Atlanta TEMPO Monthly - No clouds';
        end
        fpath = fullfile('/','Users','Josh','Documents','MATLAB','BEHR','Workspaces','Wind speed', p);
        F = dir(fullfile(fpath,'TEMPO_BEHR-SIM*.mat'));
        D = load(fullfile(fpath, F(1).name));
        Data = D.Data;
        fns = fieldnames(Data);
        
        if ~exist('fieldname','var')
            input('Choose the field to average from the following list (enter to continue):','s');
            for a=1:numel(fns)
                fprintf('\t%d - %s\n',a,fns{a});
            end
            fn_num = input(sprintf('Enter selection between 1 and %d: ', numel(fns)),'s');
            
            while true
                if strcmpi(fn_num,'q')
                    fprintf('Aborting\n')
                    return
                else
                    fn_num = str2double(fn_num);
                    if fn_num >= 1 || fn_num <= numel(fns)
                        break
                    end
                end
            end
            
            fieldname = fns{fn_num};
        end
            
        
        sdate = datenum(start_date);
        edate = datenum(end_date);
        
        % Because the TEMPO pixels are always the same size and location
        % (as far as my simulation goes!) they will overlap exactly and we
        % don't need to worry about weighting them in any way.
        
        for a=1:numel(F)
            [s,e] = regexp(F(a).name,'\d\d\d\d\d\d\d\d');
            fdate = datenum(F(a).name(s:e),'yyyymmdd');
            if fdate >= sdate && fdate <= edate
                D=load(fullfile(fpath,F(a).name));
                this_mat = cat(3,D.Data.(fieldname));
                avg_mat = cat(4, avg_mat, this_mat);
            end
        end
        
        avg_mat = nanmean(avg_mat,4);
        lon = D.Data(1).Longitude;
        lat = D.Data(1).Latitude;
    end

end

