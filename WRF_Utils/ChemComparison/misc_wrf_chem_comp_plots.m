function [ varargout ] = misc_wrf_chem_comp_plots( plttype, varargin )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
%
% Many of these subfunctions allow you to call them with a structure to
% skip the questions used to subset the data. This structure can have any
% or all of the following fields:
%
%       match_file - the file name containing the MatchDC3 structure to
%       load. These must be located in the directory specified by the
%       matchdir variable in this function.
%
%       quantity - which quantity to plot. Can be 'NOx', 'NO2', 'MPN',
%       'HNO3', 'NO2:NO', 'lnox_total', 'PHOTR_NO2'. Note that some may not
%       be in a particular match file.
%
%       lats_filter - can be one of 'all', 'midlats', 'subtrop'. Midlats
%       are > 35 N only, subtrop 20--35 N.
%
%       strat_filter - must be a scalar logical. If true, removes aircraft
%       data with O3/CO < 1.5 which is indicative of stratospheric
%       influence.
%
%       fresh_filter - must be a scalar logical. If true, removes aircraft
%       data with NOx / HNO3 < 5 which is indicative of fresh (unaged)
%       emissions.

E = JLLErrors;
matchdir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/WRF/Chem Comparison';

switch lower(plttype)
    case 'grand-scatter'
        grand_scatterplot(varargin{:});
    case 'grand-profile'
        if nargout == 0
            pres_binned_grand_profile(varargin{:});
        else
            varargout = cell(1,6);
            [varargout{:}] = pres_binned_grand_profile(varargin{:});
        end
    case 'make-hybrid'
        [varargout{1}, varargout{2}] = make_hybrid_profiles(varargin{:});
    case 'save-hybrid'
        save_hybrid_profs();
    case 'amf-sens'
        if nargout > 0
            varargout{1} = amf_test(nargout);
        else
            amf_test(nargout);
        end
    case 'scatter-diff'
        scatter_by_pres_time
    case 'timeser'
        avg2gridcell_timeser
    case 'two-hist'
        two_histograms()
    case 'flashcounts'
        figs = compare_flashcounts(varargin{:});
        if nargout > 0
            varargout = {figs};
        end
    case 'make-lightning-met-profs'
        make_lightning_met_profs();
    case 'plot-lightning-met-profs'
        figs = plot_lightning_temp_diffs(varargin{:});
        if nargout > 0
            varargout = {figs};
        end
    case 'bc-test'
        compare_geos_moz_bcs()
end

    function varargout = grand_scatterplot(Opts)
        if ~exist('Opts','var')
            Opts = struct;
        end
        Match = get_match_struct(Opts);
        [air_data, wrf_data, xstr] = read_data(Match, Opts);
        wrf_xx = unique_wrf_element(Match);
        lats_xx = subset_latitudes(Match, Opts);
        trop_xx = filter_strat_influence(Match, Opts);
        fresh_xx = filter_fresh_lnox(Match, Opts);
        
        data_log = lats_xx & trop_xx & fresh_xx;
        
        if nargout == 0
            figure;
            l=gobjects(2,1);
            l(1) = line(air_data(data_log), Match.data.pres(data_log), 'color', [0.5 0.5 0.5], 'linestyle', 'none', 'marker', '.');
            % remove duplicate WRF data points to declutter the scatter
            % plots somewhat
            l(2) = line(wrf_data(wrf_xx & data_log), Match.wrf.pres(wrf_xx & data_log), 'color', 'r', 'linestyle', 'none', 'marker', '.');
            
            legend(l, {'Aircraft Data', 'WRF'});
            set(gca,'fontsize',16);
            xlabel(xstr);
            ylabel('Pressure (hPa)');
            set(gca,'ydir','reverse');
        else
            varargout{1} = air_data(data_log);
            varargout{2} = Match.data.pres(data_log);
            varargout{3} = wrf_data(data_log); % removed wrf_xx here because we want the sampling statistics to be the same in the output
            if isfield(Match.wrf, 'pres')
                varargout{4} = Match.wrf.pres(data_log);
            else
                varargout{4} = (Match.wrf.P + Match.wrf.PB)/100;
            end
            varargout{5} = xstr;
        end
    end

    function [air_data_bins, air_data_means, air_pres_bins, wrf_data_bins, wrf_data_means, wrf_pres_bins] = pres_binned_grand_profile(Opts)
        if ~exist('Opts','var')
            Opts = struct;
            do_error_bars = ask_yn('Include error bars on the means?');
        else
            do_error_bars = false;
        end
        [air_data, air_pres, wrf_data, wrf_pres, xstr] = grand_scatterplot(Opts);
        [air_data_bins, air_pres_bins, air_quants] = bin_omisp_pressure(air_pres, air_data);
        [wrf_data_bins, wrf_pres_bins, wrf_quants] = bin_omisp_pressure(wrf_pres, wrf_data);
        [air_data_means, ~, air_data_std] = bin_omisp_pressure(air_pres, air_data, 'mean');
        [wrf_data_means, ~, wrf_data_std] = bin_omisp_pressure(wrf_pres, wrf_data, 'mean');
        
        if nargout == 0
            figure;
            l=gobjects(4,1);
            l(1) = line(air_data_bins, air_pres_bins, 'color', 'b', 'linestyle', '-', 'marker', 'o');
            l(2) = line(wrf_data_bins, wrf_pres_bins, 'color', 'r', 'linestyle', '-', 'marker', 'v');
            l(3) = line(air_data_means, air_pres_bins, 'color', 'b', 'linestyle', '--', 'marker', 's');
            l(4) = line(wrf_data_means, wrf_pres_bins, 'color', 'r', 'linestyle', '--', 'marker', '^');

            if do_error_bars
                scatter_errorbars(air_data_means, air_pres_bins, air_data_std, 'color', 'b');
                scatter_errorbars(wrf_data_means, wrf_pres_bins, wrf_data_std, 'color', 'r');
            end
            
            legend(l, {'Aircraft Data Med.', 'WRF Med.', 'Aircraft Mean','WRF Mean'});
            set(gca,'fontsize',16);
            xlabel(xstr);
            ylabel('Pressure (hPa)');
            set(gca,'ydir','reverse');
        end
    end

    function save_hybrid_profs()
        save_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/LNOx-AMFs/Profiles';
        match_files = {'DC3-Comparison-no_lnox-fixedBC-using_wrfgridcorners.mat',...
            'DC3-Comparison-flashrate-1-molflash-500-iccg-2-newprofile-fixedBC-using_wrfgridcorners.mat',...
            'DC3-Comparison-500mol-nudge-test.mat',...
            'DC3-Comparison-flashrate-1-molflash-665-iccg-2-newprofile-fixedBC-using_wrfgridcorners.mat'};
        save_files = {'DC3-WRF-0-Hybrid-Profiles.mat', 'DC3-WRF-500-Hybrid-Profiles.mat',...
            'DC3-WRF-500-nudge-Hybrid-Profiles.mat', 'DC3-WRF-665-Hybrid-Profiles.mat'};
        
        Opts.quantity = 'NO2';
        Opts.lats_filter = 'all';
        Opts.strat_filter = false;
        Opts.fresh_filter = false;
        for a=1:numel(match_files)
            Opts.match_file = match_files{a};
            profs = make_hybrid_profiles(Opts);
            save(fullfile(save_dir, save_files{a}), 'profs');
        end
    end

    function [profs, iswrf] = make_hybrid_profiles(varargin)
        [ ~, dc3_prof, dc3_pres, ~, wrf_prof, wrf_pres] = pres_binned_grand_profile(varargin{:});
        profs.dc3_ft_wrf_bl = dc3_prof;
        iswrf.dc3_ft_wrf_bl = false(size(dc3_prof));
        iswrf.dc3_ft_wrf_bl(wrf_pres > 750) = true;
        profs.dc3_ft_wrf_bl(dc3_pres > 750) = wrf_prof(iswrf.dc3_ft_wrf_bl);
        
        profs.wrf_ft_dc3_bl = wrf_prof;
        iswrf.wrf_ft_dc3_bl = true(size(wrf_prof));
        iswrf.wrf_ft_dc3_bl(dc3_pres > 750) = false;
        profs.wrf_ft_dc3_bl(wrf_pres > 750) = dc3_prof(~iswrf.wrf_ft_dc3_bl);
        
        profs.dc3_mt_wrf_bl_ut = wrf_prof;
        iswrf.dc3_mt_wrf_bl_ut = true(size(wrf_prof));
        iswrf.dc3_mt_wrf_bl_ut(wrf_pres < 750 & wrf_pres > 375) = false;
        profs.dc3_mt_wrf_bl_ut(~iswrf.dc3_mt_wrf_bl_ut) = dc3_prof(~iswrf.dc3_mt_wrf_bl_ut);
        
        profs.wrf_mt_dc3_bl_ut = dc3_prof;
        iswrf.wrf_mt_dc3_bl_ut = false(size(dc3_prof));
        iswrf.wrf_mt_dc3_bl_ut(dc3_pres < 750 & dc3_pres > 375) = true;
        profs.wrf_mt_dc3_bl_ut(iswrf.wrf_mt_dc3_bl_ut) = wrf_prof(iswrf.wrf_mt_dc3_bl_ut);
        
        profs.dc3_rl_wrf_bl_ft = wrf_prof;
        iswrf.dc3_rl_wrf_bl_ft = true(size(wrf_prof));
        iswrf.dc3_rl_wrf_bl_ft(wrf_pres >= 700 & wrf_pres <= 875) = false;
        profs.dc3_rl_wrf_bl_ft(~iswrf.dc3_rl_wrf_bl_ft) = dc3_prof(~iswrf.dc3_rl_wrf_bl_ft);
        
        profs.dc3_prof = dc3_prof;
        profs.dc3_pres = dc3_pres;
        profs.wrf_prof = wrf_prof;
        profs.wrf_pres = wrf_pres;
    end

    function varargout = amf_test(req_var_out, varargin)
        % For a given run, this will compute the AMF for four cases:
        %  1) DC3 profile
        %  2) WRF profile
        %  3) DC3 boundary layer, WRF free trop
        %  4) WRF boundary layer, DC3 free trop
        if numel(varargin) > 0
            Opts = varargin{1};
        end
        Opts.quantity = 'NO2'; % For the AMF tests, we should ONLY be using the NO2 profile
        
        profs = make_hybrid_profiles(Opts);
        % This were separate variables
        dc3_prof = profs.dc3_prof;
        dc3_pres = profs.dc3_pres;
        wrf_prof = profs.wrf_prof;
        wrf_pres = profs.wrf_pres;
        dc3_ft_wrf_bl = profs.dc3_ft_wrf_bl;
        wrf_ft_dc3_bl = profs.wrf_ft_dc3_bl;
        dc3_mt_wrf_bl_ut = profs.dc3_mt_wrf_bl_ut;
        dc3_rl_wrf_bl_ft = profs.dc3_rl_wrf_bl_ft;
        
        if any(dc3_pres ~= wrf_pres)
            E.callError('different_pres','DC3 and WRF profiles on different pressure levels');
        end
        
        figure; subplot(3,2,1);
        line(wrf_prof, wrf_pres, 'color', 'r', 'marker', 'x');
        set(gca,'ydir','reverse','fontsize',14);
        ylabel('Pressure (hPa)'); title('WRF avg. profile');
        
        subplot(3,2,2)
        line(dc3_prof, dc3_pres, 'color', 'b', 'marker', 's');
        set(gca,'ydir','reverse','fontsize',14);
        title('DC3 avg. profile')
        
        subplot(3,2,3)
        line(wrf_prof(wrf_pres > 750 | wrf_pres < 375), wrf_pres(wrf_pres > 750 | wrf_pres < 375), 'color', 'r', 'marker', 'x');
        line(dc3_prof(wrf_pres < 750 & wrf_pres > 375), dc3_pres(wrf_pres < 750 & wrf_pres > 375), 'color', 'b', 'marker', 's');
        line(dc3_mt_wrf_bl_ut, wrf_pres, 'color', 'k', 'linestyle', '--');
        set(gca,'ydir','reverse','fontsize',14);
        ylabel('Pressure (hPa)'); title('WRF avg. profile w/ DC3 mid trop');
        
        subplot(3,2,4)
        line(wrf_prof(wrf_pres > 875 | wrf_pres < 700), wrf_pres(wrf_pres > 875 | wrf_pres < 700), 'color', 'r', 'marker', 'x');
        line(dc3_prof(wrf_pres <= 875 & wrf_pres >= 700), dc3_pres(wrf_pres <= 875 & wrf_pres >= 700), 'color', 'b', 'marker', 's');
        line(dc3_rl_wrf_bl_ft, wrf_pres, 'color', 'k', 'linestyle', '--');
        set(gca,'ydir','reverse','fontsize',14);
        ylabel('Pressure (hPa)'); title('WRF avg. profile w/ DC3 mid trop');
        
        subplot(3,2,5)
        line(dc3_prof(dc3_pres > 750), dc3_pres(dc3_pres > 750), 'color', 'b', 'marker', 's');
        line(wrf_prof(wrf_pres < 750), wrf_pres(wrf_pres < 750), 'color', 'r', 'marker', 'x');
        line(wrf_ft_dc3_bl, wrf_pres, 'color', 'k', 'linestyle', '--'); % just a check that we are plotting the same profile
        set(gca,'ydir','reverse','fontsize',14);
        xlabel('[NO_2] (ppmv)'); ylabel('Pressure (hPa)'); title('WRF FT, DC3 BL');
        
        subplot(3,2,6)
        line(dc3_prof(dc3_pres < 750), dc3_pres(dc3_pres < 750), 'color', 'b', 'marker', 's');
        line(wrf_prof(wrf_pres > 750), wrf_pres(wrf_pres > 750), 'color', 'r', 'marker', 'x');
        line(dc3_ft_wrf_bl, dc3_pres, 'color', 'k', 'linestyle', '--'); % just a check that we are plotting the same profile
        set(gca,'ydir','reverse','fontsize',14);
        xlabel('[NO_2] (ppmv)'); title('DC3 FT, WRF BL');
        
        
        % Alternate figure for group meeting
        figure;
        subplot(1,4,1);
        line(wrf_prof, wrf_pres, 'color', 'r', 'marker', 'x');
        set(gca,'ydir','reverse','fontsize',14);
        xlabel('[NO_2] (ppmv)'); ylabel('Pressure (hPa)'); title('WRF avg. profile');
       
        subplot(1,4,2);
        line(wrf_prof(wrf_pres > 875), wrf_pres(wrf_pres > 875), 'color', 'r', 'marker', 'x');
        line(wrf_prof(wrf_pres < 700), wrf_pres(wrf_pres < 700), 'color', 'r', 'marker', 'x');
        line(dc3_prof(wrf_pres <= 875 & wrf_pres >= 700), dc3_pres(wrf_pres <= 875 & wrf_pres >= 700), 'color', 'b', 'marker', 's');
        line(dc3_mt_wrf_bl_ut, wrf_pres, 'color', 'k', 'linestyle', '--');
        set(gca,'ydir','reverse','fontsize',14);
        xlabel('[NO_2] (ppmv)'); title('WRF avg. profile w/ DC3 mid trop');
        
        subplot(1,4,3);
        line(wrf_prof(wrf_pres > 750), wrf_pres(wrf_pres > 750), 'color', 'r', 'marker', 'x');
        line(wrf_prof(wrf_pres < 375), wrf_pres(wrf_pres < 375), 'color', 'r', 'marker', 'x');
        line(dc3_prof(wrf_pres < 750 & wrf_pres > 375), dc3_pres(wrf_pres < 750 & wrf_pres > 375), 'color', 'b', 'marker', 's');
        line(dc3_mt_wrf_bl_ut, wrf_pres, 'color', 'k', 'linestyle', '--');
        set(gca,'ydir','reverse','fontsize',14);
        xlabel('[NO_2] (ppmv)'); title('WRF avg. profile w/ DC3 mid trop');
        
        subplot(1,4,4);
        l(1)=line(dc3_prof(dc3_pres < 750), dc3_pres(dc3_pres < 750), 'color', 'b', 'marker', 's');
        l(2)=line(wrf_prof(wrf_pres > 750), wrf_pres(wrf_pres > 750), 'color', 'r', 'marker', 'x');
        set(gca,'ydir','reverse','fontsize',14);
        xlabel('[NO_2] (ppmv)'); title('DC3 FT, WRF BL');
        legend(l', {'DC3 data', 'WRF Data'})
        
        % Assume profiles at the center of the US, in June
        us_clon = -95;
        us_clat = 37.5;
        warn_state = warning('off', 'all'); % Disable nan's warning in intepPr2
        [wrf_amfs, szas, vzas, raas, albs, surfps] = amfSensitivityTest(wrf_prof, wrf_pres, us_clon, us_clat, 6);
        dc3_amfs = amfSensitivityTest(dc3_prof, dc3_pres, us_clon, us_clat, 6);
        dc3_mt_wrf_bl_ut_amfs = amfSensitivityTest(dc3_mt_wrf_bl_ut, wrf_pres, us_clon, us_clat, 6);
        wrf_ft_dc3_bl_amfs = amfSensitivityTest(wrf_ft_dc3_bl, wrf_pres, us_clon, us_clat, 6);
        dc3_ft_wrf_bl_amfs = amfSensitivityTest(dc3_ft_wrf_bl, dc3_pres, us_clon, us_clat, 6);
        dc3_rl_wrf_bl_ft_amfs = amfSensitivityTest(dc3_rl_wrf_bl_ft, wrf_pres, us_clon, us_clat, 6);
        warning(warn_state);
        
        % Print out everything
        amf_cell = cell(numel(wrf_amfs), 11);
        for a = 1:numel(wrf_amfs)
            [x1,x2,x3,x4,x5] = ind2sub(size(wrf_amfs), a);
            amf_cell(a,:) = {szas(x1), vzas(x2), raas(x3), albs(x4), surfps(x5), wrf_amfs(a), dc3_amfs(a), dc3_rl_wrf_bl_ft_amfs(a), wrf_ft_dc3_bl_amfs(a), dc3_mt_wrf_bl_ut_amfs(a), dc3_ft_wrf_bl_amfs(a)};
        end
        
        amf_table = cell2table(amf_cell, 'VariableNames', {'SZA','VZA','RAA','ALB','SurfP','WRF_AMF','DC3_AMF','DC3_ResidLayer_WRF_BL_FT_AMFS','WRF_FT_DC3_BL_AMF','WRF_BL_UT_DC3_MT_AMF','DC3_FT_WRF_BL_AMF'});
        if req_var_out == 0
            % print the table
            amf_table
        else
            varargout{1} = amf_table;
        end
    end

    function avg2gridcell_timeser()
        Match = get_match_struct();
        [air_data, wrf_data, xstr] = read_data(Match);
        wrf_xx = unique_wrf_element(Match);
        
        % Average the aircraft data to wrf resolution (time and space)
        air_avg_data = nan(size(wrf_xx));
        air_std_data = nan(size(wrf_xx));
        air_pres = nan(size(wrf_xx));
        air_time = nan(size(wrf_xx));
        wrf_avg_data = wrf_data(wrf_xx);
        wrf_pres = Match.wrf.pres(wrf_xx);
        for a=1:numel(wrf_xx)
            i = wrf_xx(a);
            if a < numel(wrf_xx)
                j = wrf_xx(a+1)-1;
            else
                j = numel(wrf_xx);
            end
            air_avg_data(a) = nanmean(air_data(i:j));
            air_std_data(a) = nanstd(air_data(i:j));
            air_pres(a) = nanmean(Match.data.pres(i:j));
            air_time(a) = nanmean(Match.data.datevec(i:j));
        end
        
        l = gobjects(2,1);
        figure; 
        hold on
        l(1) = scatter(air_time, air_avg_data(:) - wrf_avg_data(:), [], air_pres);
%        scatter_errorbars(air_time, air_avg_data, air_std_data, 'color', 'b')
%        l(2) = scatter(air_time, wrf_avg_data, [], wrf_pres, 'marker', 'o');
%        legend(l, {'Aircraft', 'WRF'});
        ylabel(sprintf('%s (ppmv)', xstr));
        datetick('x');
    end

    function scatter_by_pres_time()
        Match = get_match_struct();
        [air_data, wrf_data, xstr] = read_data(Match);
        difftype = ask_multichoice('Absolute or percent difference?',{'a','p'});
        absdiff_bool = strcmpi(difftype, 'a');
        
        time_inds = unique(Match.indicies.time);
        time_inds(isnan(time_inds)) = [];
        
        air_wrf_avgdiff = nan(size(air_data));
        time_vals = nan(size(air_data));
        pres_vals = nan(size(air_data));
        
        i=1;
        for a=1:numel(time_inds)
            tt = Match.indicies.time == time_inds(a);
            
            air_subset = air_data(tt);
            wrf_subset = wrf_data(tt);
            wrf_pres_subset = Match.wrf.pres(tt);
            bt_inds = unique(Match.indicies.bottom_top(tt));
            
            for b=1:numel(bt_inds)
                bb = Match.indicies.bottom_top(tt) == bt_inds(b);
                time_vals(i) = Match.wrf.time(time_inds(a));
                pres_vals(i) = nanmean(wrf_pres_subset(bb));
                if absdiff_bool
                    air_wrf_avgdiff(i) = nanmean(wrf_subset(bb)) - nanmean(air_subset(bb)); 
                else
                    air_wrf_avgdiff(i) = reldiff(nanmean(wrf_subset(bb)), nanmean(air_subset(bb)))*100; 
                end
                i=i+1;
            end
        end
        
        time_vals(i:end) = [];
        pres_vals(i:end) = [];
        air_wrf_avgdiff(i:end) = [];
        
        clim = [-max(abs(air_wrf_avgdiff)), max(abs(air_wrf_avgdiff))];
        
        figure;
        scatter(time_vals, pres_vals, [], air_wrf_avgdiff, 'filled', 'marker', 's');
        colormap(blue_red_only_cmap);
        caxis(clim);
        cb=colorbar;
        if absdiff_bool
            cb.Label.String = sprintf('%s (WRF - Obs, ppmv)', xstr);
        else
            cb.Label.String = sprintf('%s ((WRF - Obs)/Obs * 100)', xstr);
        end
        set(gca,'fontsize',14,'ydir','reverse');
        datetick('x');
    end

    function two_histograms()
        Match = get_match_struct();
        [air_data, wrf_data, chemstr] = read_data(Match);
        xx = unique_wrf_element(Match);
        wrf_xx = false(size(wrf_data));
        wrf_xx(xx) = true;
        
        bottom_pres = ask_number('Enter the bottom pressure to consider', 'default', 1020);
        top_pres = ask_number(sprintf('Enter the top pressure to consider (< %f)', bottom_pres), 'default', 0, 'testfxn', @(p) p < bottom_pres, 'testmsg', sprintf('Must be < %f', bottom_pres));
        
        air_pp = Match.data.pres <= bottom_pres & Match.data.pres >= top_pres;
        wrf_pp = Match.wrf.pres <= bottom_pres & Match.wrf.pres >= top_pres;
        
        low_bin = min([min(air_data(air_pp)), min(wrf_data(wrf_pp))]);
        high_bin = max([max(air_data(air_pp)), max(wrf_data(wrf_pp))]);
        bins = linspace(low_bin, high_bin, 50);
        
        air_n = hist(air_data(air_pp), bins);
        wrf_n = hist(wrf_data(wrf_pp), bins);
        
        figure; 
        bar(bins, air_n, 0.8, 'facecolor','none', 'edgecolor','b');
        hold on
        bar(bins, wrf_n, 0.4, 'facecolor', 'none', 'edgecolor','r');
        legend('Aircraft','WRF');
        title(sprintf('%s between %f and %f hPa', chemstr, bottom_pres, top_pres));
    end

    function compare_geos_moz_bcs()
        gc_bc_dir = '/Volumes/share-wrf1/Tests/GeosBC-Test/GeosBCs';
        moz_bc_dir = '/Volumes/share-wrf1/Tests/GeosBC-Test/MozBCs';
        %moz_bc_dir = '/Volumes/share-wrf1/Tests/Cheyenne-Test-w-GeosBCs'; % use this as the "mozart" bc to test if cheyenne's WRF run is close to the yellowstone run
        
        GC = dir(fullfile(gc_bc_dir, 'wrfout*'));
        MOZ = dir(fullfile(moz_bc_dir, 'wrfout*'));
        wrf_sz = get_wrf_array_size(fullfile(gc_bc_dir, GC(1).name));
        
        gci = ncinfo(fullfile(gc_bc_dir, GC(1).name));
        mozi = ncinfo(fullfile(moz_bc_dir, MOZ(1).name));
        
        gc_var = {gci.Variables.Name};
        moz_var = {mozi.Variables.Name};
        vv = find_common_elements(gc_var, moz_var);
        var_dims = zeros(size(vv));
        for a=1:numel(var_dims)
            var_dims(a) = length(gci.Variables(vv(a)).Dimensions);
        end
        wrf_vars = gc_var(vv);
        wrf_vars = wrf_vars(var_dims == 4);
        
        msg = sprintf('What level''s history should we examine? (1-%d)', wrf_sz(3));
        err_msg = sprintf('Choose a single number between 1 and %s', wrf_sz(3));
        lev = ask_number(msg, 'testfxn', @(x) isscalar(x) && x >= 1 && x <= wrf_sz(3), 'testmsg', err_msg);
        
        plot_var = ask_multichoice('Choose a variable to plot:', wrf_vars, 'list', true);
        
        nc_count = [Inf, Inf, 1, Inf];
        nc_start = [1, 1, lev, 1];
        
        gc_dates = date_from_wrf_filenames(GC);
        moz_dates = date_from_wrf_filenames(MOZ);
        
        [gc_xx, moz_xx] = find_common_elements(gc_dates, moz_dates);
        n_times = numel(gc_xx);
        
        gc_vals = nan(wrf_sz(1), wrf_sz(2), n_times);
        moz_vals = nan(wrf_sz(1), wrf_sz(2), n_times);
        for a=1:n_times
            fprintf('Reading file %d of %d\n', a, n_times);
            
            gc_file = fullfile(gc_bc_dir, GC(gc_xx(a)).name);
            moz_file = fullfile(moz_bc_dir, MOZ(moz_xx(a)).name);
            
            if a==1
                xlon = ncread(gc_file, 'XLONG');
                xlat = ncread(gc_file, 'XLAT');
            end
            
            gc_vals(:,:,a) = ncread(gc_file, plot_var, nc_start, nc_count);
            moz_vals(:,:,a) = ncread(moz_file, plot_var, nc_start, nc_count);
            unit_str = ncreadatt(moz_file, plot_var, 'units');
        end
        
        % Make the GUI based plots
        [slon, slat] = state_outlines('not', 'ak', 'hi');
        plot_slice_gui(gc_vals, xlon, xlat, slon, slat, sprintf('%s with GC BCs, level %d (%s)', plot_var, lev, unit_str));
        plot_slice_gui(moz_vals, xlon, xlat, slon, slat, sprintf('%s with MOZ BCs, level %d (%s)', plot_var, lev, unit_str));
        plot_slice_gui(gc_vals - moz_vals, xlon, xlat, slon, slat, sprintf('%s, GC - MOZ, level %d (%s)', plot_var, lev, unit_str));
        plot_slice_gui(reldiff(gc_vals, moz_vals, 'avg')*100, xlon, xlat, slon, slat, sprintf('%s, percent diff GC vs MOZ, level %d', plot_var, lev));
    end

    function make_lightning_met_profs(plotvar)
        allowed_plotvars = {'TT','QVAPOR'};
        if ~exist('plotvar','var')
            plotvar = ask_multichoice('Which variable to plot?', allowed_plotvars, 'list',true);
        else
            if ~ismember(plotvar, allowed_plotvars)
                E.badinput('PLOTVAR must be one of %s', strjoin(allowed_plotvars, ', '));
            end
        end
        nudge_flash_dir = '/Volumes/share2/USERS/LaughnerJ/WRF/DC3/flashcount-nudge';
        nonudge_flash_dir = '/Volumes/share2/USERS/LaughnerJ/WRF/DC3/flashcount-nonudge';
        %nudge_tt_dir = '/Volumes/share2/USERS/LaughnerJ/WRF/DC3/DC3-500_mol_flash-run_to_get_flashrates-nudge';
        nudge_tt_dir = '/Volumes/share2/USERS/LaughnerJ/WRF/DC3/DC3-500_mol_flash-2x_flashrate_nudge';
        %nonudge_tt_dir = '/Volumes/share2/USERS/LaughnerJ/WRF/DC3/DC3-500_mol_flash-run_to_get_flashrates-nonudge';
        nonudge_tt_dir = '/Volumes/share2/USERS/LaughnerJ/WRF/DC3/DC3-500_mol_flash-1x_flashrate_nonudge';
        if strcmpi(plotvar, 'TT')
            save_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/WRF/Lightning Temperature Profs';
        elseif strcmpi(plotvar, 'QVAPOR')
            save_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/WRF/Lightning Qvapor Profs';
        else
            E.notimplemented('No save directory specified for plotvar = %s', plotvar);
        end
        F = dir(fullfile(nudge_tt_dir, 'wrfout*'));
        F2 = dir(fullfile(nonudge_tt_dir, 'wrfout*'));
        F3 = dir(fullfile(nudge_flash_dir, 'wrfout*'));
        F4 = dir(fullfile(nonudge_flash_dir, 'wrfout*'));
        nudge_names = {F.name};
        ff2 = find_common_elements(nudge_names, {F2.name});
        ff3 = find_common_elements(nudge_names, {F3.name});
        ff4 = find_common_elements(nudge_names, {F4.name});
        
        if ~isequal(ff2, ff3)
            E.callError('Different files in %s and %s', nonudge_tt_dir, nudge_flash_dir)
        elseif ~isequal(ff2,ff4)
            E.callError('Different files in %s and %s', nonudge_tt_dir, nonudge_flash_dir)
        end
        
        F = F(ff2);
        wrf_dates = date_from_wrf_filenames(nudge_names(ff2));
        wrf_lon = ncread(fullfile(nudge_flash_dir, F(1).name), 'XLONG');
        wrf_lat = ncread(fullfile(nudge_flash_dir, F(1).name), 'XLAT');
        for a=1:numel(F)
            fprintf('File %d of %d (%s)\n', a, numel(F), F(a).name);
            nudge_ic_flash = ncread(fullfile(nudge_flash_dir, F(a).name), 'IC_FLASHCOUNT');
            nudge_cg_flash = ncread(fullfile(nudge_flash_dir, F(a).name), 'CG_FLASHCOUNT');
            nonudge_ic_flash = ncread(fullfile(nonudge_flash_dir, F(a).name), 'IC_FLASHCOUNT');
            nonudge_cg_flash = ncread(fullfile(nonudge_flash_dir, F(a).name), 'CG_FLASHCOUNT');
            if a == 1
                % flash counts must be differenced, so need at least the
                % second file
                nudge_last_flashes = nudge_ic_flash + nudge_cg_flash;
                nonudge_last_flashes = nonudge_ic_flash + nonudge_cg_flash;
                continue
            end
                
            nudge_flash_map = nudge_ic_flash + nudge_cg_flash - nudge_last_flashes;
            nonudge_flash_map = nonudge_ic_flash + nonudge_cg_flash - nonudge_last_flashes;
            nudge_last_flashes = nudge_ic_flash + nudge_cg_flash;
            nonudge_last_flashes = nonudge_ic_flash + nonudge_cg_flash;
            
            xx = nonudge_flash_map > 0;
            nudge_tt_profs = ncread(fullfile(nudge_tt_dir, F(a).name), plotvar);
            nudge_tt_profs = permute(nudge_tt_profs, [3 1 2]);
            nudge_tt_profs = nudge_tt_profs(:,xx);
            
            nonudge_tt_profs = ncread(fullfile(nonudge_tt_dir, F(a).name), plotvar);
            nonudge_tt_profs = permute(nonudge_tt_profs, [3 1 2]);
            nonudge_tt_profs = nonudge_tt_profs(:,xx);
            
            TT_Profs = struct('date', datestr(wrf_dates(a), 'yyyy-mm-dd HH:MM:SS'),...
                              'nudge_flash_map', nudge_flash_map,...
                              'nonudge_flash_map', nonudge_flash_map,...
                              'map_lons', wrf_lon,...
                              'map_lats', wrf_lat,...
                              'flash_lons', wrf_lon(xx),...
                              'flash_lats', wrf_lat(xx),...
                              'nudge_tt_profs', nudge_tt_profs,...
                              'nonudge_tt_profs', nonudge_tt_profs,...
                              'nudge_flashes', nudge_flash_map(xx),...
                              'nonudge_flashes', nonudge_flash_map(xx)...
                              );
            save_name = sprintf('Lightning_%s_%s.mat', plotvar, datestr(wrf_dates(a), 'yyyy-mm-dd_HH-MM-SS'));
            save(fullfile(save_dir, save_name), 'TT_Profs');
        end
    end

    function fig = plot_lightning_temp_diffs(plotvar)
        allowed_plotvars = {'TT','QVAPOR'};
        if ~exist('plotvar','var')
            plotvar = ask_multichoice('Which variable to plot?', allowed_plotvars, 'list',true);
        else
            if ~ismember(plotvar, allowed_plotvars)
                E.badinput('PLOTVAR must be one of %s', strjoin(allowed_plotvars, ', '));
            end
        end
        if strcmpi(plotvar,'TT')
            data_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/WRF/Lightning Temperature Profs';
            unit_name = 'K';
        elseif strcmpi(plotvar,'QVAPOR')
            data_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/WRF/Lightning Qvapor Profs';
            unit_name = 'kg/kg';
        else
            E.notimplemented('No data path given for plotvar = %s', plotvar);
        end
        F = dir(fullfile(data_dir, 'Lightning*.mat'));
        mean_prof_nudge = nan(29,numel(F));
        mean_prof_nonudge = nan(29,numel(F));
        mean_prof_diffs = nan(29,numel(F));
        mean_prof_reldiffs = nan(29,numel(F));
        weights = nan(29,numel(F));
        
        fprintf('Loading files: ');
        for a=1:numel(F)
            if mod(a,100)==1
                pdone = a/numel(F)*100;
                fprintf('%.1f%% ',pdone);
            end
            D = load(fullfile(data_dir, F(a).name));
            TT_Profs = D.TT_Profs;
            mean_prof_nudge(:,a) = nanmean(TT_Profs.nudge_tt_profs,2);
            mean_prof_nonudge(:,a) = nanmean(TT_Profs.nonudge_tt_profs,2);
            mean_prof_diffs(:,a) = nanmean(TT_Profs.nonudge_tt_profs - TT_Profs.nudge_tt_profs,2);
            mean_prof_reldiffs(:,a) = nanmean(reldiff(TT_Profs.nonudge_tt_profs, TT_Profs.nudge_tt_profs)*100, 2);
            weights(:,a) = size(TT_Profs.nudge_tt_profs,2);
        end
        fprintf('\n');
        
        mean_nudge = nansum2(mean_prof_nudge .* weights, 2) ./ nansum2(weights,2);
        std_nudge = sqrt(var(mean_prof_nudge, weights(1,:), 2));
        mean_nonudge = nansum2(mean_prof_nonudge .* weights, 2) ./ nansum2(weights,2);
        std_nonudge = sqrt(var(mean_prof_nonudge, weights(1,:), 2));
        %mean_diff = nansum2(mean_prof_diffs .* weights,2) ./ nansum2(weights,2);
        %std_diff = sqrt(var(mean_prof_diffs, weights(1,:), 2));
        %mean_reldiff = nansum2(mean_prof_reldiffs .* weights,2) ./ nansum2(weights,2);
        %std_reldiff = sqrt(var(mean_prof_reldiffs, weights(1,:), 2));
        mean_diff = mean_nonudge - mean_nudge;
        std_diff = sqrt(std_nonudge.^2 + std_nudge.^2);
        [mean_reldiff, std_reldiff] = reldiff(mean_nonudge, mean_nudge, std_nonudge, std_nudge);
        mean_reldiff = mean_reldiff * 100;
        std_reldiff = std_reldiff * 100;
        y = (1:numel(mean_diff))';
        
        fig(1) = figure;
        plot(mean_nonudge, y, 'k-', 'linewidth', 2);
        scatter_errorbars(mean_nonudge, y, std_nonudge, 'color', 'k', 'linewidth', 2, 'direction', 'x');
        xlabel(sprintf('%s_{nonudge} (%s)', plotvar, unit_name));
        ylabel('Model level');
        title(sprintf('%s without nudging', plotvar));
        
        fig(2) = figure;
        plot(mean_nudge, y, 'k-', 'linewidth', 2);
        scatter_errorbars(mean_nudge, y, std_nudge, 'color', 'k', 'linewidth', 2, 'direction', 'x');
        xlabel(sprintf('%s_{nudge} (%s)', plotvar, unit_name));
        ylabel('Model level');
        title(sprintf('%s with nudging', plotvar));
        
        fig(2) = figure;
        plot(mean_diff, y, 'k-', 'linewidth', 2);
        scatter_errorbars(mean_diff, y, std_diff, 'color', 'k', 'linewidth', 2, 'direction', 'x');
        xlabel(sprintf('%1$s_{nonudge} - %1$s_{nudge} (%2$s)', plotvar, unit_name));
        ylabel('Model level');
        title(sprintf('%s difference', plotvar));
        
        fig(4) = figure;
        plot(mean_reldiff, y, 'k-', 'linewidth', 2);
        scatter_errorbars(mean_reldiff, y, std_reldiff, 'color', 'k', 'linewidth', 2, 'direction', 'x');
        xlabel(sprintf('(%1$s_{nonudge} - %1$s_{nudge})/%1$s_{nudge} \\times 100%', plotvar));
        ylabel('Model level');
        title(sprintf('%s percent difference', plotvar));
    end


%%%% INTERNAL FUNCTIONS %%%%


    function Match = get_match_struct(Opts)
        if ~exist('Opts','var')
            Opts = struct;
        end
        matchfiles = dir(fullfile(matchdir, '*.mat'));
        matchfiles = {matchfiles.name};
        if ~isfield(Opts,'match_file')
            fname = ask_multichoice('Choose comparison file:', matchfiles, 'list', true);
        else
            if ~ismember(Opts.match_file,matchfiles)
                E.badinput('Match file %s specified in Opts structure is not present in %s', Opts.match_file, matchdir);
            end
            fname = Opts.match_file;
        end
        M = load(fullfile(matchdir, fname));
        matchfield = glob(fieldnames(M), 'Match.*');
        if numel(matchfield)>1
            matchfield = ask_multichoice('Choose match variable', matchfield, 'list', true);
        else
            matchfield = matchfield{1};
        end
        Match = M.(matchfield);
    end

    function [air_data, wrf_data, xstr] = read_data(Match, Opts)
        if ~exist('Opts','var')
            Opts = struct;
        end
        
        comparisons = {'NOx','NO2','MPN','HNO3','NO2:NO'};
        if isfield(Match.wrf,'o3') && isfield('Match.data','o3')
            comparisons{end+1} = 'o3';
        end
        if isfield(Match.wrf,'lnox_total')
            comparisons{end+1} = 'Aircraft NOx vs model LNOx';
        end
        if isfield(Match.wrf, 'PHOTR_NO2')
            comparisons{end+1} = 'NO2 Photolysis';
        end
        if isfield(Match.wrf, 'TT')
            comparisons{end+1} = 'Temperature';
        end
        if isfield(Match.wrf, 'QVAPOR')
            comparisons{end+1} = 'Water vapor mixing ratio';
        end
        if ~isfield(Opts,'quantity')
            quantity = ask_multichoice('Which species to compare?', comparisons, 'list', true);
        else
            if ~ismember(Opts.quantity,comparisons)
                E.badinput('%s is not a valid quantity for this match file; valid quantities are %s', Opts.quantity, strjoin(comparisons, ', '));
            end
            quantity = Opts.quantity;
        end
        % The DC3 NO2 and MPN data have already been corrected for MPN
        % dissociation in the NO2 channel.
        switch lower(quantity)
            case 'nox'
                air_data = (Match.data.no + Match.data.no2);
                wrf_data = (Match.wrf.no + Match.wrf.no2);
                xstr = 'NO_x';
            case 'no2'
                air_data = Match.data.no2;
                wrf_data = Match.wrf.no2;
                xstr = 'NO_2';
            case 'mpn'
                air_data = Match.data.mpn;
                wrf_data = Match.wrf.mpn;
                xstr = 'MPN';
            case 'hno3'
                % HNO3 in Match.data was already taken as the average of
                % the two HNO3 measurements during DC3
                air_data = Match.data.hno3;
                wrf_data = Match.wrf.hno3;
                xstr = 'HNO_3';
            case 'no2:no'
                air_data = Match.data.no2 ./ Match.data.no;
                wrf_data = Match.wrf.no2 ./ Match.wrf.no;
                xstr = 'NO2/NO ratio';
            case 'o3'
                air_data = Match.data.o3;
                wrf_data = Match.wrf.o3;
                xstr = 'O_3';
            case 'aircraft nox vs model lnox'
                air_data = (Match.data.no + Match.data.no2);
                wrf_data = Match.wrf.lnox_total;
                xstr = 'Aircraft NOx, WRF lnox total';
            case 'no2 photolysis'
                air_data = Match.data.PHOTR_NO2;
                wrf_data = Match.wrf.PHOTR_NO2;
                xstr = 'jNO_2 (min^{-1})';
            case 'temperature'
                air_data = Match.data.TT;
                wrf_data = Match.wrf.TT;
                xstr = 'Temperature (K)';
            case 'water vapor mixing ratio'
                air_data = Match.data.QVAPOR;
                wrf_data = Match.wrf.QVAPOR;
                xstr = 'Water vapor mixing ratio (kg/kg)';
        end
        
        % Also make sure that NaNs in one vector are NaNs in the other
        nans = isnan(air_data) | isnan(wrf_data);
        air_data(nans) = nan;
        wrf_data(nans) = nan;
    end

    function figs = compare_flashcounts(nudge_case)
        allowed_cases = {'Long run', 'Nudged DC3'};
        if ~exist('nudge_case', 'var')
            nudge_case = ask_multichoice('Which simulation to use for the nudged case?', allowed_cases, 'list', true);
        else
            if ~ismember(nudge_case, allowed_cases)
                E.badinput('nudge_case must be one of %s', strjoin(allowed_cases, ', '));
            end
        end
        nonudge_a = ncinfo('/Volumes/share2/USERS/LaughnerJ/WRF/DC3/flashcount-nonudge/wrfout_subset_d01_2012-05-13_00-00-00');
        nonudge_b = ncinfo('/Volumes/share2/USERS/LaughnerJ/WRF/DC3/flashcount-nonudge/wrfout_subset_d01_2012-06-24_00-00-00');
        switch lower(nudge_case)
            case 'long run'
                nudge_a = ncinfo('/Volumes/share-wrf1/Outputs/2012/05/wrfout_d01_2012-05-13_00-00-00');
                nudge_b = ncinfo('/Volumes/share-wrf1/Outputs/2012/06/wrfout_d01_2012-06-24_00-00-00');
                nudge_str = 'Long run domain, with nudging';
            case 'nudged dc3'
                nudge_a = ncinfo('/Volumes/share2/USERS/LaughnerJ/WRF/DC3/flashcount-nudge/wrfout_subset_d01_2012-05-13_00-00-00');
                nudge_b = ncinfo('/Volumes/share2/USERS/LaughnerJ/WRF/DC3/flashcount-nudge/wrfout_subset_d01_2012-06-24_00-00-00');
                nudge_str = 'DC3 domain, with nudging';
            otherwise
                E.notimplemented(nudge_case)
        end
        
        nonudge_lon = ncread(nonudge_a.Filename, 'XLONG');
        nonudge_lat = ncread(nonudge_a.Filename, 'XLAT');
        nudge_lon = ncread(nudge_a.Filename, 'XLONG');
        nudge_lat = ncread(nudge_a.Filename, 'XLAT');
        
        nonudge_a_fc = ncread(nonudge_a.Filename, 'IC_FLASHCOUNT') + ncread(nonudge_a.Filename, 'CG_FLASHCOUNT');
        nonudge_b_fc = ncread(nonudge_b.Filename, 'IC_FLASHCOUNT') + ncread(nonudge_b.Filename, 'CG_FLASHCOUNT');
        nonudge_del = nonudge_b_fc - nonudge_a_fc;
        nudge_a_fc = ncread(nudge_a.Filename, 'IC_FLASHCOUNT') + ncread(nudge_a.Filename, 'CG_FLASHCOUNT');
        nudge_b_fc = ncread(nudge_b.Filename, 'IC_FLASHCOUNT') + ncread(nudge_b.Filename, 'CG_FLASHCOUNT');
        nudge_del = nudge_b_fc - nudge_a_fc;
        if ~isequal(nudge_lon, nonudge_lon) || ~isequal(nudge_lat, nonudge_lat)
            F = scatteredInterpolant(double(nudge_lon(:)), double(nudge_lat(:)), double(nudge_del(:)));
            lrdel_interp = F(double(nonudge_lon), double(nonudge_lat));
        else
            lrdel_interp = nudge_del;
        end
        
        figs(1) = figure;
        pcolor(nonudge_lon, nonudge_lat, nonudge_del);
        shading flat
        caxis([0 6000]);
        cb = colorbar;
        cb.Label.String = 'Total flashes';
        set(gca,'fontsize',16)
        title('DC3 domain, no nudging');
        state_outlines('k');
        
        figs(2) = figure;
        pcolor(nudge_lon, nudge_lat, nudge_del);
        shading flat
        caxis([0 6000]);
        cb = colorbar;
        cb.Label.String = 'Total flashes';
        set(gca,'fontsize',16)
        title(nudge_str);
        state_outlines('k');
        
        figs(3) = figure;
        pcolor(nonudge_lon, nonudge_lat, lrdel_interp ./ nonudge_del);
        shading flat
        cb = colorbar;
        cb.Label.String = 'flashcount ratio (nudge / no nudge)';
        set(gca,'fontsize',16)
        title('Ratio');
        state_outlines('k');
        
        xx = nudge_lon > min(nonudge_lon(:)) & nudge_lon < max(nonudge_lon(:));
        yy = nudge_lat > min(nonudge_lat(:)) & nudge_lat < max(nonudge_lat(:));
        
        boxvec = cat(1, nonudge_del(:), nudge_del(xx & yy));
        gvec = cat(1,...
                    repmat({'No nudging'},numel(nonudge_del),1),...
                    repmat({'With nudging'},sum(xx(:)&yy(:)),1));
        figs(4) = figure; 
        boxplot(boxvec, gvec);
        set(gca,'fontsize',16);
    end
end

function lats_xx = subset_latitudes(Match, Opts)
E = JLLErrors;
if ~exist('Opts', 'var')
    Opts = struct; 
end
allowed_vals = {'All','Midlatitude (>35 N)','Subtropics (20-35 N)'};
short_vals = {'all','midlats','subtrop'};
if ~isfield(Opts,'lats_filter')
    lats_choice = ask_multichoice('Which latitudes to use?', allowed_vals, 'list', true);
else
    xx = strcmpi(Opts.lats_filter, short_vals);
    if sum(xx) ~= 1
        E.badinput('Opts.lats_filter must be one of %s (case insensitive)', strjoin(short_vals, ', '));
    end
    lats_choice = allowed_vals{xx};
end
switch lats_choice
    case 'All'
        lats_xx = true(size(Match.data.lat));
    case 'Midlatitude (>35 N)'
        lats_xx = Match.data.lat > 35;
    case 'Subtropics (20-35 N)'
        lats_xx = Match.data.lat > 20 & Match.data.lat <= 35;
    otherwise
        E.notimplemented(lats_choice);
end
end

function trop_xx = filter_strat_influence(Match, Opts)
E = JLLErrors;
if ~exist('Opts', 'var')
    Opts = struct; 
end
if ~isfield(Match.data, 'o3') || ~isfield(Match.data, 'co')
    fprintf('O3 and CO from DC3 not present in Match file, cannot filter stratospheric influence\n');
    trop_xx = true(size(Match.data.lon));
    return
end

if ~isfield(Opts,'strat_filter')
    strat_choice = strcmpi(ask_multichoice('Remove stratospheric influence (O3/CO > 1.5)?', {'y','n'}),'y');
else
    if ~isscalar(Opts.strat_filter) || ~islogical(Opts.strat_filter)
        E.badinput('Opts.strat_filter must be a scalar logical, if given')
    end
    strat_choice = Opts.strat_filter;
end
if strat_choice
    trop_xx = Match.data.o3 ./ Match.data.co < 1.5;
else
    trop_xx = true(size(Match.data.o3));
end
end

function fresh_xx = filter_fresh_lnox(Match,Opts)
E = JLLErrors;
if ~exist('Opts', 'var')
    Opts = struct; 
end
if ~isfield(Opts,'fresh_filter')
    fresh_choice = strcmpi(ask_multichoice('Remove fresh lightning NOx/freshly emitted NOx (NOx/HNO3 > 5)?',{'y','n'}),'y');
else
    if ~isscalar(Opts.fresh_filter) || ~islogical(Opts.fresh_filter)
        E.badinput('Opts.fresh_filter must be a scalar logical, if given')
    end
    fresh_choice = Opts.fresh_filter;
end
if fresh_choice
    fresh_xx = (Match.data.no + Match.data.no2) ./ Match.data.hno3 < 5;
else
    fresh_xx = true(size(Match.data.no));
end
end

function uu = unique_wrf_element(Match)
% Returns xx, a vector of indicies which will, when used to index any of
% the Match.wrf data fields, return only the first instance of the WRF
% value for a particular grid cell and time.

% This can be nicely handled by unique() as soon as we make a matrix where
% the rows are the indicies for each data element.

i_we = Match.indicies.west_east(:);
i_sn = Match.indicies.south_north(:);
i_bt = Match.indicies.bottom_top(:);
i_time = Match.indicies.time(:);

i_mat = [i_we, i_sn, i_bt, i_time];
nans = find(any(isnan(i_mat),2));

[~, xx] = unique([i_we, i_sn, i_bt, i_time], 'rows', 'stable');
xx(ismember(xx,nans)) = [];
uu = false(1, length(i_we));
uu(xx) = true;
end

