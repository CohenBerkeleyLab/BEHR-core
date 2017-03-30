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
            varargout{4} = Match.wrf.pres(data_log);
            varargout{5} = xstr;
        end
    end

    function [air_data_bins, air_data_means, air_pres_bins, wrf_data_bins, wrf_data_means, wrf_pres_bins] = pres_binned_grand_profile(Opts)
        if ~exist('Opts','var')
            Opts = struct;
        end
        [air_data, air_pres, wrf_data, wrf_pres, xstr] = grand_scatterplot(Opts);
        [air_data_bins, air_pres_bins, air_quants] = bin_omisp_pressure(air_pres, air_data);
        [wrf_data_bins, wrf_pres_bins, wrf_quants] = bin_omisp_pressure(wrf_pres, wrf_data);
        air_data_means = bin_omisp_pressure(air_pres, air_data, 'mean');
        wrf_data_means = bin_omisp_pressure(wrf_pres, wrf_data, 'mean');
        
        if nargout == 0
            figure;
            l=gobjects(4,1);
            l(1) = line(air_data_bins, air_pres_bins, 'color', 'b', 'linestyle', '-', 'marker', 'o');
            l(2) = line(wrf_data_bins, wrf_pres_bins, 'color', 'r', 'linestyle', '-', 'marker', 'v');
            l(3) = line(air_data_means, air_pres_bins, 'color', 'b', 'linestyle', '--', 'marker', 's');
            l(4) = line(wrf_data_means, wrf_pres_bins, 'color', 'r', 'linestyle', '--', 'marker', '^');

            legend(l, {'Aircraft Data Med.', 'WRF Med.', 'Aircraft Mean','WRF Mean'});
            set(gca,'fontsize',16);
            xlabel(xstr);
            ylabel('Pressure (hPa)');
            set(gca,'ydir','reverse');
        end
    end

    function varargout = amf_test(req_var_out)
        % For a given run, this will compute the AMF for four cases:
        %  1) DC3 profile
        %  2) WRF profile
        %  3) DC3 boundary layer, WRF free trop
        %  4) WRF boundary layer, DC3 free trop
        [ ~, dc3_prof, dc3_pres, ~, wrf_prof, wrf_pres] = pres_binned_grand_profile();
        if any(dc3_pres ~= wrf_pres)
            E.callError('different_pres','DC3 and WRF profiles on different pressure levels');
        end
        dc3_ft_wrf_bl = dc3_prof;
        dc3_ft_wrf_bl(dc3_pres > 750) = wrf_prof(wrf_pres > 750);
        wrf_ft_dc3_bl = wrf_prof;
        wrf_ft_dc3_bl(wrf_pres > 750) = dc3_prof(dc3_pres > 750);
        dc3_mt_wrf_bl_ut = wrf_prof;
        dc3_mt_wrf_bl_ut(wrf_pres < 750 & wrf_pres > 375) = dc3_prof(dc3_pres < 750 & dc3_pres > 375);
        dc3_rl_wrf_bl_ft = wrf_prof;
        dc3_rl_wrf_bl_ft(wrf_pres >= 700 & wrf_pres <= 875) = dc3_prof(dc3_pres >= 700 & dc3_pres <= 875);
        
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
        end
        
        % Also make sure that NaNs in one vector are NaNs in the other
        nans = isnan(air_data) | isnan(wrf_data);
        air_data(nans) = nan;
        wrf_data(nans) = nan;
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

