function [ varargout ] = misc_amf_lnox_plots( plttype, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E=JLLErrors;
homedir = getenv('HOME');
workspace_dir = fullfile(homedir, 'Documents','MATLAB','BEHR','Workspaces','LNOx-AMFs');
share_dir = fullfile('/Volumes', 'share2', 'USERS', 'LaughnerJ', 'WRF', 'DC3');

switch lower(plttype)
    case 'factor-hists'
        figout = make_factor_histograms(varargin{:});
        if nargout > 0
            varargout{1} = figout;
        end
    case 'contours'
        [figout1, figout2] = make_amf_contour_plots(varargin{:});
        if nargout > 0
            varargout = {figout1, figout2};
        end
    case '1factor-series'
        make_1_factor_series()
    case 'dc3-comp-table'
        make_dc3_comparison_table();
    case 'no2-profiles'
        figout = make_dc3_matched_profiles('NO2');
        if nargout > 0
            varargout{1} = figout;
        end
    case 'nox-profiles'
        make_dc3_matched_profiles('NOx');
    case 'allwrf-no2-profiles'
        varargout{1} = make_avg_wrf_profiles();
    case 'behr-avgs'
        do_behr_avgs();
    case 'plot-behr-avgs';
        figout = plot_behr_avgs(varargin{:});
        if nargout > 0
            varargout{1} = figout;
        end
    case 'sw-plots'
        figout = make_sw_plots(varargin{:});
        if nargout > 0
            varargout{1} = figout;
        end
    otherwise
        E.badinput('PLTTYPE %s not recognized', plttype)
end

    function fig = make_factor_histograms(amf_file)
        if ~exist('amf_file', 'var')
            F = dir(fullfile(workspace_dir, '*.mat'));
            fnames = {F.name};
            amf_file = ask_multichoice('Which AMF file should be used?', fnames, 'list', true);
        end
        D = load(fullfile(workspace_dir, amf_file));
        fns = sort(fieldnames(D));
        factors = {'SZA','VZA','RAA','Alb','Surf P'};
        means = nan(numel(fns), numel(factors));
        stddevs = nan(numel(fns), numel(factors));
        
        x = 1:numel(factors);
        markers = {'o','s','^'};
        cols = {'k','b','r'};
        l = gobjects(3,1);
        fig = figure;
        for a=1:numel(fns)
            offset = (a-2)*0.15;
            for b=1:numel(factors)
                spreads = range(D.(fns{a}).AMFs,b);
                means(a,b) = mean(spreads(:));
                stddevs(a,b) = std(spreads(:));
            end
            l(a) = line(x + offset, means(a,:), 'marker', markers{a}, 'color', cols{a}, 'linestyle', 'none', 'linewidth', 2, 'markersize', 12);
            scatter_errorbars(x+offset, means(a,:), stddevs(a,:), 'color', cols{a}, 'linewidth',2);
        end
        set(gca,'xtick',x,'xticklabel',factors)
        set(gca,'fontsize', 16)
        
        legstr = cell(1,numel(fns));
        for a=1:numel(fns)
            mol_flash = str2double(regexp(fns{a}, '\d\d\d', 'match', 'once'));
            if mol_flash == 0
                legstr{a} = 'No Lightning';
            else
                legstr{a} = sprintf('%d mol flash^{-1}', mol_flash);
            end
        end
        
        legend(l, legstr)
        ylabel('AMF Range')
    end

    function [fig1, fig2] = make_amf_contour_plots(amf_file, clear_or_cloudy, difference)
        % make contour plots of the difference in AMFs with and without
        % lightning.
        if ~exist('amf_file', 'var')
            F = dir(fullfile(workspace_dir, '*.mat'));
            fnames = {F.name};
            amf_file = ask_multichoice('Which AMF file should be used?', fnames, 'list', true);
        end
        if ~exist('clear_or_cloudy', 'var')
            try_reg = regexpi(amf_file, 'clear|cloudy', 'once', 'match');
            if ~isempty(try_reg)
                clear_or_cloudy = lower(try_reg);
            else
                clear_or_cloudy = 'clear';
            end
        else
            if ~ismember(lower(clear_or_cloudy), {'clear','cloudy'})
                E.badinput('If given, CLEAR_OR_CLOUDY must be ''clear'' or ''cloudy''')
            end
        end
        allowed_diffs = {'500-0', '665-500'};
        if ~exist('difference', 'var')
            difference = ask_multichoice('Which difference to plot?', allowed_diffs, 'list', true);
        else
            if ~ismember(difference, allowed_diffs)
                E.badinput('DIFFERENCE type "%s" not valid (permitted values: %s)', difference, strjoin(allowed_diffs));
            end
        end
        
        D = load(fullfile(workspace_dir, amf_file));
        
        [sza_vec, vza_vec, ~, alb_vec, surfp_vec] = param_vecs(D.amfs_wrf_000lnox);
        [def_sza, def_vza, def_alb, def_surfp] = lut_def_params(clear_or_cloudy, alb_vec);
        
        % First make plots for clear sky conditions, SZA and VZA first,
        % then ALB and SurfP
        if strcmpi(difference, '500-0')
            perdiff = reldiff(D.amfs_wrf_500lnox.AMFs, D.amfs_wrf_000lnox.AMFs)*100;
        elseif strcmpi(difference, '665-500')
            perdiff = reldiff(D.amfs_wrf_665lnox.AMFs, D.amfs_wrf_500lnox.AMFs)*100;
        else
            E.notimplemented('Difference type "%s" not implemented', difference);
        end
        
        meandiff = squeeze(mean(perdiff,3));
        subdiff = meandiff(:,:,def_alb,def_surfp);
        tstr = sprintf('Albedo = %.2f, Surf. P = %.1f', alb_vec(def_alb), surfp_vec(def_surfp));
        fig1 = do_contour(sza_vec, vza_vec, subdiff, 'SZA', 'VZA', tstr);
        
        subdiff = squeeze(meandiff(def_sza, def_vza, :, :));
        tstr = sprintf('SZA = %.1f, VZA = %.1f', sza_vec(def_sza), vza_vec(def_vza));
        fig2 = do_contour(alb_vec, surfp_vec, subdiff, 'Albedo', 'Surface Pressure', tstr);
        set(gca,'ydir','reverse');
    end

    function fig = make_dc3_matched_profiles(quantity)
        allowed_quantities={'NOx','NO2'};
        if ~ismember(quantity, allowed_quantities)
            E.badinput('QUANTITY must be one of %s', strjoin(allowed_quantities, ', '));
        end
        Opts.quantity = quantity;
        Opts.lats_filter = 'all';
        Opts.strat_filter = false;
        Opts.fresh_filter = false;
        
        Profs(1).file = 'DC3-Comparison-no_lnox-fixedBC-using_wrfgridcorners.mat';
        Profs(1).linestyle = '-.';
        Profs(1).legend = 'WRF - No LNO_x';
        Profs(1).color = 'r';
        Profs(2).file = 'DC3-Comparison-flashrate-1-molflash-500-iccg-2-newprofile-fixedBC-using_wrfgridcorners.mat';
        Profs(2).linestyle = '-';
        Profs(2).legend = 'WRF - 500 mol LNO_x/flash';
        Profs(2).color = 'r';
        Profs(3).file = 'DC3-Comparison-flashrate-1-molflash-665-iccg-2-newprofile-fixedBC-using_wrfgridcorners.mat';
        Profs(3).linestyle = '--';
        Profs(3).legend = 'WRF - 665 mol LNO_x/flash';
        Profs(3).color = 'r';
        Profs(4).file = 'DC3-Comparison-500mol-nudge-test-wTemperature.mat';
        Profs(4).legend = 'WRF - 500 mol LNOx/fl - nudged';
        Profs(4).linestyle = '--';
        Profs(4).color = 'k';
        Profs(5).file = 'DC3-Comparison-500mol-nudge-2xflashrate-wTTandQ.mat';
        Profs(5).legend = 'WRF - 500 mol/fl, 2x flashrate - nudged';
        Profs(5).linestyle = ':';
        Profs(5).color = 'k';
        
        n = numel(Profs);
        l = gobjects(n+1,1);
        fig=figure;
        for a=1:n
            Opts.match_file = Profs(a).file;
            [Profs(a).dc3_med,Profs(a).dc3_mean,Profs(a).dc3_pres,Profs(a).wrf_med,Profs(a).wrf_mean,Profs(a).wrf_pres] = misc_wrf_chem_comp_plots('grand-profile',Opts);
            l(a) = line(Profs(a).wrf_mean*1e6, Profs(a).wrf_pres, 'color', Profs(a).color, 'linewidth', 2, 'linestyle', Profs(a).linestyle);
        end
        l(n+1) = line(Profs(1).dc3_mean*1e6, Profs(a).dc3_pres, 'color', 'b', 'linewidth', 2, 'linestyle', '-');
        set(gca,'fontsize',14,'ydir','reverse');
        xlabel(sprintf('[%s] (pptv)',quantity));
        ylabel('Pressure (hPa)');
        lstr = {Profs.legend, 'DC3'};
        legend(l,lstr);
    end

    function profs = make_avg_wrf_profiles()
        wrf_dirs = {fullfile(share_dir,'lnox_off-fixed_BCs'),...
                    fullfile(share_dir,'iccg_eq_2-fr_factor_1-mol_flash_500_newprof-fixedBC'),...
                    fullfile(share_dir,'iccg_eq_2-fr_factor_1-mol_flash_665-fixedBC')};
        prof_ids = {'lightning_000mol', 'lightning_500mol', 'lightning_665mol'};
        profs = make_empty_struct_from_cell(prof_ids);
        start_date = '2012-05-18';
        
        for a=1:numel(wrf_dirs)
            F = dir(fullfile(wrf_dirs{a}, 'wrfout_subset*'));
            wrf_dnums = date_from_wrf_filenames(F);
            xx = wrf_dnums >= datenum(start_date) & hour(wrf_dnums) >= 17 & hour(wrf_dnums) <= 22 & minute(wrf_dnums) == 0;
            F = F(xx);
            
            [no2, pres] = read_wrf_vars(wrf_dirs{a}, F, {'no2', 'pres'});
            [profs.(prof_ids{a}), profs.pressures] = bin_omisp_pressure(pres, no2, 'mean');
%             pvec = 1:ndims(no2);
%             pvec(3) = [];
%             pvec = [3, pvec];
%             no2 = permute(no2, pvec);
%             pres = permute(pres, pvec);
%             sz = size(no2);
%             n = prod(sz(2:end));
%             pressures = behr_pres_levels();
%             pressures = pressures(:); % ensure is column vector
%             no2_interp = nan(numel(pressures), n);
%             for b=1:n
%                 no2_tmp = interp1(log(pres(:,b)), log(no2(:,b)), log(pressures));
%                 no2_interp(:,b) = exp(no2_tmp);
%             end
%             
%             profs.(prof_ids{a}) = nanmean(no2_interp,2);
        end
    end

    function fig = make_sw_plots(amf_file)
        if ~exist('amf_file', 'var')
            F = dir(fullfile(workspace_dir, '*.mat'));
            fnames = {F.name};
            amf_file = ask_multichoice('Which AMF file should be used?', fnames, 'list', true);
            
        end
        
        D = load(fullfile(workspace_dir, amf_file));
        [sza_vec, vza_vec, raa_vec, alb_vec, surfp_vec] = param_vecs(D.amfs_wrf_000lnox);
        [def_sza, def_vza, def_alb, def_surfp] = lut_def_params(amf_file, alb_vec);
        
        xlstr = 'Normalized scattering weights';
        ylstr = 'Pressure (hPa)';
        
        fig=figure;
        fig.Position(3:4) = fig.Position(3:4)*2;
        subplot(2,2,1);
        make_indiv_sw_plot(sza_vec, vza_vec(def_vza), raa_vec(end), alb_vec(def_alb), surfp_vec(def_surfp));
        label_axis_with_letter('(a)', 'xshift', 0.15);
        ylabel(ylstr);
        subplot(2,2,2);
        make_indiv_sw_plot(sza_vec(def_sza), vza_vec, raa_vec(end), alb_vec(def_alb), surfp_vec(def_surfp));
        label_axis_with_letter('(b)', 'xshift', 0.15);
        subplot(2,2,3);
        make_indiv_sw_plot(sza_vec(def_sza), vza_vec(def_vza), raa_vec(end), alb_vec, surfp_vec(def_surfp));
        label_axis_with_letter('(c)', 'xshift', 0.15);
        xlabel(xlstr); ylabel(ylstr);
        subplot(2,2,4);
        make_indiv_sw_plot(sza_vec(def_sza), vza_vec(def_vza), raa_vec(end), alb_vec(def_alb), surfp_vec);
        label_axis_with_letter('(d)', 'xshift', 0.15);
        xlabel(xlstr);
        
        
    end

    function make_indiv_sw_plot(sza, vza, raa, alb, surfp)
        n = [numel(sza), numel(vza), numel(raa), numel(alb), numel(surfp)];
        maxn = max(n);
        if sum(n>1) > 1
            E.badinput('Must input a vector for only one parameter');
        end
        
        params = {sza, vza, raa, alb, surfp};
        param_names = {'SZA', 'VZA', 'RAA', 'ALB', 'SurfP'};
        titlestr = '';
        for a=1:numel(params)
            if n(a) == 1
                if ~isempty(titlestr)
                    sep = ', ';
                else
                    sep = '';
                end
                titlestr = [titlestr, sprintf('%s%s = %.2f', sep, param_names{a}, params{a})];
                params{a} = repmat(params{a},maxn,1);
            else
                params{a} = params{a}(:);
                varied_name = param_names{a};
                varied_ind = a;
            end
        end
        
        
        fileDamf = fullfile(BEHR_paths('amf_tools_dir'),'damf.txt');
        plevs = behr_pres_levels;
        lobjs = gobjects(maxn,1);
        lstr = cell(1,maxn);
        plims = [Inf, -Inf];
        for a=1:maxn
            sc_wts = rDamf2(fileDamf, plevs, params{1}(a), params{2}(a), params{3}(a), params{4}(a), params{5}(a));
            color = map2colmap(a,1,maxn,'jet');
            
            pp = plevs < params{5}(a);
            %plims = [min(plims(1), min(plevs(pp))), max(plims(2), max(plevs(pp)))];
            lobjs(a) = line(sc_wts(pp) ./ sc_wts(end), plevs(pp), 'linewidth', 2, 'color', color);
            if strcmpi(varied_name, 'ALB')
                lstr{a} = sprintf('%s = %.2g', varied_name, params{varied_ind}(a));
            else
                lstr{a} = sprintf('%s = %.1f', varied_name, params{varied_ind}(a));
            end
        end
        legend(lobjs, lstr, 'location', 'northwest');
        title(titlestr);
        set(gca,'fontsize',12,'ydir','reverse');
        ylim([min(plevs), max(plevs)]);
    end

    function do_behr_avgs()
        work_dir_root = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/LNOx-AMFs';
        subdirs = {'BEHR_noLNOx', 'BEHR_500LNOx', 'BEHR_665LNOx'};
        cloud_types = {'omi', 'modis'};
        %data_fields = {'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly'};
        data_fields = {'BEHRAMFTrop', 'BEHRAMFTropVisOnly'};
        
        lonlim = [-108, -87];
        latlim = [27 43];
        start_date = '2012-05-18';
        end_date = '2012-06-23';
        
        for a=1:numel(subdirs)
            behrdir = fullfile(work_dir_root, subdirs{a});
            for b=1:numel(cloud_types)
                for c=1:numel(data_fields)
                    fprintf('***** Averaging %s in %s using %s clouds *****\n', data_fields{c}, subdirs{a}, cloud_types{b});
                    [~, NO2_GRID, LON_GRID, LAT_GRID, COUNT, opts] = no2_column_map_2014(start_date, end_date, lonlim, latlim,...
                        'behrdir', behrdir, 'fileprefix', 'OMI_BEHR_v2-1C_', 'clouds', cloud_types{b},...
                        'cloudfraccrit', 0.2, 'mapfield', data_fields{c}, 'makefig', false); %#ok<ASGLU>
                    savename = sprintf('BEHR_avg_%s-%s_cloud.mat', data_fields{c}, cloud_types{b});
                    save(fullfile(behrdir, savename), 'NO2_GRID', 'LON_GRID', 'LAT_GRID', 'COUNT', 'opts');
                end
            end
        end
    end

    function figs = plot_behr_avgs(plot_var, cloud_type)
        allowed_plot_vars = {'BEHRColumnAmountNO2Trop', 'BEHRColumnAmountNO2TropVisOnly', 'BEHRAMFTrop', 'BEHRAMFTropVisOnly'};
        if ~exist('plot_var', 'var')
            plot_var = ask_multichoice('Use averages of which variable?', allowed_plot_vars, 'list', true);
        else
            if ~ismember(plot_var, allowed_plot_vars)
                E.badinput('PLOT_VAR must be one of %s', strjoin(allowed_plot_vars, ', '));
            end
        end
        
        allowed_cloud_types = {'omi', 'modis'};
        if ~exist('cloud_type', 'var')
            cloud_type = ask_multichoice('Use averages filtered by which cloud type?', allowed_cloud_types, 'list', true);
        else
            if ~ismember(cloud_type, allowed_cloud_types)
                E.badinput('CLOUD_TYPE must be one of %s', strjoin(allowed_cloud_types, ', '));
            end
        end
        
        fname = sprintf('BEHR_avg_%s-%s_cloud.mat', plot_var, cloud_type);
        lnox_000 = load(fullfile(workspace_dir, 'BEHR_noLNOx', fname));
        lnox_500 = load(fullfile(workspace_dir, 'BEHR_500LNOx', fname));
        lnox_665 = load(fullfile(workspace_dir, 'BEHR_665LNOx', fname));
        
        % Pieces of the label - which variable and if it is visible only or
        % not
        var_name = regexp(plot_var, '(ColumnAmount|AMF)', 'match', 'once');
        switch lower(var_name)
            case 'columnamount'
                var_name = 'VCD';
        end
        vis_or_total = regexp(plot_var, 'VisOnly', 'once');
        if isempty(vis_or_total)
            vis_or_total = '';
        else
            vis_or_total = 'Visible ';
        end
        
        dels = {reldiff(lnox_500.NO2_GRID, lnox_000.NO2_GRID)*100,...
                reldiff(lnox_665.NO2_GRID, lnox_500.NO2_GRID)*100,...
                lnox_500.NO2_GRID - lnox_000.NO2_GRID,...
                lnox_665.NO2_GRID - lnox_500.NO2_GRID};
        del_names = {'%%\\Delta %s%s 500 vs 0 mol flash^{-1}',...
                     '%%\\Delta %s%s 665 vs 500 mol flash^{-1}',...
                     '\\Delta %s%s 500 vs 0 mol flash^{-1}',...
                     '\\Delta %s%s 665 vs 500 mol flash^{-1}'};
        del_names = cprintf(del_names, vis_or_total, var_name);
        lon = lnox_000.LON_GRID;
        lat = lnox_000.LAT_GRID;
        
        figs = gobjects(numel(dels),1);
        
        for a=1:numel(dels)
            figs(a) = figure;
            pcolor(lon, lat, dels{a});
            shading flat;
            cb = colorbar;
            cb.Label.String = del_names{a};
            caxis([-max(abs(cb.Ticks)), max(abs(cb.Ticks))])
            colormap(blue_red_cmap);
            set(gca,'fontsize',16);
            state_outlines('k');
        end
    end
end

function [sza_vec, vza_vec, raa_vec, alb_vec, surfp_vec] = param_vecs(amfs)
sza_vec = reshape(amfs.SZAs(:,1,1,1,1),[],1);
vza_vec = reshape(amfs.VZAs(1,:,1,1,1),[],1);
raa_vec = reshape(amfs.RAAs(1,1,:,1,1),[],1);
alb_vec = reshape(amfs.ALBs(1,1,1,:,1),[],1);
surfp_vec = reshape(amfs.SurfPs(1,1,1,1,:),[],1);
end

function [def_sza, def_vza, def_alb, def_surfp] = lut_def_params(clear_or_cloudy, alb_vec)
E = JLLErrors;
try_reg = regexpi(clear_or_cloudy, 'clear|cloudy', 'once', 'match');
if ~isempty(try_reg)
    clear_or_cloudy = lower(try_reg);
else
    E.badinput('Cannot tell if clear or cloudy');
end

def_sza = 5;
def_vza = 5;
if strcmpi(clear_or_cloudy, 'clear')
    [~,def_alb] = min(abs(alb_vec - 0.04));
    def_surfp = 1;
elseif strcmpi(clear_or_cloudy, 'cloudy')
    [~,def_alb] = min(abs(alb_vec - 0.8));
    def_surfp = 5;
end
end

function fig=do_contour(xvec,yvec,subdiff,xstr,ystr,tstr)
[Y,X] = meshgrid(yvec, xvec);
fig=figure;
contourf(X, Y, subdiff);
xlabel(xstr); ylabel(ystr);
cb=colorbar; cb.Label.String = '%\Delta AMF 500 mol flash^{-1} - No LNO_x';
title(tstr);
set(gca,'fontsize',16)
end