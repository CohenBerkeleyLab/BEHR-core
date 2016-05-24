function [  ] = fit_var_plotting(plottype, varargin )
%FIT_VAR_PLOTTING(PLOTTYPE, FFIXED, F, ...) Misc. plotting functions for output of fit_line_density_var.m
%   Collection of miscellaneous plotting functions that takes the Ffixed
%   and F or Nfixed and N outputs from fit_line_density_variation.m to show
%   how the goodness of fit or the other fitting parameters change when one
%   fitting parameter is fixed. Specify the plot desired using the PLOTTYPE
%   input string. Options are:
%
%       'residuals' - plot the change in residuals for each value of the
%       fixed parameters, relative to their value at the identified
%       minimum. Takes Ffixed and F (or Nfixed and N) as the first two
%       inputs and the wind category (as a string) for the third. The wind
%       category will be used in the title.
%
%       'crosseffects' - plot the change in each of the other four fitting
%       parameters when one is fixed. Takes the same inputs as 'residuals'
%
%       'show_fits' - will plot the actual fits against the real data for
%       each value the fixed parameter takes on. Takes the same first three
%       arguments as the 'residuals' plot, but you must also specify which
%       fixed parameter you want to plot for. This must be a string that
%       matches the string in Ffixed.fixed_par.

E=JLLErrors;

first_warning = true;

if ~exist('plottype','var')
    plottype = ask_multichoice('Which plot to make?', {'fits'},'list',true);
end

switch lower(plottype)
    case 'fits'
        plot_3_fits()
    case 'residuals'
        residuals(varargin{:})
    case 'residuals-onepar'
        residuals_onepar(varargin{:})
    case 'crosseffects'
        cross_effects(varargin{:})
    case 'showfits'
        show_fits(varargin{:})
    case 'mchist'
        mchist(varargin{:});
    case 'stats'
        stats_table();
    otherwise
        fprintf('Plot type not recognized\n')
end

    function plot_3_fits
        [data_file, data_path] = uigetfile('*SimpleFits*.mat','Select the fits file for the conditions to plot');
        if data_file == 0
            fprintf('Cancelled\n')
            return
        end
        SF = load(fullfile(data_path,data_file));
        ldfile = regexprep(data_file,'SimpleFits(-ssresid|-unexvar)*','LineDensities');
        if ~exist(fullfile(data_path,ldfile),'file')
            E.filenotfound('simple fits file');
        else
            LD = load(fullfile(data_path,ldfile));
        end
        
        % Get city name
        [s,e] = regexp(data_file,'(Atlanta|Birmingham|Montgomery)');
        city_name = data_file(s:e);
        
        % Get wind conditions
        WC = load(sprintf('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/%s-Wind-Conditions-1900UTC-5layers.mat',city_name));
        
        % Determine whether to plot the WRF ones too
        plot_wrf = false;
        if all(isfield(LD,{'no2x_wrffast','no2x_wrfslow'}))
            plot_wrf_ans = questdlg('Include WRF line densities?');
            if strcmpi(plot_wrf_ans,'Yes')
                plot_wrf = true;
            elseif strcmpi(plot_wrf_ans,'Cancel')
                return
            end
        end
        
        figure; hold on
        lstr = {'Coarse monthly data','Coarse monthly fit', 'Fine monthly data', 'Fine monthly fit', 'Hybrid data', 'Hybrid fit'};
        plot(LD.no2x_mn108slow, LD.no2ld_mn108slow, 'go', 'linewidth', 2);
        plot(LD.no2x_mn108slow, SF.f_mn108slow.emgfit, '--', 'color', [0 0.5 0], 'linewidth', 2);
        plot(LD.no2x_mnslow, LD.no2ld_mnslow, 'o', 'linewidth', 2, 'color', [1 0.5 0]);
        plot(LD.no2x_mnslow, SF.f_mnslow.emgfit, 'r--', 'linewidth', 2);
        plot(LD.no2x_hyslow, LD.no2ld_hyslow, 'co', 'linewidth', 2);
        plot(LD.no2x_hyslow, SF.f_hyslow.emgfit, 'b--', 'linewidth', 2);
        if plot_wrf
            plot(LD.no2x_wrfslow, LD.no2ld_wrfslow, 'ko', 'linewidth',2);
            plot(LD.no2x_wrfslow, SF.f_wrfslow.emgfit, 'k--','linewidth',2)
            lstr = [lstr, {'WRF data','WRF fit'}];
        end
        legend(lstr{:});
        xlabel('Distance from city (km)')
        ylabel('Line density (mol km^{-1})')
        set(gca,'fontsize',20)
        if isfield(LD,'wind_crit')
            title(sprintf('%s: wind < %.1f',city_name,LD.wind_crit))
        else
            title(sprintf('%s: slow wind',city_name))
        end
        
        figure; hold on
        plot(LD.no2x_mn108fast, LD.no2ld_mn108fast, 'go', 'linewidth', 2);
        plot(LD.no2x_mn108fast, SF.f_mn108fast.emgfit, '--', 'color', [0 0.5 0], 'linewidth', 2);
        plot(LD.no2x_mnfast, LD.no2ld_mnfast, 'o', 'linewidth', 2, 'color', [1 0.5 0]);
        plot(LD.no2x_mnfast, SF.f_mnfast.emgfit, 'r--', 'linewidth', 2);
        plot(LD.no2x_hyfast, LD.no2ld_hyfast, 'co', 'linewidth', 2);
        plot(LD.no2x_hyfast, SF.f_hyfast.emgfit, 'b--', 'linewidth', 2);
        if plot_wrf
            plot(LD.no2x_wrffast, LD.no2ld_wrffast, 'ko', 'linewidth',2);
            plot(LD.no2x_wrffast, SF.f_wrffast.emgfit, 'k--','linewidth',2)
            %lstr = [lstr, {'WRF data','WRF fit'}]; already done in slow
        end
        legend(lstr{:});
        xlabel('Distance from city (km)')
        ylabel('Line density (mol km^{-1})')
        set(gca,'fontsize',20)
        if isfield(LD,'wind_crit')
            title(sprintf('%s: wind \\geq %.1f',city_name,LD.wind_crit))
        else
            title(sprintf('%s: fast wind',city_name))
        end
        
        % Make the fitting parameters into a matrix (to print in latex
        % format) and a table (to look at visually)
        fns = {'f_mn108fast','f_mnfast','f_hyfast','f_mn108slow','f_mnslow','f_hyslow'};
        if plot_wrf
            E.notimplemented('plot_wrf needs updated with uncertainty')
            A = cat(1, struct2array(SF.f_mn108fast.ffit), struct2array(SF.f_mnfast.ffit), struct2array(SF.f_hyfast.ffit), struct2array(SF.f_wrffast.ffit),...
                struct2array(SF.f_mn108slow.ffit), struct2array(SF.f_mnslow.ffit), struct2array(SF.f_hyslow.ffit), struct2array(SF.f_wrfslow.ffit))';
            varnames = {'Mn108Fast','MnFast','HyFast','WRFFast','Mn108Slow','MnSlow','HySlow','WRFSlow'};
        else
            %Aold = cat(1, struct2array(SF.f_mn108fast.ffit), struct2array(SF.f_mnfast.ffit), struct2array(SF.f_hyfast.ffit),...
            %    struct2array(SF.f_mn108slow.ffit), struct2array(SF.f_mnslow.ffit), struct2array(SF.f_hyslow.ffit))';
            A = [];
            
            for a=1:numel(fns)
                vals = struct2array(SF.(fns{a}).ffit);
                apriori = fns{a}(3:end); % remove the "f_" bit
                uncertainty = calc_total_fit_uncert(SF.(fns{a}), LD, apriori);
                a_row = nan(1,10);
                a_row(1:2:9) = vals;
                a_row(2:2:10) = uncertainty;
                A = cat(1, A, a_row);
            end
            A = A';
            varnames = {'Mn108Fast','MnFast','HyFast','Mn108Slow','MnSlow','HySlow'};
        end
        fns = {'a (mol)','uncert. a','x_0 (km)','uncert. x','mu_x (km)','uncert. mu','sigma_x (km)','uncert. sigma','B (mol)','uncert. B'};
        latex_fns = {'$a$ (mol \chem{NO_2})';'';'$x_0$ (km)';'';'$\mu_x$ (km)';'';'$\sigma_x$ (km)';'';'$B$ (mol \chem{NO_2} km$^{-1}$)';''};
        if isfield(LD,'wind_crit')
            wc = WC.windvel >= LD.wind_crit;
            [emis, uncert_E, tau, uncert_tau] = compute_emg_emis_tau(A(1,:),A(2,:),A(3,:),A(4,:), 'vec', WC.windvel(wc));
            
            A = cat(1,A,emis,uncert_E,tau,uncert_tau);
            fns = cat(2, fns, 'E (Mg NOx/h)','uncert. E', 'tau (h)','uncert. tau');
            latex_fns = cat(1, latex_fns, '$E$ (Mg \chem{NO_x} h$^{-1}$)', {''}, '$\tau_{\mathrm{eff}}$ (h)',{''});
        else
            fprintf('Cannot compute emissions and tau without wind criterion\n');
        end
        
        % print the table
        T = array2table(A,'RowNames',fns,'VariableNames',varnames)
        % Print the latex output
        L = cat(2, latex_fns, num2cell(A));
        fprintf('\nLatex format:\n\n');
        %mat2latex(L,'%#.3g',1)
        mat2latex(L,'u',1)
    end

    function residuals(resid_type) 
        % Let's first just plot the relative sum of squared residuals vs. the fixed
        % value for each parameter or the R-square value. The last option
        % selects which one to use; input 'rel' or nothing to select the
        % ratio of SSresid to optimal and 'r2' to plot R-square values.
        
        resid_type_in = exist('resid_type','var');
        if resid_type_in
            [Ffixed, F, ~, prof_wind_type] = resid_user_input(~resid_type_in);
        else
            [Ffixed, F, resid_type, prof_wind_type] = resid_user_input(~resid_type_in);
        end
        
        %%%%% MAIN FXN %%%%%
        params = {Ffixed(1,:).fixed_par};
        fns = fieldnames(F.ffit);
        for a=1:numel(params)
            x=[Ffixed(:,a).fixed_val];
            switch lower(resid_type)
                case 'resid'
                    y=[Ffixed(:,a).ssresid];
                    yopt = F.ssresid;
                    ytext = 'Sum of squared resid fixed';
                case 'rel'
                    y=[Ffixed(:,a).ssresid] ./ F.ssresid;
                    yopt = 1;
                    ytext = 'Sum of squared resid fixed / sum SR unfixed';
                case 'r2'
                    y = nan(size(Ffixed,1),1);
                    for b=1:numel(y)
                        y(b) = Ffixed(b,a).stats.r2;
                    end
                    yopt = F.stats.r2;
                    ytext = 'R^2';
                case 'r'
                    y = nan(size(Ffixed,1),1);
                    for b=1:numel(y)
                        y(b) = Ffixed(b,a).stats.r;
                    end
                    yopt = F.stats.r;
                    ytext = 'R';
                otherwise
                    fprintf('Residual type not recognized');
                    return
            end
            figure;
            plot(x,y,'ko','markersize',8,'linewidth',2)
            line(F.ffit.(fns{a}),yopt,'linestyle','none','marker','x','markersize',8,'color','k','linewidth',2)
            xlabel(sprintf('Fixed value of %s',params{a}));
            ylabel(ytext)
            set(gca,'fontsize',16)
            title(sprintf('%s residuals - %s fixed, %d points',prof_wind_type,params{a},numel(x)))
        end
        tilefigs;
    end

    function residuals_onepar(par, varargin)
        % A variation on the last one that allows you to plot the residual
        % change with one parameter over many cases. Inputs: par must be a
        % the index of the parameter to plot (a = 1, x_0 = 2, etc). The
        % rest of the inputs must be pairs of Ffixed and F structures.
        if ismember(varargin{end},{'rel','r2','r'})
            resid_type = varargin{end};
            varargin(end) = [];
        else
            resid_type = 'rel';
        end
        params = {varargin{1}(1,:).fixed_par};
        fns = fieldnames(varargin{2}.ffit);
        for a=1:numel(varargin)/2
            Ffixed = varargin{a*2-1};
            F = varargin{a*2};
            x=[Ffixed(:,par).fixed_val];
            switch resid_type
                case 'rel'
                    y=[Ffixed(:,par).ssresid] ./ F.ssresid;
                    yopt = 1;
                    ytext = 'Sum of squared resid fixed / sum SR unfixed';
                case 'r2'
                    y = nan(size(Ffixed,1),1);
                    for b=1:numel(y)
                        y(b) = Ffixed(b,par).stats.r2;
                    end
                    yopt = F.stats.r2;
                    ytext = 'R^2';
                case 'r'
                    y = nan(size(Ffixed,1),1);
                    for b=1:numel(y)
                        y(b) = Ffixed(b,par).stats.r;
                    end
                    yopt = F.stats.r;
                    ytext = 'R';
                otherwise
                    fprintf('Residual type not recognized');
                    return
            end
            figure;
            plot(x,y,'ko','markersize',8,'linewidth',2)
            line(F.ffit.(fns{par}),yopt,'linestyle','none','marker','x','markersize',8,'color','k','linewidth',2)
            xlabel(sprintf('Fixed value of %s',params{par}));
            ylabel(ytext)
            set(gca,'fontsize',16)
            title(sprintf('Residuals - %s fixed, %d points',params{par},numel(x)))
        end
    end

    function cross_effects()
        % Next let's see how the other parameters vary as we fix each
        % parameter in turn
        [Ffixed, F, ~, prof_wind_type] = resid_user_input(false);
        params = {Ffixed(1,:).fixed_par};
        fns = fieldnames(F.ffit);
        for a=1:numel(params)
            figure;
            subplot(2,2,1);
            ind = 1:numel(params);
            ind(ind==a) = [];
            
            x = [Ffixed(:,a).fixed_val];
            for b=1:numel(ind)
                subplot(2,2,b)
                y = nan(1,size(Ffixed,1));
                for c=1:size(Ffixed,1)
                    y(c) = Ffixed(c,a).ffit.(fns{ind(b)});
                end
                line(x,y,'linestyle','none','marker','o','color','k','linewidth',2,'markersize',8);
                line(F.ffit.(fns{a}), F.ffit.(fns{ind(b)}), 'marker','x','color','k','linestyle','none','markersize',8,'linewidth',2)
                xlabel(sprintf('%s (fixed)',Ffixed(c,a).fixed_par))
                ylabel(fns{ind(b)})
                set(gca,'fontsize',16)
            end
            suptitle(sprintf('%s: %s fixed',prof_wind_type,params{a}));
        end
    end

    function show_fits(Ffixed,F,prof_wind_type,fixed_var)
        params = {Ffixed(1,:).fixed_par};
        xx = strcmp(fixed_var,params);
        if sum(xx) == 0
            fprintf('%s not recognized as a parameter in Ffixed\n',fixed_var)
            return
        end
        figure;
        l(1)=line(Ffixed(1,xx).no2x, Ffixed(1,xx).no2ld, 'marker','o','color','k','linestyle','none');
        l(2)=line(F.no2x, F.emgfit, 'color','k','linewidth',2,'linestyle','--');
        lstr{1} = 'Data'; lstr{2} = 'Free fit';
        cols = {'b','r',[0 0.5 0],'m',[1 0.5 0], [0 0.5 0.5]};
        for a=1:size(Ffixed,1)
            cind = mod(a-1,6)+1;
            l(a+2) = line(Ffixed(a,xx).no2x, Ffixed(a,xx).emgfit, 'color', cols{cind}, 'linestyle', '-.', 'linewidth',1);
            lstr{a+2} = sprintf('%s = %.3g', Ffixed(a,xx).fixed_par, Ffixed(a,xx).fixed_val);
        end
        legend(l',lstr)
        title(sprintf('%s - %s fixed',prof_wind_type, fixed_var))
    end

    function mchist(mcpoints, ffit)
        params = fieldnames(ffit);
        for a=1:5
            figure; hist(mcpoints(:,a),50);
            y = get(gca,'ylim');
            line([ffit.(params{a}), ffit.(params{a})], y, 'color', 'r', 'linewidth', 2);
        end
    end

    function stats_table
        % Will load all the simple fit files in the selected directory and
        % print out tables of the requested statistic.
        fit_dir = uigetdir('.','Choose the directory with the simple fit files to query');
        F_ss = dir(fullfile(fit_dir,'*SimpleFits-ssresid*'));
        F_unex = dir(fullfile(fit_dir,'*SimpleFits-unexvar*'));
        
        % Load one file to get some information from
        D = load(fullfile(fit_dir,F_ss(1).name));
        fns = fieldnames(D);
        stats = fieldnames(D.(fns{1}).stats);
        
        % Ask the user which statistic to summarize in the table
        stat = ask_multichoice('Which statistic to summarize?',stats,'list',true);
        
        [tstr_ss, rownames_ss] = make_stats_table(fit_dir, F_ss, stat);
        [tstr_unex, rownames_unex] = make_stats_table(fit_dir, F_unex, stat);
        
        fprintf('SSRESID table:')
        struct2table(tstr_ss, 'RowNames', rownames_ss)
        fprintf('UNEXVAR table:')
        struct2table(tstr_unex, 'RowNames', rownames_unex)
    end

    function uncert = calc_total_fit_uncert(sfit, LD, apriori)
        % Calculate the overall uncertainty of the parameters, following
        % the supplement in Beirle et. al. They assign the following
        % uncertainties:
        %   VCD = 30%, 25% in Lu et al. 2015 NOx:NO2 ratio = 10% Fit
        %   confidence interval = ~10-50% SME (standard mean error?) of fit
        %   results = ~10-40% Choice of b (across wind integration
        %   distance) = 10% Choice of wind fields = 30%
        % I will ignore the SME as Lu et al 2015 does, because that was for
        % Beirle's average over 8 different sectors, rather than the method
        % of aligning wind directions. Lu goes on to specify that NOx:NO2
        % ratio only matters for emissions (E ~ a*w/x0 where w is avg. wind
        % speed) but not burdens (a) or lifetime (tau = x0 / w) and the
        % uncertainty in VCDs does not matter for lifetime. However, I
        % would argue that the uncertainty in VCDs DOES matter for lifetime
        % as a spatial bias in the VCDs could introduce essentially a
        % second lifetime that would convolve with the true lifetime. So
        % this will assume that all the sources of uncertainty count for
        % all the fit parameters, except for the NOx:NO2 ratio.
        num_obs_fn = sprintf('num_obs_%s',apriori);
        if ~isfield(LD,num_obs_fn)
            s_vcd = 0.25;
            if first_warning
                first_warning = false;
                warning('No number of valid observations given in the line density structure, setting VCD uncertainty to 25%')
            end
        else
            % If there is a number of valid observations variable, then we
            % will reduce the uncertainty in the line density by the
            % smallest number present to be conservative.
            s_vcd = 0.25 / sqrt(min(LD.(num_obs_fn)(LD.(num_obs_fn)>0)));
        end
        s_b = 0.1;
        s_wind = 0.3;
        per_uncert = sqrt((sfit.stats.percent_ci95/100).^2 + s_vcd.^2 + s_b.^2 + s_wind.^2);
        uncert = abs(struct2array(sfit.ffit)) .* per_uncert';
    end
end

function [tstr, rownames] = make_stats_table(fit_dir, F, stat)
D = load(fullfile(fit_dir,F(1).name));
fns = fieldnames(D);
apriori = fns;
for a=1:numel(apriori)
    % clean up a priori names for column headers
    apriori{a} = strrep(apriori{a},'f_','');
end
% Prep the structure and cell array to hold the row names
tstr = make_empty_struct_from_cell(apriori,'');
tstr = repmat(tstr,numel(F),1);
rownames = cell(numel(F),1);


for a=1:numel(F)
    D = load(fullfile(fit_dir,F(a).name));
    % Get the city name and wind sep speed
    [s,e] = regexp(F(a).name,'(Atlanta|Birmingham|Montgomery)');
    city_name = F(a).name(s:e);
    [s,e] = regexp(F(a).name,'\dpt\d');
    wind_spd = F(a).name(s:e);
    rownames{a} = sprintf('%s_%s',city_name,wind_spd);
    for f=1:numel(fns)
        tstr(a).(apriori{f}) = nanmean(D.(fns{f}).stats.(stat));
    end
end
end

function [Ffixed, F, resid_type, prof_wind_type] = resid_user_input(ask_resid_type)
%%%%% USER INPUT %%%%%
[filename, filepath] = uigetfile('*VariationalFits*.mat','Choose the file to plot');
if isnumeric(filename) && filename == 0
    return
end

% Get file and apriori to plot
D = load(fullfile(filepath, filename));
fns = fieldnames(D);
apriori = cell(numel(fns)/2,1);
a = 1;
for f=1:numel(fns)
    if ~isempty(strfind(fns{f},'F_'))
        apriori{a} = fns{f}(3:end);
        a=a+1;
    end
end
[sel,ok] = listdlg('ListString',apriori,'SelectionMode','single','PromptString','Choose which to plot');
if ok == 1
    fn1 = strcat('F_',apriori{sel});
    fn2 = strcat('Ffix_',apriori{sel});
    F = D.(fn1);
    Ffixed = D.(fn2);
else
    return
end

% Which residual type to plot
if ask_resid_type
    allowed_resid = {'resid','rel','r2','r'};
    resid_type = ask_multichoice('Which residual type to use?',allowed_resid,'list',true);
else 
    resid_type = '';
end

% Make plot name
[s,e] = regexp(filename,'(Atlanta|Birmingham|Montgomery)');
city_name = filename(s:e);
[s,e] = regexp(filename,'\dpt\d');
wind_spd = regexprep(filename(s:e),'pt','\.');
prof_wind_type = sprintf('%s (%s)',city_name,wind_spd);
end



