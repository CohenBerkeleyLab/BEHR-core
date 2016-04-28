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
    otherwise
        fprintf('Plot type not recognized\n')
end

    function plot_3_fits
        [data_file, data_path] = uigetfile('*LineDensities*.mat','Select the line density file for the conditions to plot');
        LD = load(fullfile(data_path,data_file));
        fitfile = strrep(data_file,'LineDensities','SimpleFits');
        if ~exist(fullfile(data_path,fitfile),'file')
            E.filenotfound('simple fits file');
        else
            SF = load(fullfile(data_path,fitfile));
        end
        
        % Get city name
        [s,e] = regexp(data_file,'(Atlanta|Birmingham|Montgomery)');
        city_name = data_file(s:e);
        
        % Get wind conditions
        WC = load(sprintf('/Users/Josh/Documents/MATLAB/BEHR/Workspaces/Wind speed/%s-Wind-Conditions-1900UTC-5layers.mat',city_name));
        
        figure; hold on
        plot(LD.no2x_mn108slow, LD.no2ld_mn108slow, 'go', 'linewidth', 2);
        plot(LD.no2x_mn108slow, SF.f_mn108slow.emgfit, '--', 'color', [0 0.5 0], 'linewidth', 2);
        plot(LD.no2x_mnslow, LD.no2ld_mnslow, 'o', 'linewidth', 2, 'color', [1 0.5 0]);
        plot(LD.no2x_mnslow, SF.f_mnslow.emgfit, 'r--', 'linewidth', 2);
        plot(LD.no2x_hyslow, LD.no2ld_hyslow, 'co', 'linewidth', 2);
        plot(LD.no2x_hyslow, SF.f_hyslow.emgfit, 'b--', 'linewidth', 2);
        legend('Coarse monthly data','Coarse monthly fit', 'Fine monthly data', 'Fine monthly fit', 'Hybrid data', 'Hybrid fit');
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
        legend('Coarse monthly data','Coarse monthly fit', 'Fine monthly data', 'Fine monthly fit', 'Hybrid data', 'Hybrid fit');
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
        A = cat(1, struct2array(SF.f_mn108fast.ffit), struct2array(SF.f_mnfast.ffit), struct2array(SF.f_hyfast.ffit),...
            struct2array(SF.f_mn108slow.ffit), struct2array(SF.f_mnslow.ffit), struct2array(SF.f_hyslow.ffit))';
        
        fns = {'a (mol)','x_0 (km)','mu_x (km)','sigma_x (km)','B (mol)'};
        latex_fns = {'$a$ (mol \chem{NO_2})';'$x_0$ (km)';'$\mu_x$ (km)';'$\sigma_x$ (km)';'$B$ (mol \chem{NO_2} km$^{-1}$)'};
        if isfield(LD,'wind_crit')
            avg_wind = mean(WC.windvel(WC.windvel >= LD.wind_crit))/1000*3600; % windvel in m/s, convert to km/h
            emis = 1.32 .* A(1,:) .* avg_wind ./ A(2,:);
            
            % Calculate assumed mass of NOx if NOx:NO2 ratio is 1.32
            % MM NO = 30.01 g/mol
            % MM NO2 = 46.01 g/mol
            mol2Mg = (1/1.32 * 46.01 + (1-1/1.32)*30.01)*1e-6;
            emis = emis * mol2Mg;
            
            tau = A(2,:) ./ avg_wind;
            A = cat(1,A,emis,tau);
            fns = cat(2, fns, 'E (Mg NOx/h)', 'tau (h)');
            latex_fns = cat(1, latex_fns, '$E$ (Mg \chem{NO_x} h$^{-1}$)', '$\tau_{\mathrm{eff}}$ (h)');
        else
            fprintf('Cannot compute emissions and tau without wind criterion\n');
        end
        
        % print the table
        T = array2table(A,'RowNames',fns,'VariableNames',{'Mn108Fast','MnFast','HyFast','Mn108Slow','MnSlow','HySlow'})
        % Print the latex output
        L = cat(2, latex_fns, mat2cell_simple(A));
        fprintf('\nLatex format:\n\n');
        mat2latex(L,'%#.3g')
    end

    function residuals(Ffixed,F,prof_wind_type,resid_type)
        % Let's first just plot the relative sum of squared residuals vs. the fixed
        % value for each parameter or the R-square value. THe last option
        % selects which one to use; input 'rel' or nothing to select the
        % ratio of SSresid to optimal and 'r2' to plot R-square values.
        if ~exist('resid_type','var')
            resid_type = 'rel';
        end
        params = {Ffixed(1,:).fixed_par};
        fns = fieldnames(F.ffit);
        for a=1:numel(params)
            x=[Ffixed(:,a).fixed_val];
            switch resid_type
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

    function cross_effects(Ffixed,F,prof_wind_type)
        % Next let's see how the other parameters vary as we fix each
        % parameter in turn
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
end

