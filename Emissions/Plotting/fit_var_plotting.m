function [  ] = fit_var_plotting(plottype, Ffixed, F, varargin )
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

switch lower(plottype)
    case 'residuals'
        residuals(Ffixed,F,varargin{:})
    case 'crosseffects'
        cross_effects(Ffixed,F, varargin{:})
    case 'showfits'
        show_fits(Ffixed, F, varargin{:})
    otherwise
        fprintf('Plot type not recognized\n')
end


    function residuals(Ffixed,F,prof_wind_type)
        % Let's first just plot the relative sum of squared residuals vs. the fixed
        % value for each parameter
        params = {Ffixed(1,:).fixed_par};
        fns = fieldnames(F.ffit);
        for a=1:numel(params)
            x=[Ffixed(:,a).fixed_val];
            y=[Ffixed(:,a).ssresid] ./ F.ssresid;
            figure;
            plot(x,y,'ko','markersize',8,'linewidth',2)
            line(F.ffit.(fns{a}),1,'linestyle','none','marker','x','markersize',8,'color','k','linewidth',2)
            xlabel(sprintf('Fixed value of %s',params{a}));
            ylabel('Sum of squared resid fixed / sum SR unfixed')
            set(gca,'fontsize',16)
            title(sprintf('%s residuals - %s fixed, %d points',prof_wind_type,params{a},numel(x)))
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
end

