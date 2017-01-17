function [ output_args ] = behr_prod_indiv_scatter( indiv_stats, field_to_plot, plot_mode, clim )
%BEHR_PROD_INDIV_SCATTER Plots day by day scatter plots where % diff > 0.5
%   BEHR_PROD_INDIV_SCATTER( indiv_stats, field_to_plot ) makes plots for
%   FIELD_TO_PLOT from INDIV_STATS returned by BEHR_PROD_TEST.
%
%   BEHR_PROD_INDIV_SCATTER( ___, 'perdiff', [cmin, cmax] ) allows you to
%   specify a color axis range.
%
%   BEHR_PROD_INDIV_SCATTER( ___, 'diff' ) plots absolute differences
%   again, still for pixels with a percent difference > 0.5%.
%
%   BEHR_PROD_INDIV_SCATTER( ___, 'diff', [cmin, cmax] )  allows you to
%   specify a color axis range.
%
%   BEHR_PROD_INDIV_SCATTER( ___, 'values' ) will make two plots per day,
%   with new and old values.
%
%   BEHR_PROD_INDIV_SCATTER( ___, 'values', [cmin, cmax] )  allows you to
%   specify a color axis range
%
%   Any of the plotting modes (perdiff, diff, values) can have "all" added
%   as a prefix to show all differences, not just those > 0.5%.

E = JLLErrors;

n = numel(indiv_stats);
if ~exist('plot_mode','var')
    plot_mode = 'perdiff';
end
if ~exist('clim', 'var')
    clim = [];
end

for a=1:n
    if ~regexp(plot_mode,'all')
        xx = abs(indiv_stats(a).(field_to_plot).percent_differences)>0.5;
    else
        xx = true(size(indiv_stats(a).(field_to_plot).percent_differences));
    end
    if ismember(lower(plot_mode), {'perdiff','diff','allperdiff','alldiff'})
        f=figure;
        if ismember(lower(plot_mode), {'perdiff','allperdiff'})
            delta = indiv_stats(a).(field_to_plot).percent_differences(xx);
        elseif ismember(lower(plot_mode), {'diff','alldiff'})
            delta = indiv_stats(a).(field_to_plot).differences(xx);
        end
        scatter(indiv_stats(a).(field_to_plot).Longitude(xx), indiv_stats(a).(field_to_plot).Latitude(xx), 8, delta)
        colorbar;
        if ~isempty(clim)
            caxis(clim)
        end
        states('k','cont');
        title(indiv_stats(a).(field_to_plot).date)
    elseif ismember(lower(plot_mode), {'values','allvalues'})
        f(1) = figure;
        scatter(indiv_stats(a).(field_to_plot).Longitude(xx), indiv_stats(a).(field_to_plot).Latitude(xx), 8, indiv_stats(a).(field_to_plot).value_pairs(xx,1));
        colorbar;
        if ~isempty(clim)
            caxis(clim)
        end
        states('k','cont');
        title(sprintf('New - %s',indiv_stats(a).(field_to_plot).date));
        f(1).Position(1) = f(1).Position(1) - f(1).Position(3)/2; % move to left
        
        f(2) = figure;
        scatter(indiv_stats(a).(field_to_plot).Longitude(xx), indiv_stats(a).(field_to_plot).Latitude(xx), 8, indiv_stats(a).(field_to_plot).value_pairs(xx,2));
        colorbar;
        if ~isempty(clim)
            caxis(clim)
        end
        states('k','cont');
        title(sprintf('Old - %s',indiv_stats(a).(field_to_plot).date));
        f(2).Position(1) = f(2).Position(1) + f(2).Position(3)/2; % move to right
    else
        E.badinput('Plotting mode ''%s'' not recognized', plot_mode);
    end
    drawnow
    fprintf('Press any key for next day\n')
    pause
    for b=1:numel(f)
        close(f(b))
    end
end


end

