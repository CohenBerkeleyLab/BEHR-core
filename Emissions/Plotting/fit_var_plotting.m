function [ varargout ] = fit_var_plotting(plottype, varargin )
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
HOME_DIR = getenv('HOME');

if ~exist('plottype','var')
    plottype = ask_multichoice('Which plot to make?', {'fits'},'list',true);
end

if nargout > 0
    vout = cell(1, nargout);
else
    vout = {};
end

switch lower(plottype)
    case 'fits'
        plot_3_fits(varargin{:})
    case 'varfits'
        plot_fit_envelopes()
    case 'e-tau'
        emis_tau_table();
    case 't-tests'
        [vout{:}] = t_tests(nargout, varargin{:});
    case 'xwind'
        xwind_emis_tau();
    case 'sectors'
        plot_sectors();
    case 'residuals'
        residuals(varargin{:})
    case 'residuals-onepar'
        residuals_onepar_wrapper();
    case 'residuals-onepar-manual'
        residuals_onepar(varargin{:})
    case 'crosseffects'
        cross_effects(varargin{:})
    case 'showfits'
        show_fits(varargin{:})
    case 'mchist'
        mchist(varargin{:});
    case 'stats'
        stats_table(varargin{:});
    case 'autocorr'
        autocorrelation_test(varargin{:});
    otherwise
        fprintf('Plot type not recognized\n')
end

varargout = vout;

    function plot_3_fits(filename)
        if ~exist('filename','var')
            [data_file, data_path] = uigetfile('*SimpleFits*.mat','Select the fits file for the conditions to plot');
            if data_file == 0
                fprintf('Cancelled\n')
                return
            end
        else
            [data_path, data_file, ext] = fileparts(filename);
            data_file = [data_file, ext];
        end
        SF = load(fullfile(data_path,data_file));
        ldfile = regexprep(data_file,'SimpleFits(-ssresid|-unexvar)*','LineDensities');
        if ~exist(fullfile(data_path,ldfile),'file')
            E.filenotfound('line density file');
        else
            LD = load(fullfile(data_path,ldfile));
        end
        
        % Get city name
        [s,e] = regexp(data_file,'(Atlanta|Birmingham|Montgomery)');
        city_name = data_file(s:e);
        
        % Get wind conditions
        WC = load(fullfile(HOME_DIR,'Documents','MATLAB','BEHR','Workspaces','Wind speed',sprintf('%s-Wind-Conditions-1900UTC-5layers-earthrel.mat',city_name)));
        
        % Determine whether to plot the WRF ones too
        plot_wrf = false;
        if all(isfield(LD,{'no2x_wrffast','no2x_wrfslow'})) && isDisplay
            plot_wrf_ans = questdlg('Include WRF line densities?');
            if strcmpi(plot_wrf_ans,'Yes')
                plot_wrf = true;
            elseif strcmpi(plot_wrf_ans,'Cancel')
                return
            end
        end
        
        figure; hold on
        lstr = {'Line densities using coarse monthly a priori','Fit - coarse monthly', 'Line densities using fine monthly a priori', 'Fit - fine monthly', 'Line densities using fine daily a priori', 'Fit - fine daily'};
        plot(LD.no2x_mn108slow, LD.no2ld_mn108slow, 'ko', 'linewidth', 2);
        plot(LD.no2x_mn108slow, SF.f_mn108slow.emgfit, '--', 'color', 'k', 'linewidth', 2);
        plot(LD.no2x_mnslow, LD.no2ld_mnslow, 'o', 'linewidth', 2, 'color', [1 0.5 0]);
        plot(LD.no2x_mnslow, SF.f_mnslow.emgfit, 'r--', 'linewidth', 2);
        plot(LD.no2x_hyslow, LD.no2ld_hyslow, 'co', 'linewidth', 2);
        plot(LD.no2x_hyslow, SF.f_hyslow.emgfit, 'b--', 'linewidth', 2);
        if plot_wrf
            plot(LD.no2x_wrfslow, LD.no2ld_wrfslow, 'go', 'linewidth',2);
            plot(LD.no2x_wrfslow, SF.f_wrfslow.emgfit, '--','linewidth',2,'color',[0 0.5 0])
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
        plot(LD.no2x_mn108fast, LD.no2ld_mn108fast, 'ko', 'linewidth', 2);
        plot(LD.no2x_mn108fast, SF.f_mn108fast.emgfit, '--', 'color', 'k', 'linewidth', 2);
        plot(LD.no2x_mnfast, LD.no2ld_mnfast, 'o', 'linewidth', 2, 'color', [1 0.5 0]);
        plot(LD.no2x_mnfast, SF.f_mnfast.emgfit, 'r--', 'linewidth', 2);
        plot(LD.no2x_hyfast, LD.no2ld_hyfast, 'co', 'linewidth', 2);
        plot(LD.no2x_hyfast, SF.f_hyfast.emgfit, 'b--', 'linewidth', 2);
        if plot_wrf
            plot(LD.no2x_wrffast, LD.no2ld_wrffast, 'go', 'linewidth',2);
            plot(LD.no2x_wrffast, SF.f_wrffast.emgfit, '--','linewidth',2,'color',[0 0.5 0])
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
        %mat2latex(L,'%#.2g',1)
        mat2latex(L,'u10',1)
    end

    function plot_fit_envelopes
        [data_files, data_path, varfns_to_keep, simpfns_to_keep, param_inds] = input_varfit_params;
        
        for f=1:numel(data_files)
            % Load the selected file and separate the variables into the varied
            % fits and the simple fits
            V = load(fullfile(data_path, data_files{f}));
            varfns = glob(fieldnames(V),'Ffix.*');
            simpfns = glob(fieldnames(V),'F_.*');
            
            % Extract city name and wind speed from file name
            city_name = regexp(data_files{f}, '(Atlanta|Birmingham|Montgomery)','match','once');
            wind_spd = regexp(data_files{f}, '\dpt\d','match','once');
            
            % Just in case the order differs in other files, we'll actually
            % compare the variables to keep
            varfns(~ismember(varfns, varfns_to_keep)) = [];
            simpfns(~ismember(simpfns, simpfns_to_keep)) = [];
            
            % For each apriori and each fitting parameter, plot the line
            % densities and the best fit, then plot the fit as each parameter
            % is tweaked by +/- 0.5 and +/- 1 sigma, represented by which
            % indices to look at
            
            for a=1:numel(varfns)
                Ffix = V.(varfns{a});
                n_sd = size(Ffix,1);
                ideal_frac_sigma = [-1 -0.26 0.26 1];
                inds = round((ideal_frac_sigma*0.5 + 0.5)*n_sd);
                inds = min(max(inds,1),n_sd);
                true_frac_sigma = (2*inds - n_sd)./n_sd;
                
                cols = {'r','r','b','b'};
                lstyles = {'-','--','--','-'};
                n_inds = numel(inds);
                
                for b=param_inds %size(Ffix,2) % the fixed parameter changes across the second dimension
                    figure;
                    plot(Ffix(1,b).no2x, Ffix(1,b).no2ld, 'ko'); % line densities same for all
                    hold on
                    plot(Ffix(1,b).no2x, V.(simpfns{a}).emgfit, 'k--','linewidth',2);
                    
                    % Get the varied fits
                    l = gobjects(n_inds,1);
                    vals = nan(n_inds,1);
                    for c=1:n_inds
                        i = inds(c);
                        l(c) = line(Ffix(i,b).no2x, Ffix(i,b).emgfit, 'color', cols{c}, 'linewidth', 2, 'linestyle', lstyles{c});
                        vals(c) = Ffix(i,b).fixed_val;
                    end
                    legend(l, sprintfmulti('Fit %1$.2f\\sigma_{%2$s} (%2$s = %3$.3g)', true_frac_sigma, Ffix(1,b).fixed_par, vals));
                    set(gca,'fontsize',14)
                    xlabel('Dist (km)')
                    ylabel('Line dens. (mol km^{-1})');
                    apri_name = strrep(simpfns{a},'F_','');
                    title(sprintf('%s: %s, %s: %s', Ffix(1,b).fixed_par, city_name, wind_spd, apri_name))
                end
            end
        end
    end

    function emis_tau_table
        % This function allows you to create a table of just emissions and
        % lifetimes (with uncertainty) for fast winds. Will divide all
        % selected files by city, a priori, and wind speed bin. Note that
        % it will not include the WRF sum of NEI emissions (computed by
        % misc_wrf_wind_plots('sum-emis')) so you'll need to add that
        % yourself.
        [fnames, pname] = uigetfile('*SimpleFits*.mat','Choose all files to use','multiselect','on');
        table_uncert = questdlg('Use confidence intervals or standard deviations?','CI or SD','CI','SD','CI');
        table_uncert = lower(table_uncert);
        
        % First find all cities represented, then all the wind speeds.
        [cities, wind_speeds, wind_labels] = get_cities_and_winds_from_fnames(fnames);
        
        % the apriori cases to look at
        apri = {'f_mn108fast','f_mnfast','f_hyfast'};
        
        [emis, tau] = tabulate_emis_tau(fnames, pname, cities, wind_speeds, apri, table_uncert);
        
        latex_header1 = '& ';
        table_header = cell(1,numel(cities)*3);
        for a=1:numel(cities)
            latex_header1 = [latex_header1, sprintf('& \\multicolumn{3}{c}{%s} ',cities{a})];
            for c=1:numel(apri)
                i = (a-1)*3+c;
                table_header{i} = sprintf('%s_%s',cities{a},apri{c}(3:end));
            end
        end
        latex_header1 = [latex_header1, '\\'];
        latex_header2 = '& Wind speed bin & \parbox[t]{1.5cm}{Monthly\\108 km} & \parbox[t]{1.5cm}{Monthly\\12 km} & \parbox[t]{1.5cm}{Daily\\12 km} & \parbox[t]{1.5cm}{Monthly\\108 km} & \parbox[t]{1.5cm}{Monthly\\12 km} & \parbox[t]{1.5cm}{Daily\\12 km} \\ \middlehline';
        
        
        emis_rownames = cell(size(wind_labels));
        emis_rownames{1} = '\multirow{3}{*}{E (Mg \chem{NO_x} h$^{-1}$)}';
        emis_rownames(2:end) = {''};
        emis_rownames = cat(2, emis_rownames, wind_labels);
        
        tau_rownames = cell(size(wind_labels));
        tau_rownames{1} = '\multirow{3}{*}{$\tau$ (h)}';
        tau_rownames(2:end) = {''};
        tau_rownames = cat(2, tau_rownames, wind_labels);
        
        if strcmp(table_uncert, 'ci')
            table_rownames = wind_labels;
            for a=2:2:numel(table_rownames)
                table_rownames(a) = {sprintf('%s uncert.', wind_labels{a-1})};
            end
        else
            table_rownames = cell(size(emis,1),1);
            for a=1:3:numel(table_rownames)
                wl = wind_labels{(a-1)*2/3+1};
                table_rownames{a} = wl;
                table_rownames{a+1} = ['sigma ',wl];
                table_rownames{a+2} = ['N(DOF) ',wl];
            end
        end
        
        emis_table = array2table(emis, 'VariableNames', table_header, 'RowNames', table_rownames)
        tau_table = array2table(tau, 'VariableNames', table_header, 'RowNames', table_rownames)
        
        if strcmp(table_uncert,'ci')
            emis_cell = cat(2, emis_rownames, num2cell(emis));
            tau_cell = cat(2, tau_rownames, num2cell(tau));
            
            fprintf('\n\n Latex format \n\n')
            
            fprintf('%s\n',latex_header1);
            fprintf('%s\n',latex_header2);
            mat2latex(emis_cell,'u10',1);
            fprintf('\b\b \\middlehline\n');
            mat2latex(tau_cell,'u10',1);
        end
        
    end

    function varargout = t_tests(nout, fnames)
        
        t_type = ask_multichoice('Which t-test to use?', {'Two sample','Paired'}, 'list', true);
        
        if ~exist('fnames','var')
            [fnames, pname] = uigetfile('*SimpleFits*.mat','Choose all files to use','multiselect','on');
        else
            if ~iscell(fnames)
                fnames = {fnames};
            end
            pname = '.';
        end
        [cities, wind_speeds] = get_cities_and_winds_from_fnames(fnames);
        
        if ~iscell(fnames)
            fnames = {fnames};
        end
        
        % the apriori cases to look at
        apri = {'f_mn108fast','f_mnfast','f_hyfast'};
        
        [emis, tau] = tabulate_emis_tau(fnames, pname, cities, wind_speeds, apri, 'sd');
        
        % Prep the array that the t-values will go into
        nsets = numel(cities)*numel(wind_speeds)*numel(apri);
        S = repmat(struct('city','','wind_speed','','apri',''), numel(cities)*numel(wind_speeds)*numel(apri), 1);
        set_names = cell(1,nsets);
        emis_tarray = -inf(nsets, nsets); % using inf rather than nan so that I can tell the difference between unfilled value and t = nan
        tau_tarray = -inf(nsets, nsets);
        
        emis_vals = nan(1,nsets);
        emis_uncert = nan(1,nsets);
        emis_dofs = nan(1,nsets);
        
        tau_vals = nan(1,nsets);
        tau_uncert = nan(1,nsets);
        tau_dofs = nan(1,nsets);
        
        i=1;
        for a=1:numel(cities)
            for b=1:numel(wind_speeds)
                for c=1:numel(apri)
                    S(i).city = cities{a};
                    S(i).wind_speed = wind_speeds{b};
                    S(i).apri = apri{c};
                    set_names{i} = sprintf('%s_%s_%s', S(i).city, S(i).wind_speed, S(i).apri(3:end));
                    
                    % Need to translate indicies to get the values from the
                    % emis and tau arrays which alternate value, SD, and
                    % nDOF down the first dimension for each wind speed and
                    % have a priori in order specified grouped by city
                    % across the second dimension.
                    x = (b-1)*3 + 1; % points to the row for the wind speed
                    y = (a-1)*numel(apri) + c; % points to the column for the city and apriori
                    
                    emis_vals(i) = emis(x,y);
                    emis_uncert(i) = emis(x+1,y);
                    emis_dofs(i) = emis(x+2,y);
                    tau_vals(i) = tau(x,y);
                    tau_uncert(i) = tau(x+1,y);
                    tau_dofs(i) = tau(x+2,y);
                    
                    i=i+1;
                end
            end
        end
        
        if strcmpi(t_type, 'Two sample')
            % Now the painful bit, loop over every set and compare with every
            % other set. Not duplicating calculations. Yay loops. I will use
            % the method from section 4-3 case 2 of Harris "Quantitative
            % Chemical Analysis" 8th edition.
            
            emis_dof_tablecell = cell(size(emis_tarray));
            tau_dof_tablecell = cell(size(tau_tarray));
            
            for a=1:nsets
                for b=a:nsets
                    % Emissions first
                    s_pooled = sqrt( (emis_uncert(a)^2 * emis_dofs(a) + emis_uncert(b)^2 * emis_dofs(b)) / (emis_dofs(a) + emis_dofs(b)) );
                    % Since this is the number of measurements, not degrees of
                    % freedom, put back in the 5 that were "taken up" by
                    % fitting 5 parameters.
                    n_pooled = sqrt( (emis_dofs(a)+5) * (emis_dofs(b)+5) / (emis_dofs(a) + 5 + emis_dofs(b) + 5) );
                    emis_tarray(a,b) = abs(emis_vals(a) - emis_vals(b))/s_pooled * n_pooled;
                    emis_dof_tablecell{a,b} = sprintf('%d (t_crit95 = %g)', emis_dofs(a) + emis_dofs(b), tinv(0.975, emis_dofs(a) + emis_dofs(b)));
                    
                    % Lifetimes second
                    s_pooled = sqrt( (tau_uncert(a)^2 * tau_dofs(a) + tau_uncert(b)^2 * tau_dofs(b)) / (tau_dofs(a) + tau_dofs(b)) );
                    n_pooled = sqrt( (tau_dofs(a)+5) * (tau_dofs(b)+5) / (tau_dofs(a) + 5 + tau_dofs(b) + 5) );
                    tau_tarray(a,b) = abs(tau_vals(a) - tau_vals(b))/s_pooled * n_pooled;
                    tau_dof_tablecell{a,b} = sprintf('%d (t_crit95 = %g)', tau_dofs(a) + tau_dofs(b), tinv(0.975, tau_dofs(a) + tau_dofs(b)));
                end
            end
            
            % Make up the tables, print them to screen if no output requested
            emis_table = array2table(emis_tarray, 'VariableNames', set_names, 'RowNames', set_names);
            emis_dof_table = cell2table(emis_dof_tablecell, 'VariableNames', set_names, 'RowNames', set_names);
            tau_table = array2table(tau_tarray, 'VariableNames', set_names, 'RowNames', set_names);
            tau_dof_table = cell2table(tau_dof_tablecell, 'VariableNames', set_names, 'RowNames', set_names);
            
            if nout < 1
                emis_table %#ok<*NOPRT>
                emis_dof_table
                tau_table
                tau_dof_table
            else
                varargout{1} = emis_table;
                varargout{2} = emis_dof_table;
                varargout{3} = tau_table;
                varargout{4} = tau_dof_table;
            end
        elseif strcmpi(t_type, 'Paired')
            % Since the choice of a priori is the testing parameter, we
            % will look at all possible combinations of two a priori,
            % taking each city/wind division pair as a separate test. While
            % not truly independent, it is a sort of bootstrapping in a
            % sense, where we choose a subset of the data and see if we get
            % the same answer.
            %
            % The tabulation returns 0s if the particular city/wind
            % division wasn't included. This happens if you pick a
            % different number of wind divisions for one of the cities, or
            % something similar. Do not include these results.
            
            zz = emis_dofs == 0 & tau_dofs == 0;
            S(zz) = [];
            emis_vals(zz) = [];
            tau_vals(zz) = [];
            
            apri_combos = nchoosek(apri, 2);
            set_cities = {S.city};
            set_winds = {S.wind_speed};
            set_apri = {S.apri};
            
            t_calc_emis = nan(size(apri_combos,1), 1);
            t_crit_emis = nan(size(apri_combos,1), 1);
            t_calc_tau = nan(size(apri_combos,1), 1);
            t_crit_tau = nan(size(apri_combos,1), 1);
            
            varnames = {'E t_calc', 'E t_crit', 'tau t_calc', 'tau t_crit'};
            rownames = cell(1, size(apri_combos,1));
            
            for a=1:size(apri_combos,1)
                rownames{a} = sprintf('%s v. %s', apri_combos{a,1}, apri_combos{a,2});
                
                emis_dels = nan(1,numel(S) / numel(apri));
                tau_dels = nan(1,numel(S) / numel(apri));
                i = 1;
                % Find each city/wind combo that for the first a priori
                for b=1:numel(S)
                    if strcmp(S(b).apri, apri_combos{a,1})
                        % Find the result that has the same city and wind
                        % division but the second a priori
                        xx = strcmpi(set_cities, S(b).city) & strcmpi(set_winds, S(b).wind_speed) & strcmpi(set_apri, apri_combos{a,2});
                        emis_dels(i) = abs(emis_vals(b) - emis_vals(xx));
                        tau_dels(i) = abs(tau_vals(b) - tau_vals(xx));
                        i = i+1;
                    end
                end
                t_calc_emis(a) = nanmean(emis_dels) / nanstd(emis_dels) * sqrt(numel(emis_dels));
                t_crit_emis(a) = tinv(0.975, numel(emis_dels) - 1);
                t_calc_tau(a) = nanmean(tau_dels) / nanstd(tau_dels) * sqrt(numel(tau_dels));
                t_crit_tau(a) = tinv(0.975, numel(tau_dels) - 1);
            end
            
            table_out = array2table([t_calc_emis, t_crit_emis, t_calc_tau, t_crit_tau], 'VariableNames', varnames, 'RowNames', rownames);
            if nout < 1
                table_out
            else
                varargout{1} = table_out;
            end
        else
            E.notimplemented(t_type)
        end
        
    end

    function autocorrelation_test(filename)
        if ~exist('filename','var')
            [data_file, data_path] = uigetfile('*SimpleFits*.mat','Select the fits file for the conditions to plot');
            if data_file == 0
                fprintf('Cancelled\n')
                return
            end
        else
            [data_path, data_file, ext] = fileparts(filename);
            data_file = [data_file, ext];
        end
        SF = load(fullfile(data_path,data_file));
        ldfile = regexprep(data_file,'SimpleFits(-ssresid|-unexvar)*','LineDensities');
        if ~exist(fullfile(data_path,ldfile),'file')
            E.filenotfound('line density file');
        else
            LD = load(fullfile(data_path,ldfile));
        end
        
        % For each fast wind series, calculate the residuals.
        suffixes = {'hyfast','mnfast','mn108fast'};
        resids = cell(size(suffixes));
        for a = 1:numel(suffixes)
            line_dens = LD.(sprintf('no2ld_%s',suffixes{a}));
            fit_line = SF.(sprintf('f_%s',suffixes{a})).emgfit;
            resids{a} = fit_line - line_dens;
        end
        
        % Durbin-Watson statistic (d). Critical values are given in tables
        % A6 and A7 of Chatterjee and Hadi (Regression Analysis by Example,
        % 5th ed.) or A1 and A2 of
        % http://www.dm.unibo.it/~simoncin/Durbin_Watson_tables.pdf 
        %
        % If d < d_L, there is autocorrelation
        % If d_L < d < d_U, the test is inconclusive
        % If d > d_U, there is no autocorrelation
        % 
        % Also make residual plots, should see no long runs of positive and
        % negative residuals.
        
        for a=1:numel(resids)
            d = nansum2((resids{a}(2:end) - resids{a}(1:end-1)).^2)./nansum2(resids{a}.^2);
            rho = nansum2(resids{a}(2:end) .* resids{a}(1:end-1)) ./ nansum2(resids{a}.^2);
            pm = repmat('+',1,numel(resids{a}));
            pm(resids{a}<0) = '-';
            nruns = sum(resids{a}(2:end) .* resids{a}(1:end-1) < 0)+1;
            npos = sum(resids{a}>0);
            nneg = sum(resids{a}<0);
            mu = 2*npos*nneg ./ (npos+nneg) + 1;
            sigma2 = 2 * npos * nneg * (2*npos*nneg - npos - nneg) ./ ((npos + nneg).^2 * (npos + nneg - 1));
            fprintf('%s d = %.4f; rho = %.4f\n', suffixes{a}, d, rho);
            fprintf('Runs: %s\n', pm);
            fprintf('\t %d runs; expected %.2f +/- %.2f\n', nruns, mu, sigma2);
            figure;
            plot(resids{a}, 'ko-');
            xlabel('Index');
            ylabel('Residual');
            title(suffixes{a});
            line([0 numel(resids{a})], [0 0], 'color', 'r', 'linestyle', '--');
        end
        
    end

    function xwind_emis_tau()
        % This function will compare the lifetimes and emissions derived
        % from line densities using 4 different across-wind integration
        % distances for the three a priori at fast winds. Different figures
        % will be made for each wind bin defined. This assumes that each
        % across wind distance is in it's own subfolder.
        
        cities = {'Atlanta','Birmingham'};
        wind_speeds = {'3pt0','4pt0','5pt0'}; % must be as written in the file names
        apriori = {'f_mn108fast','f_mnfast','f_hyfast'};
        apri_names = {'Coarse monthly','Fine monthly','Fine daily'}; % for use in legend. should be same size as apriori
        
        xwind_subfolders = {'25km-side-earthrel','50km-side-earthrel','75km-side-earthrel','100km-side-earthrel','120km-side-earthrel'}; % should have the number of km/degrees in there somewhere
        xwind_root = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/EMG fits/Autorun/FullDaily-NumObs';
        
        emis_tables = cell(size(xwind_subfolders));
        tau_tables = cell(size(xwind_subfolders));
        
        for a=1:numel(xwind_subfolders)
            F = dir(fullfile(xwind_root, xwind_subfolders{a}, '*SimpleFits*.mat'));
            fnames = {F.name};
            pname = fullfile(xwind_root, xwind_subfolders{a});
            [emis_tables{a}, tau_tables{a}] = tabulate_emis_tau(fnames, pname, cities, wind_speeds, apriori);
        end
        
        % Now we need to reorganize these so that we have vectors for each
        % city, each a priori, each wind speed with an emission and tau for
        % each across wind distance.
        %
        % Each table is organized by wind bin down the first dimension,
        % alternating value and uncertainty. Along the second dimensions,
        % the a priori are grouped by city, in the order specified by the
        % apriori cell array.
        %
        % Basically we'll concatenate the tables together to form the
        % vectors so that they'll be organized in the same way, just with
        % vectors in each place instead of individual
        
        emis = cat(3, emis_tables{:});
        tau = cat(3, tau_tables{:});
        
        % Convert the subfolder names into numbers. Remember though, when
        % we define the box width, we actually do so in degrees and just
        % name the files by km by assuming 100 km to a degree.
        xwind_dist = regexp(xwind_subfolders, '\d*', 'match', 'once');
        for a=1:numel(xwind_dist)
            xwind_dist{a} = str2double(xwind_dist{a});
        end
        xwind_dist = cell2mat(xwind_dist)/100*2; % back in degrees now, also the whole distance (kind of diameter vs radius)
        
        % Marker and color will differentiate a priori.  Match a priori
        % color to the EMG fit plots.
        markers = {'o','^','x'};
        colors = {'k','r','b'};
        
        % Now plot, one per wind bin. Stagger in the x-direction slightly
        % to make it easier to tell the error bars apart.
        for a=1:2:size(emis,1)
            f_emis = figure;
            l_emis = gobjects(numel(apriori),1);
            
            f_tau = figure;
            l_tau = gobjects(numel(apriori),1);
            
            stagger_range = 0.05;
            stagger = linspace(-1,1,numel(apriori))*stagger_range;
            for c=1:numel(apriori)
                for b=1:numel(cities)
                    ax_emis = subplot(numel(cities),1,b,'parent',f_emis);
                    i = (b-1)*numel(apriori) + c;
                    if b==1
                        l_emis(c) = line(xwind_dist+stagger(c), squeeze(emis(a,i,:)), 'color', colors{c}, 'linewidth', 2, 'marker', markers{c}, 'markersize', 8, 'parent', ax_emis);
                    else
                        line(xwind_dist+stagger(c), squeeze(emis(a,i,:)), 'color', colors{c}, 'linewidth', 2, 'marker', markers{c}, 'markersize', 8, 'parent', ax_emis);
                    end
                    scatter_errorbars(xwind_dist+stagger(c), squeeze(emis(a,i,:)), squeeze(emis(a+1,i,:)), 'color', colors{c}, 'linewidth', 1, 'parent', ax_emis);
                    
                    if c == numel(apriori)
                        if b == 1
                            legend(ax_emis, l_emis, apri_names, 'Location', 'northwest');
                        end
                        title(ax_emis, sprintf('%s: %s',cities{b},wind_speeds{(a+1)/2}));
                        if b == numel(cities)
                            xlabel(ax_emis, 'Across wind integration distance (degrees)');
                        end
                        ylabel(ax_emis, sprintf('Emissions\n(Mg NO_x hr^{-1})'))
                        xlim(ax_emis, [min(xwind_dist)-2*stagger_range, max(xwind_dist)+2*stagger_range]);
                        set(ax_emis, 'xtick', xwind_dist);
                        set(ax_emis, 'ygrid', 'on');
                        set(ax_emis,'fontsize',14);
                    end
                    
                    ax_tau = subplot(numel(cities),1,b,'parent',f_tau);
                    if b==1
                        l_tau(c) = line(xwind_dist+stagger(c), squeeze(tau(a,i,:)), 'color', colors{c}, 'linewidth', 2, 'marker', markers{c}, 'markersize', 8, 'parent', ax_tau);
                    else
                        line(xwind_dist+stagger(c), squeeze(tau(a,i,:)), 'color', colors{c}, 'linewidth', 2, 'marker', markers{c}, 'markersize', 8, 'parent', ax_tau);
                    end
                    
                    scatter_errorbars(xwind_dist+stagger(c), squeeze(tau(a,i,:)), squeeze(tau(a+1,i,:)), 'color', colors{c}, 'linewidth', 1, 'parent', ax_tau);
                    if c == numel(apriori)
                        if b == 1
                            legend(l_tau, apri_names, 'Location', 'northwest');
                        end
                        title(ax_tau, sprintf('%s: %s',cities{b},wind_speeds{(a+1)/2}))
                        if b == numel(cities)
                            xlabel(ax_tau, 'Across wind integration distance (degrees)');
                        end
                        ylabel(ax_tau, 'Lifetime (h)');
                        xlim(ax_tau, [min(xwind_dist)-2*stagger_range, max(xwind_dist)+2*stagger_range]);
                        set(ax_tau, 'xtick', xwind_dist);
                        set(ax_tau, 'ygrid', 'on');
                        set(ax_tau,'fontsize',14);
                    end
                end
            end
            
            
            
            
        end
        
    end

    function plot_sectors
        % Very kludgy function to plot the 8 wind sectors for one city,
        % wind speed, a priori at once.
        
        [fname, pname] = uigetfile('*LineDensities*.mat','Choose a line density file', 'MultiSelect', 'off');
        if isnumeric(fname) && fname == 0
            fprintf('Cancelled\n')
            return
        end
        
        F = load(fullfile(pname, fname));
        apri_opts = {'hyfast','hyslow','mnfast','mnslow','mn108fast','mn108slow'};
        apri_ind = listdlg('liststring',apri_opts,'SelectionMode','single');
        apri = apri_opts{apri_ind};
        
        no2_x = F.(sprintf('no2x_%s',apri));
        no2_ld = F.(sprintf('no2ld_%s',apri));
        
        fns = fieldnames(no2_x);
        for f=1:numel(fns)
            figure;
            plot(no2_x.(fns{f}), no2_ld.(fns{f}), 'k');
            title(fns{f});
        end
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

    function residuals_onepar_wrapper
        % Gets the required input for residuals_onepar from user responses
        % rather than having the user manually input all the structures.
        [data_files, data_path, varfns, simpfns, par_inds] = input_varfit_params;
        
        % We can only have one parameter, so make sure that's the case.
        if numel(par_inds) ~= 1
            E.badinput('Can only select one parameter to plot')
        end
        
        % Also get which residual type we're going to plot
        allowed_resid_types = {'rel','r2','r'};
        resid_ind = listdlg('ListString', allowed_resid_types, 'PromptString', 'Choose residual type to plot', 'selectionmode','single');
        resid_type = allowed_resid_types{resid_ind};
        
        
        % Now get each variable pair from all the files requested
        vars = cell(1, numel(data_files)*numel(varfns)*2);
        i=1;
        for a=1:numel(data_files)
            V = load(fullfile(data_path, data_files{a}));
            for b=1:numel(varfns)
                vars{i} = V.(varfns{b});
                vars{i+1} = V.(simpfns{b});
                i=i+2;
            end
        end
        
        % Now call the function that actually does the plotting
        residuals_onepar(par_inds, vars{:}, resid_type);
    end

    function residuals_onepar(par, varargin)
        % A variation on the last one that allows you to plot the residual
        % change with one parameter over many cases. Inputs: par must be a
        % the index of the parameter to plot (a = 1, x_0 = 2, etc).
        % Including one of the strings 'rel', 'r2', or 'r' will change what
        % measure of residuals is plotted. The rest of the inputs must be
        % pairs of Ffixed and F structures.

        if ischar(varargin{end}) && ismember(varargin{end},{'rel','r2','r'})
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
            xu = (x - F.ffit.(fns{par})) ./ F.stats.sd(par);
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
            line(F.ffit.(fns{par}),yopt,'linestyle','none','marker','x','markersize',12,'color','r','linewidth',2)
            
            
            % Compute each x tick as a multiple of the uncertainty away
            % from the optimum value and add this information to the labels
            set(gca,'fontsize',16)
            xt = get(gca,'xtick');
            xl = get(gca,'xlim');
            set(gca,'xtickmode','manual');
            set(gca,'xlim',xl);
            set(gca,'xtick',xt); % for some reason, sometime setting this to manual will remove all ticks
            frac_u = (xt - F.ffit.(fns{par})) ./ F.stats.sd(par);
            xt_labels = get(gca,'xticklabel');
            % Check if the labels are a different power than the tick
            test = str2double(xt_labels) ./ xt';
            test = unique(test(~isnan(test)));
            if numel(test) ~= 1
                E.notimplemented('different label exponents')
            elseif test ~= 1
                p = -log10(test);
                xt_labels = sprintfmulti('%s \\times 10^{%d}\\newline (%.2f{\\itu})', xt_labels, p, frac_u);
            else
                xt_labels = sprintfmulti('%s (%.2f{\\itu})', xt_labels, frac_u);
            end
            set(gca,'xticklabel',xt_labels);
            
            xlabel(sprintf('Fixed value of %s\n(dist. from optimum)',params{par}));
            ylabel(ytext)
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

    function stats_table(dir_in)
        % Will load all the simple fit files in the selected directory and
        % print out tables of the requested statistic.
        if ~exist('dir_in','var')
            fit_dir = uigetdir('.','Choose the directory with the simple fit files to query');
        else
            fit_dir = dir_in;
        end
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
        %   VCD = 30%, 25% in Lu et al. 2015 
        %   NOx:NO2 ratio = 10% Fit
        %   confidence interval = ~10-50% 
        %   SME (standard mean error?) of fit results = ~10-40% 
        %   Choice of b (across wind integration distance) = 10% 
        %   Choice of wind fields = 30%
        %
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
        fit_params = struct2array(sfit.ffit);
        frac_uncerts = sfit.stats.percent_ci95'/100;
        if ~isfield(LD,num_obs_fn)
            uncert = calc_fit_param_uncert(fit_params, frac_uncerts, 'warn', first_warning);
            if first_warning
                first_warning = false;
            end
        else
            % If there is a number of valid observations variable, then we
            % will reduce the uncertainty in the line density by the
            % smallest number present to be conservative.
            num_obs = LD.(num_obs_fn);
            uncert = calc_fit_param_uncert(fit_params, frac_uncerts, num_obs, 'warn', first_warning);
        end
    end

    function [data_files, data_path, varapri_to_keep, simpapri_to_keep, param_inds] = input_varfit_params()
        % This subfunction will return variational fit fits files and their
        % directory, the a priori variables in those files to use, and the
        % indices of which parameters to plot (i.e. a = 1, x0 = 2...)
        [data_files, data_path] = uigetfile('*VariationalFits*.mat','Select the fits file for the conditions to plot','multiselect','on');
        if isnumeric(data_files) && data_files == 0
            E.userCancel;
        elseif ischar(data_files)
            data_files = {data_files};
        end
        
        V = load(fullfile(data_path, data_files{1}));
        varfns = glob(fieldnames(V),'Ffix.*');
        simpfns = glob(fieldnames(V),'F_.*');
        
        apris = listdlg('ListString',simpfns,'PromptString','Choose the apriori/wind bins to plot');
        if isempty(apris)
            E.userCancel;
        end
        % Assume that the apriori are in the same order in both lists and
        % cut down to just the ones we selected.
        varapri_to_keep = varfns(apris);
        simpapri_to_keep = simpfns(apris);
        
        params = {V.(varfns{1})(1,:).fixed_par};
        param_inds = listdlg('ListString',params,'PromptString','Choose which parameter(s) to vary');
        if isempty(param_inds)
            E.userCancel;
        end
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

function [emis, tau] = tabulate_emis_tau(fnames, pname, cities, wind_speeds, apri, uncert)
% -- fnames must be the file names of the simple fits files.
% -- pnames is the path name for where those files are located, it must also
% contain the line density files.
% -- cities is a cell array of city names
% -- wind_speeds is a cell array of wind speeds as strings as written in
% the file names.
% -- apri is a cell array of the variables within the simple fits files
% that contain the fits for a particular a priori - they'll be e.g.
% "f_hyfast"
% -- uncert should be 'ci' (confidence interval) or 'sd' (standard
% deviation). Defaults to 'ci' if not given.

E = JLLErrors;
HOMEDIR = getenv('HOME');

if ~exist('uncert','var')
    uncert = 'ci';
else
    uncert = lower(uncert);
    allowed_uncerts = {'ci','sd'};
    if ~ismember(uncert, allowed_uncerts)
        E.badinput('UNCERT must be one of %s', strjoin(allowed_uncerts,', '));
    end
end

if ~iscell(fnames)
    fnames = {fnames};
end

% the final table will have 3 a priori per city and will have
% uncertainty alternate with value along the first dimension
emis = nan(2*numel(wind_speeds), 3*numel(cities));
tau = nan(2*numel(wind_speeds), 3*numel(cities));

% for each city and wind speed, find the right file and load it. It
% may not exist, then skip.
for a=1:numel(cities)
    % Get the wind vector file, we'll need it for the uncertainty
    % of tau
    wind_path = fullfile(HOMEDIR,'Documents','MATLAB','BEHR','Workspaces','Wind speed');
    wind_name = sprintf('%s-Wind-Conditions-1900UTC-5layers-earthrel.mat',cities{a});
    W = load(fullfile(wind_path, wind_name));
    
    for b=1:numel(wind_speeds)
        xx = ~iscellcontents(strfind(fnames,cities{a}),'isempty') & ~iscellcontents(strfind(fnames,wind_speeds{b}),'isempty');
        if sum(xx) == 0
            continue
        elseif sum(xx) > 1
            E.toomanyfiles(sprintf('%s %s wind file',cities{a},wind_speeds{b}));
        end
        
        SF = load(fullfile(pname, fnames{xx}));
        % find the corresponding line density file; we'll need it
        % for the number of observations
        ldname = regexprep(fnames{xx},'SimpleFits(-ssresid|-unexvar)*','LineDensities');
        LD = load(fullfile(pname, ldname));
        
        for c=1:numel(apri)
            fitparams = struct2array(SF.(apri{c}).ffit)';
            if strcmp(uncert, 'ci')
                fituncert = SF.(apri{c}).stats.percent_ci95 / 100;
            elseif strcmp(uncert, 'sd')
                fituncert = SF.(apri{c}).stats.percentsd / 100;
            else
                E.notimplemented('uncertainty as %s',uncert);
            end
            
            num_obs_fn = sprintf('num_obs_%s',apri{c}(3:end));
            if isfield(LD,num_obs_fn)
                num_obs = LD.(num_obs_fn);
            else
                num_obs = 1;
            end
            fit_uncert = calc_fit_param_uncert(fitparams, fituncert);%, num_obs);
            
            windcrit = str2double(strrep(wind_speeds{b},'pt','.'));
            
            [emis_i, u_emis_i, tau_i, u_tau_i] = compute_emg_emis_tau(fitparams(1), fit_uncert(1), fitparams(2), fit_uncert(2), 'vec', W.windvel(W.windvel>=windcrit));
            i = (a-1)*3 + c;
            if strcmp(uncert,'sd')
                j = (b-1)*3 + 1; % leave an extra line for DoF
            else
                j = (b-1)*2 + 1;
            end
            emis(j,i) = emis_i;
            emis(j+1,i) = u_emis_i;
            tau(j,i) = tau_i;
            tau(j+1,i) = u_tau_i;
            if strcmp(uncert,'sd')
                ld_fn = strrep(apri{c},'f_','no2ld_');
                n_dof = sum(~isnan(LD.(ld_fn)))-5;
                emis(j+2,i) = n_dof;
                tau(j+2,i) = n_dof;
            end
        end
    end
end
end

function [cities, wind_speeds, wind_labels] = get_cities_and_winds_from_fnames(fnames)
if ~iscell(fnames)
    fnames = {fnames};
end
cities = cell(size(fnames));
wind_speeds = cell(size(fnames));
for a=1:numel(fnames)
    [s,e] = regexp(fnames{a},'(Atlanta|Birmingham|Montgomery)');
    cities{a} = fnames{a}(s:e);
    [s,e] = regexp(fnames{a},'\d*pt\d*');
    wind_speeds{a} = fnames{a}(s:e);
end
cities = sort(unique(cities));
wind_speeds = sort(unique(wind_speeds));

wind_labels = cell(numel(wind_speeds)*2,1);
for b=1:numel(wind_speeds)
    wind_labels{(b-1)*2 + 1} = sprintf('$\\geq %s$', strrep(wind_speeds{b},'pt','.'));
end
end

