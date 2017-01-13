function [  ] = run_many_emg_fits(  )
% Do all the EMG fitting steps specified for all the files specified
load_dir = '/Users/Josh/Documents/MATLAB/BEHR/Workspaces/EMG fits/Autorun/FullDaily-NumObs/50km-side-earthrel-no0to-112-lonwt14-1822UTC';
F = dir(fullfile(load_dir,'*LineDensities*.mat'));
fit_steps = {'simple'}; % simple, var, or both
fitfxn = 'ssresid'; % ssresid or unexvar
tic
for a=1:numel(F)
    D = load(fullfile(load_dir,F(a).name));
    if ismember('simple',fit_steps)
        O = simple_fitting(D, fitfxn);
        savefields(O, load_dir, make_simplefit_name(F(a).name, fitfxn));
    end
    if ismember('var',fit_steps)
        if ~isempty(strfind(F(a).name,'Montgomery'))
            % Montgomery isn't large enough to do this procedure reliably
            continue
        end
        O = var_fitting(D, fitfxn);
        savefields(O, load_dir, make_varfit_name(F(a).name, fitfxn));
    end
end
toc
end

function savefields(S, save_dir, save_name)
fns = fieldnames(S);
for f=1:numel(fns)
    eval(sprintf('%1$s = S.%1$s;',fns{f}));
end
save(fullfile(save_dir,save_name), fns{:});
end

function S = simple_fitting(In, fitfxn)
fmincon_out_level = 'none'; % could be 'none', 'final', 'iter'
[S.f_hyfast.ffit, S.f_hyfast.emgfit, S.f_hyfast.stats, S.f_hyfast.f0, ~, S.f_hyfast.fitresults] = fit_line_density(In.no2x_hyfast, In.no2ld_hyfast, fmincon_out_level,'fittype',fitfxn);
[S.f_hyslow.ffit, S.f_hyslow.emgfit, S.f_hyslow.stats, S.f_hyslow.f0, ~, S.f_hyslow.fitresults] = fit_line_density(In.no2x_hyslow, In.no2ld_hyslow, fmincon_out_level,'fittype',fitfxn);
[S.f_mnfast.ffit, S.f_mnfast.emgfit, S.f_mnfast.stats, S.f_mnfast.f0, ~, S.f_mnfast.fitresults] = fit_line_density(In.no2x_mnfast, In.no2ld_mnfast, fmincon_out_level,'fittype',fitfxn);
[S.f_mnslow.ffit, S.f_mnslow.emgfit, S.f_mnslow.stats, S.f_mnslow.f0, ~, S.f_mnslow.fitresults] = fit_line_density(In.no2x_mnslow, In.no2ld_mnslow, fmincon_out_level,'fittype',fitfxn);
[S.f_mn108fast.ffit, S.f_mn108fast.emgfit, S.f_mn108fast.stats, S.f_mn108fast.f0, ~, S.f_mn108fast.fitresults] = fit_line_density(In.no2x_mn108fast, In.no2ld_mn108fast, fmincon_out_level,'fittype',fitfxn);
[S.f_mn108slow.ffit, S.f_mn108slow.emgfit, S.f_mn108slow.stats, S.f_mn108slow.f0, ~, S.f_mn108slow.fitresults] = fit_line_density(In.no2x_mn108slow, In.no2ld_mn108slow, fmincon_out_level,'fittype',fitfxn);
if all(isfield(In,{'no2x_wrffast','no2x_wrfslow'}))
    [S.f_wrffast.ffit, S.f_wrffast.emgfit, S.f_wrffast.stats, S.f_wrffast.f0, ~, S.f_wrffast.fitresults] = fit_line_density(In.no2x_wrffast, In.no2ld_wrffast, fmincon_out_level,'fittype',fitfxn);
    [S.f_wrfslow.ffit, S.f_wrfslow.emgfit, S.f_wrfslow.stats, S.f_wrfslow.f0, ~, S.f_wrfslow.fitresults] = fit_line_density(In.no2x_wrfslow, In.no2ld_wrfslow, fmincon_out_level,'fittype',fitfxn);
end
end

function save_name = make_simplefit_name(ld_name, fitfxn)
save_name = strrep(ld_name,'LineDensities',['SimpleFits-',fitfxn]);
end

function S = var_fitting(In, fitfxn)
fprintf('Sampling the residuals using the standard deviations derived from the Hessian\n');
[S.Ffix_hyfast, S.F_hyfast] = fit_line_density_variation(In.no2x_hyfast, In.no2ld_hyfast, 30,'fittype',fitfxn);
[S.Ffix_hyslow, S.F_hyslow] = fit_line_density_variation(In.no2x_hyslow, In.no2ld_hyslow, 30,'fittype',fitfxn);
[S.Ffix_mnfast, S.F_mnfast] = fit_line_density_variation(In.no2x_mnfast, In.no2ld_mnfast, 30,'fittype',fitfxn);
[S.Ffix_mnslow, S.F_mnslow] = fit_line_density_variation(In.no2x_mnslow, In.no2ld_mnslow, 30,'fittype',fitfxn);
[S.Ffix_mn108fast, S.F_mn108fast] = fit_line_density_variation(In.no2x_mn108fast, In.no2ld_mn108fast, 30,'fittype',fitfxn);
[S.Ffix_mn108slow, S.F_mn108slow] = fit_line_density_variation(In.no2x_mn108slow, In.no2ld_mn108slow, 30,'fittype',fitfxn);
end

function save_name = make_varfit_name(ld_name, fitfxn)
save_name = strrep(ld_name,'LineDensities',['VariationalFits-',fitfxn]);
end
