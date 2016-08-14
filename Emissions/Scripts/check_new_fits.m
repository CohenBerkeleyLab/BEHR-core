function [  ] = check_new_fits(  )
%CHECK_NEW_FITS() Compare two versions of EMG fits against each other.
%   Upon running this function, it asks for two directories. First the one
%   with the new version of the fits, then the old version. It iterates
%   through all files present in both directories containing SimpleFits in
%   the name and indicates both the difference in fitting function values
%   (fval) at the optimum solution and if any of the fitting parameters
%   disagree by more than 0.1%. If so it gives the two values as new vs.
%   old.
always_print_params = false;
always_print_cis = false;

new_dir = uigetdir('.','Choose the directory with the NEW simple fits');
old_dir = uigetdir('.','Choose the directory with the OLD simple fits');

% Find all the new simple fits 
F = dir(fullfile(new_dir,'*SimpleFits*'));
% For each, find the corresponding old file, load both, and compare the
% ffit parameters, which we'll say should agree within 0.1%.

fitparams = {'a','x_0','mu_x','sigma_x','B'};

frac_diff_crit = 0.001; % the value that the fractional difference must be greater than to be considered different

for f=1:numel(F)
    if ~exist(fullfile(old_dir,F(f).name),'file')
        fprintf('No old file for %s\n', F(f).name);
        continue
    else
        fprintf('%s:\n',F(f).name);
    end
    
    N = load(fullfile(new_dir, F(f).name));
    O = load(fullfile(old_dir, F(f).name));
    
    fns = fieldnames(N);
    for a=1:numel(fns)
        % Compare the fval values
        newfval = N.(fns{a}).fitresults.fval;
        oldfval = O.(fns{a}).fitresults.fval;
        fprintf('%s: New fval = %.3g, old fval = %.3g\n',fns{a},newfval,oldfval);
        
        % Compare the actual fit parameter values
        for b=1:numel(fitparams)
            new = N.(fns{a}).ffit.(fitparams{b});
            old = O.(fns{a}).ffit.(fitparams{b});
            del = (new - old) * 2 / (new + old); % percent diff vs average value
            if abs(del) > frac_diff_crit || always_print_params
                if abs(del) > frac_diff_crit
                    agree_str = sprintf('disagree by >%g%%',frac_diff_crit*100);
                else
                    agree_str = sprintf('agree to within %g%%', frac_diff_crit*100);
                end
                fprintf('\t%s.%s %s: new = %.3g vs old = %.3g\n', fns{a}, fitparams{b}, agree_str, new, old);
            end
        end
        
        % Also check the uncertainties
        new = N.(fns{a}).stats.ci95;
        old = O.(fns{a}).stats.ci95;
        del = (new - old) * 2 ./ (new + old);
        xx = abs(del) > frac_diff_crit;
        if always_print_cis
            fprintf('\n\t95%% C.I.s percent difference = %s%%\n', mat2str(del*100));
        elseif any(xx)
            badparams = strjoin(fitparams(xx), ', ');
            fprintf('\n\t95%% C.I.s for %s in %s disagree by >%g%%\n', badparams, fns{a}, frac_diff_crit*100);
        end
    end
    fprintf('\n');
end
end

