function [  ] = check_new_fits(  )
%CHECK_NEW_FITS() Compare two versions of EMG fits against each other.
%   Upon running this function, it asks for two directories. First the one
%   with the new version of the fits, then the old version. It iterates
%   through all files present in both directories containing SimpleFits in
%   the name and indicates both the difference in fitting function values
%   (fval) at the optimum solution and if any of the fitting parameters
%   disagree by more than 0.1%. If so it gives the two values as new vs.
%   old.
new_dir = uigetdir('.','Choose the directory with the NEW simple fits');
old_dir = uigetdir('.','Choose the directory with the OLD simple fits');

% Find all the new simple fits 
F = dir(fullfile(new_dir,'*SimpleFits*'));
% For each, find the corresponding old file, load both, and compare the
% ffit parameters, which we'll say should agree within 0.1%.

fitparams = {'a','x_0','mu_x','sigma_x','B'};

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
        for b=1:numel(fitparams)
            new = N.(fns{a}).ffit.(fitparams{b});
            old = O.(fns{a}).ffit.(fitparams{b});
            del = (new - old) * 2 / (new + old); % percent diff vs average value
            if abs(del) > 0.001
                fprintf('\t%s.%s disagree by >0.1%%: %.3g vs %.3g\n', fns{a}, fitparams{b}, new, old);
            end
        end
    end
    fprintf('\n');
end
end

