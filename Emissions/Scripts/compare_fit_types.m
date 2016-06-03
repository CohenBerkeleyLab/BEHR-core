function [  ] = compare_fit_types(  )
%COMPARE_FIT_TYPES() Checks if the ssresid and unexvar fits differ
%   On running this function, it will open a dialogue box asking you to
%   choose a directory with simple EMG fits in it. For each ssresid file it
%   finds, it looks for the corresponding unexvar file and checks if any of
%   the fitting parameters disagree by more than 0.1%, if so, it prints
%   that information.

fit_dir = uigetdir('.','Choose the directory with the simple fits');

% Find all the new simple fits 
F = dir(fullfile(fit_dir,'*SimpleFits-ssresid*'));
% For each, find other fit files
Funex = dir(fullfile(fit_dir,'*SimpleFits-unexvar*'));


fitparams = {'a','x_0','mu_x','sigma_x','B'};

for f=1:numel(F)
    if ~exist(fullfile(fit_dir,Funex(f).name),'file')
        fprintf('No unexvar file for %s\n', F(f).name);
        continue
    else
        fprintf('%s:\n',F(f).name);
    end
    
    N = load(fullfile(fit_dir, F(f).name));
    O = load(fullfile(fit_dir, Funex(f).name));
    
    fns = fieldnames(N);
    for a=1:numel(fns)
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

