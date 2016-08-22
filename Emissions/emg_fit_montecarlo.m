function [ x, f0, nacc ] = emg_fit_montecarlo(no2_x, no2_ld, varargin)
%EMG_FIT_MONTECARLO(NO2_X, NO2_LD) Uses Monte Carlo sampling to examine distribution of fitting parameters.
%   One way to calculate the relationships between the fitting parameters
%   in the EMG fit would be to use Monte Carlo sampling to see how
%   frequently combinations of the values for various fitting parameters
%   appear together. This function is essentially a wrapper for the
%   built-in MHSAMPLE function that sets up the probability function and
%   step distribution functions to feed into MHSAMPLE. This duplicates the
%   fitting function in fit_line_density, but returns Inf if the fitting
%   parameters exceed the bounds or violate the constraints specified to
%   fmincon in fit_line_density.
%
%   [X, F0, NACC] = EMG_FIT_MONTECARLO(NO2_X, NO2_LD) - must pass in the
%   x-coordinates and line density values resulting from a line density fit
%   to wind-aligned data (from calc_line_density() function). Outputs are:
%
%       X = the array of points sampled by MHSAMPLE. Each fitting parameter
%       is a column, in order: a, x0, mu_x, sigma_x, B.
%
%       f0 = the initial point chosen randomly by the algorithm.
%
%       nacc = the fraction of times a step is accepted.
%
%   [X, F0, NACC] = EMG_FIT_MONTECARLO(NO2_X, NO2_LD, EMGTYPE) - allows you
%   to select either 'lu' (default) or 'defoy' as the EMG function.  These
%   differ only in the definition of the prefactor, a. 'lu' explicitly
%   includes x0 in the prefactor, while 'defoy' includes it in a.
%
%   MHSAMPLE is set here to get 50,000 points with a burn-in of 1000 points
%   and a thinning of 10 points. The step size was set by trial and error
%   to give a fraction of accepted steps < 0.5, which is recommended to
%   avoid oversampling low-probability areas by ensuring the step is large
%   enough that an uphill move is reasonably unlikely to be allowed.
%
%   Josh Laughner <joshlaugh5@gmail.com> Feb 2016

p=inputParser;
p.addParameter('emgtype','lu');
p.addParameter('f0',[]);
p.parse(varargin{:});
pout=p.Results;
emgtype = pout.emgtype;
f0in = pout.f0;

    lb = [0, 1.6, min(no2_x), 2.5, 0];
    [~,m] = max(no2_ld);
    ub = [Inf, Inf, max(no2_x), no2_x(m) - min(no2_x), max(no2_ld)];
    function emg = emgfxn_defoy(f)
        % Reject bounds
        if any(f < lb) || any(f > ub) || f(2) + f(3) >= max(no2_x) || nonlin_constr(f) > 0
            emg = Inf(size(no2_x));
            return
        end
        emg = f(1)/2 .* exp( (f(4)^2 / (2 * f(2)^2)) - (no2_x - f(3)) ./ f(2) )...
            .* ( 1 - erf( (f(4)^2 - f(2).*(no2_x - f(3)))./(sqrt(2) * f(4) * f(2)) ) ) + f(5);
        emg(isnan(emg)) = Inf;
    end

    function emg = emgfxn_lu(f)
        % Reject bounds
        if any(f < lb) || any(f > ub) || f(2) + f(3) >= max(no2_x) || nonlin_constr(f) > 0
            emg = Inf(size(no2_x));
            return
        end
        emg = f(1)/(2 * f(2)) .* exp( f(3) / f(2) + f(4).^2 / (2*f(2).^2) - no2_x/f(2) )...
        .* erfc( -1/sqrt(2) * ((no2_x-f(3))/f(4) - f(4)/f(2)) ) + f(5);
        emg(isnan(emg)) = Inf;
    end

switch lower(emgtype)
    case 'lu'
        emgfxn = @emgfxn_lu;
    case 'defoy'
        emgfxn = @emgfxn_defoy;
    otherwise
        E.badinput('EMGTYPE must be ''lu'' or ''defoy''')
end

%fitfxn = @(x) nansum((emgfxn(x) - no2_ld).^2) / nansum(no2_ld.^2);
fitfxn = @(x) nansum((emgfxn(x) - no2_ld).^2) / nansum((no2_ld - nanmean(no2_ld)).^2);
%fitfxn = @(x) nansum((emgfxn(x) - no2_ld).^2);

P = @(f) exp(-fitfxn(f));
if isempty(f0in)
    f0 = [rand*9e5+1e5, rand*100, 2*(rand-.5)*50, rand*100, rand*5e3];
else
    f0 = f0in;
end
f0 = max(cat(1,f0,lb),[],1);
f0 = min(cat(1,f0,ub),[],1);
%delta = [1e5, 20, 20, 20, 5e3]; % max step size, found by trial and error to give fraction of accepted steps < 0.5.
delta = [1, 2e-4, 2e-4, 2e-4, 0.5]; % max step size, found by trial and error to give fraction of accepted steps < 0.5.
proprnd = @(f) f + rand(1,5) .* 2 .* delta - delta;

[x,nacc] = mhsample(f0,50000,'pdf',P,'proprnd',proprnd,'symmetric','sym','burnin',1000,'thin',10);

    function [c,ceq] = nonlin_constr(f)
        % Nonlinear constraint requires that the value of the exponential
        % term not exceed 20.  Since the error function is bound to [-1, 1]
        % and 1-erf is always on [0 2], this will hopefully keep the
        % function from going to NaN which happens when the exponetial term
        % goes to infinity while the 1-erf term goes to 0. This constraint
        % is purely to force fmincon to behave, the only possible physical
        % justification is that the a term should control the amount of
        % mass in the plume and not the exponetial itself.
        c = max((f(4)^2 / (2 * f(2)^2)) - (no2_x - f(3)) ./ f(2) - log(20));
        ceq = [];
    end

end

