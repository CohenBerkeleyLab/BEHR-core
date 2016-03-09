function [ Ffixed, F, Nfixed, N ] = fit_line_density_variation( no2_x, no2_ld, n_var, varargin )
%FIT_LINE_DENSITY_VARIATION(NO2_X, NO2_LD, N_VAR) Test the effect on an EMG fit if the fit parameters change by +/- 1 sigma.
%   [Ffixed, F, Nfixed, N] = FIT_LINE_DENSITY_VARIATION(NO2_X, NO2_LD, N_VAR) 
%   takes as input the x-coordinates and line density values (NO2_X and
%   NO2_LD). It will do an initial fit using fit_line_density to get the
%   best values of the parameters (from fmincon) and the standard deviation
%   of those fitting parameters (from nlinfit). It will then fix the value
%   of each fitting parameter in turn within a range of +/- 1 standard
%   deviation from its best value and see how the other parameters and the
%   goodness of fit (defined by the sum of squared residuals) change. N_VAR
%   determines how many values between +/- 1 s.d. each variable takes on.
%   More obviously gives a better idea of the shape of the residuals around
%   the minimum but takes longer. Outputs are:
%
%       Ffixed - an array of structures with dimensions n_var by 5 (the
%       number of parameters). Going across the second dimension varies
%       which variable is being fixed. Going down the first dimension, each
%       structure contains a different fixed value for the parameter fixed
%       in that column.
%
%       F - contains information about the best fit from fmincon.
%
%       Nfixed, N - like Ffixed and F but using nlinfit instead of fmincon
%       to optimize the parameters.
%
%   Note that it is possible to set the fixed parameter to a value outside
%   those allowed by the upper and lower bounds set in fit_line_density.
%   When this happens, that value is essentially discarded, and most values
%   in the Fffixed and Nfixed output structures will be NaNs. The
%   exceptions are those like the fixed value itself that do not need to be
%   returned from fmincon.
%
%   Parameter arguments:
%
%   'emgtype' - sets which EMG function to use. It defaults to 'lu' which
%   includes x0 in the prefactor while the option 'defoy' does not.
%
%   'sd' - allows you to pass in a 1x5 vector of standard deviations
%   corresponding to the five fitting parameters. This can be used to take
%   in standard deviations derived from a Monte Carlo sampling if nlinfit
%   fails for the fit.
%
%   Josh Laughner <joshlaugh5@gmail.com> Feb 2016

p = inputParser;
p.addParameter('emgtype','lu');
p.addParameter('sd',[]);
p.addParameter('fittype','ssresid');
p.parse(varargin{:});
pout = p.Results;

emgtype = pout.emgtype;
fittype = pout.fittype;
sd = pout.sd;

E = JLLErrors;

if ~ismember(emgtype,{'lu','defoy'})
    E.badinput('The parameter ''emgtype'' must be ''lu'' or ''defoy''')
end
if ~ismember(fittype,{'ssresid','unexvar'})
    E.badinput('fittype must be one of ''ssresid'' or ''unexvar''')
end
if ~isempty(sd) && (~isnumeric(sd) || numel(sd) ~= 5)
    E.badinput('If giving standard deviations, they must be given as a 5-element vector')
end

params = {'a','x0','mux','sx','B'};
fns_params = {'a','x_0','mu_x','sigma_x','B'};


    [F.ffit, F.emgfit, F.stats, F.f0, ~, ~, N] = fit_line_density(no2_x, no2_ld,'none','emgtype',emgtype,'fittype',fittype);
    F.ssresid = nansum2((F.emgfit - no2_ld).^2);
    F.no2x = no2_x;
    F.no2ld = no2_ld;

    % These will let us use N in the same plotting function as F
    if numel(fieldnames(N)) > 0
        for a=1:numel(fns_params)
            N.ffit.(fns_params{a}) = N.beta(a);
        end
        N.ssresid = nansum2((N.emg - no2_ld).^2);
        N.no2x = no2_x;
        N.no2ld = no2_ld;
        N.emgfit = N.emg;
    end
    
if isempty(sd)
    % Take the standard deviations from the stats field in F.
    sd = F.stats.sd;
end


% Create the structures that will define both which parameters are fixed
% and take the output from fit_line_density.

Ffixed = struct('fixed_par',params,'fixed_val',[],'ffit',struct('a',nan,'x_0',nan,'mu_x',nan,'sigma_x',nan,'B',nan),'emgfit',nan(size(no2_x)),'stats',struct('sd',0,'ci95',0,'r',0,'r2',0),'f0',nan,'ssresid',nan,'no2x',no2_x,'no2ld',no2_ld);
Ffixed = repmat(Ffixed,n_var,1);
Nfixed = struct('fixed_par',params,'fixed_val',[],'ffit',struct('a',nan,'x_0',nan,'mu_x',nan,'sigma_x',nan,'B',nan),'MSE',nan,'emgfit',nan(size(no2_x)),'ssresid',nan,'no2x',no2_x,'no2ld',no2_ld);
Nfixed = repmat(Nfixed,n_var,1);
for a=1:numel(sd)
    vals = linspace(F.ffit.(fns_params{a})-sd(a),F.ffit.(fns_params{a})+sd(a),n_var);
    for b=1:n_var
        Ffixed(b,a).fixed_val = vals(b);
        Nfixed(b,a).fixed_val = vals(b);
        
        try
            [Ffixed(b,a).ffit, Ffixed(b,a).emgfit, Ffixed(b,a).stats, Ffixed(b,a).f0, ~, ~, Ntmp] = fit_line_density(no2_x, no2_ld, 'none', 'fixed_param',Ffixed(b,a).fixed_par,'fixed_val',Ffixed(b,a).fixed_val, 'emgtype', emgtype,'fittype',fittype);
            Ffixed(b,a).ssresid = nansum2((Ffixed(b,a).emgfit - no2_ld).^2);
        catch err
            if strcmp(err.identifier,'fit_line_density:fixed_val_out_of_range')
                continue
            else
                rethrow(err);
            end
        end
        
        inds = 1:5;
        inds(inds>a) = inds(inds>a) - 1;
        
        if numel(fieldnames(Ntmp)) > 0
            for c=1:numel(fns_params)
                if c == a
                    Nfixed(b,a).ffit.(fns_params{c}) = vals(b);
                else
                    Nfixed(b,a).ffit.(fns_params{c}) = Ntmp.beta(inds(c));
                end
            end
            Nfixed(b,a).MSE = Ntmp.MSE;
            Nfixed(b,a).emgfit = Ntmp.emg;
            Nfixed(b,a).ssresid = nansum2((Nfixed(b,a).emgfit - no2_ld).^2);
        else
            for c=1:numel(fns_params)
                Nfixed(b,a).ffit.(fns_params{c}) = nan;
            end
            Nfixed(b,a).MSE = nan;
            Nfixed(b,a).emgfit = nan(size(no2_x));
            Nfixed(b,a).ssresid = nan;
        end
    end
end



end

