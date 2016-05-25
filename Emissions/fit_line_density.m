function [ ffit, emgfit, param_stats, f0, history, fitresults, N, L ] = fit_line_density( no2_x, no2_ld, varargin )
%FIT_LINE_DENSITY Fits an exponentiall modified Gaussian function to a line density
%   Expoentially modified Gaussian (EMG) functions work very well to fit
%   line densities of NO2 plumes observed from space because the two
%   components allow it to capture the plume build up or at slow winds
%   (Gaussian) and its decay downwind at faster winds (exponential decay).
%   c.f. de Foy et al., Atmos. Environ., 2014, pp. 66-77.
%
%   This function requires as input the x-coordinates in kilometers and
%   line densities (moles / km recommended, though other units *should*
%   work). A third optional input controls the level of console output
%   fmincon provides: it defaults to 'iter' meaning that fmincon will
%   output information about each step. Another option is 'final' which
%   will only output a final report.
%
%   This function will attempt to fit an EMG function to the line density
%   using fmincon. The five fitting terms (a, x_0, mu_x, sigma_x, and B) as
%   described in de Foy 2014 are returned as the first five outputs. The
%   sixth output is the EMG fit itself as a vector that will have the same
%   x-coordinates as input.
%
%   There are several additional outputs useful for debugging poor fits.
%   The seventh output (f0) is the initial guess of the five fitting
%   parameters. The eighth output (history) will contain a record of how
%   the fitting parameters were iterated during the minimization. The ninth
%   output contains additional information returned from fmincon.
%
%   There are several parameters available: 
%
%       'fixed_param' - a string that sets which fitting parameter should
%       be fixed. Must be used with 'fixed_val' to set the value. String
%       must be one of 'a', 'x0', 'mux', 'sx', or 'B'.
%
%       'fixed_val' - sets the value of the parameter fixed by fixed_param.
%
%       'f0' - allows you to define the initial point for the fit
%       minimization. A 5-element vector with the values in the order of a,
%       x0, mu_x, sigma_x, B. If not set, the function will make a best
%       guess as to a reasonable initial point.
%
%       'lb' and 'ub' - allow you to specify the lower and upper bounds of
%       each of the five fitting parameters again as 5 element vectors in
%       the same order as f0. If not given, these will be set within this
%       function to physically realistic values.
%
%       'emgtype' - allows you to vary which one of 3 EMG functions are
%       used. There were slight differences between de Foy 2014 and Lu 2015
%       regarding the form of the function; 'defoy' and 'lu' respectively
%       select each paper's respective function. The only practical
%       difference is in the definition of a: a_defoy = a_lu / x0. Defaults
%       to 'lu.'
%
%       'fittype' - allows you to decide the normalization of the fit
%       function. 'ssresid' (default) makes the fit function simply the sum
%       of squared residuals between the line density and EMG fit:
%
%           f(a,x0,mux,sigma,B) = nansum( (no2_ld - emgfit(no2_x) ).^2 )
%
%       The other option is 'unexvar' which stands for 'unexplained
%       variance.' It normalizes the sum of squared residuals by the
%       variance in the line density:
%
%           f_unexvar() = f / nansum( (no2_ld - nanmean(no2_ld)).^2 )
%
%       These should not give different answers for the actual fitting
%       parameters, but they do affect the uncertainty, as naturally the
%       unexplained variance has much shallower minima, and so the Hessian
%       matrix (which the uncertainties are derived from) reflects the
%       weaker curvature.
%
%       'nattempts' - defaults to 10, how many times the fitting function
%       should run fmincon. This function will always use the best guess
%       for the initial point first, then randomize f0 nattempts - 1 times.
%       The randomization is such that it will uniformly select a value
%       between the upper and lower bounds unless one bound is infinity,
%       then it will use a bound around 10x larger than the typical value
%       for that parameter.
%
%       DEBUG_LEVEL - integer controlling the amount of messages printed to
%       the console. Set to 0 to disable.
%
%   Josh Laughner <joshlaugh5@gmail.com> 5 Feb 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addOptional('fmincon_output','iter',@(x) ismember(x,{'none','iter','final'}));
p.addParameter('fixed_param','',@ischar);
p.addParameter('fixed_val',[]);
p.addParameter('f0',[]);
p.addParameter('lb',[]);
p.addParameter('ub',[]);
p.addParameter('emgtype','lu');
p.addParameter('fittype','ssresid');
p.addParameter('nattempts',10);
p.addParameter('DEBUG_LEVEL',1);
p.parse(varargin{:});
pout=p.Results;

fmincon_output = pout.fmincon_output;
fixed_param = pout.fixed_param;
fixed_val = pout.fixed_val;
f0in = pout.f0;
ubin = pout.ub;
lbin = pout.lb;
emgtype = lower(pout.emgtype);
fittype = lower(pout.fittype);
nattempts = pout.nattempts;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

E = JLLErrors;
E.addCustomError('fixed_val_out_of_range','The fixed value given is outside the upper and lower bounds allowed for that value.');

if ~isvector(no2_x) || ~isnumeric(no2_x)
    E.badinput('no2_x must be a vector of numeric inputs')
end
if ~isvector(no2_ld) || ~isnumeric(no2_ld) || any(size(no2_ld) ~= size(no2_x))
    E.badinput('no2_ld must be a numeric vector the same size as no2_ld')
end
if ~isempty(fixed_param) && ~ismember(fixed_param,{'a','x0','mux','sx','B'})
    E.badinput('fixed_param (if given) must be one of a, x0, mux, sx, or B')
end
if ~isempty(fixed_val) && ~ischar(fixed_val) && ~isscalar(fixed_val)
    E.badinput('fixed_val must be a scalar number')
end
if xor(isempty(fixed_param),isempty(fixed_val))
    E.badinput('Both or neither of fixed_param and fixed_val must be set')
end
if ~isempty(f0in) && numel(f0in) ~= 5
    E.badinput('If an f0 is specified as input, it must have five elements')
end
if ~isempty(ubin) && numel(ubin) ~= 5
    E.badinput('If an ub is specified as input, it must have five elements')
end
if ~isempty(lbin) && numel(lbin) ~= 5
    E.badinput('If an lb is specified as input, it must have five elements')
end
if ~ismember(emgtype,{'defoy','lu'})
    E.badinput('The parameter ''emgtype'' must be one of defoy or lu')
end
if ~ismember(fittype,{'ssresid','unexvar'})
    E.badinput('The parameter ''fittype'' must be one of ssresid or unexvar')
end
if ~isscalar(nattempts) || ~isnumeric(nattempts) || nattempts < 1
    E.badinput('nattempts must be a scalar >= 1');
end
if ~isscalar(DEBUG_LEVEL) || ~isnumeric(DEBUG_LEVEL) || DEBUG_LEVEL < 0
    E.badinput('DEBUG_LEVEL must be a scalar >= 0');
end

% Define the fit function, which although is physically a function of x, is
% to be minimized by the variation of five parameters: a, x_0, sigma_x,
% mu_x, and B. This function is the sum of squared residuals between the
% EMG fit and the actual line density. The five variational parameters will
% be set as the five elements of the vector f, necessary for inserting them
% into the fmincon routine.

    function e = emgfxn_defoy(f,x)
%emgfxn = @(f,x) 
        e = f(1)/2 .* exp( (f(4)^2 / (2 * f(2)^2)) - (x - f(3)) ./ f(2) )...
        .* ( 1 - erf( (f(4)^2 - f(2).*(x - f(3)))./(sqrt(2) * f(4) * f(2)) ) ) + f(5);
        e(isnan(e)) = Inf;
    end

    function e = emgfxn_lu(f,x)
        e = f(1)/(2 * f(2)) .* exp( f(3) / f(2) + f(4).^2 / (2*f(2).^2) - x/f(2) )...
        .* erfc( -1/sqrt(2) * ((x-f(3))/f(4) - f(4)/f(2)) ) + f(5);
        e(isnan(e)) = Inf;
    end

switch emgtype
    case 'defoy'
        emgfxn = @emgfxn_defoy;
    case 'lu'
        emgfxn = @emgfxn_lu;
end

if strcmp(fittype,'ssresid')
    fitfxn = @(f) nansum((emgfxn(f,no2_x) - no2_ld).^2);
elseif strcmp(fittype,'unexvar')
    fitfxn = @(f) nansum((emgfxn(f,no2_x) - no2_ld).^2) / nansum((no2_ld - nanmean(no2_ld)).^2);
else
    E.notimplemented(sprintf('fittype = %s',fittype));
end

history.x = [];
opts = optimoptions('fmincon','Display',fmincon_output,'OutputFcn',@outfun);%,'MaxFunEvals',10000,'MaxIter',5000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Best guess of initial conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(f0in)
    f0 = f0in;
else
    f0 = ones(1,5);

    % f(2) is x_0 or length scale for chemical decay. Assume a lifetime of 3 hr
    % and 5 m/s winds = 54 km for 1 lifetime
    f0(2) = 54;
    % f(1) is a or the mass of NO2 in one length scale. Tried considering this
    % the integral of line density from ~0 to one x_0 downwind, but that tends
    % to be too large by a factor of 10. 
    %      i = find(no2_x > 0, 1, 'first');
    %       j = find(no2_x > no2_x(i) + f0(2), 1, 'first');
    %       f0(1) = trapz(no2_x(i:j), no2_ld(i:j));
    % Instead, let's just try starting out as half the maximum value of the
    % line density.
    f0(1) = max(no2_ld)/2;

    % f(3) is mu_x or where the average for the Gaussian is, this should be
    % around the maximum.
    [~,m] = max(no2_ld);
    f0(3) = no2_x(m);

    % f(4) is sigma_x, the std. dev. of the gaussian. Can be defined as
    % FWHM/2.355.
    halfmax = (max(no2_ld) - no2_ld(1))/2 + no2_ld(1);
    mpre = find(no2_ld(1:m) < halfmax, 1, 'last');
    fwhm = abs(interp1(no2_ld(mpre:mpre+1), no2_x(mpre:mpre+1), halfmax));
    f0(4) = fwhm / 2.355;

    % f(5) is the constant term B, or background. So let's just set to be the
    % left-most line density.
    f0(5) = min(no2_ld);
end
%%%%%%%%%%%%%%%%%
% Set up bounds %
%%%%%%%%%%%%%%%%%
% Remember, f = [a, x0, mu_x, sigma_x, B]
f_lb = nan(1,5);
f_ub = nan(1,5);

% a (related to plume mass) must be > 0 to be physically meaningful
% 100 mol/km is a little more than 
f_lb(1) = 100; f_ub(1) = Inf; %f_ub(1) = max(no2_ld)*1.5;

% x0 (length scale of chemical decay) must be > 0 to be physically
% meaningful. This term tends to cause problems in the fit when the
% algorithm thinks it sees a minimum where x0 gets very small, so I chose
% one-third of the grid size as the minimum. This assumes that after 3
% lifetimes the line density would be within 5% of background, so if we
% consider that to be essentially equal to background, the smallest
% lifetime distance distinguishable would be 1/3rd of the distance between
% grid cells, as any shorter lifetime would also result in complete decay
% to background at smaller distances that the resolution of the data. This
% basically says that the decay must be observable at the resolution of the
% data. Realistically, the lifetime distances should never even be this
% short, or if they are this short, then then wind is so slow that the
% Gaussian character should take over.
f_lb(2) = 1.6; f_ub(2) = Inf;

% mu_x should be close to the position of the source, but what we can say
% for certain is that it must lie somewhere on no2_x. A stronger condition
% might be that it must lie within 1 std. dev. of the maximum of the line
% density, that may be useful if fitting continues to be foolish.
f_lb(3) = min(no2_x); f_ub(3) = max(no2_x);
%f_lb(3) = -50; f_ub(3) = 50;

% sigma_x describes the width of the Gaussian; it must be positive (and not
% just technically positive - it should have at least some width. Therefore
% similarly to x0, we will say that the narrowest observable Gaussian is
% one that consists of only 3 data points, i.e. the half width at the base
% is equal to the grid resolution of 0.05 deg ~ 5 km. The standard
% deviation is half that width (4 sd = full width at base) so we will
% require the value of sigma to be greater than half the grid resolution.
% We will also assume that the Gaussian build up is wholly contained in the
% domain, and so its maximum width is the distance from the left side of
% the domain to the point of maximum NO2.
[~,m] = max(no2_ld);
f_lb(4) = 2.5; f_ub(4) = no2_x(m) - min(no2_x);

% B is the background value. It must be > 0 to be physically meaningful,
% and should really not exceed the maximum line density.
f_lb(5) = 0; f_ub(5) = max(no2_ld);

% Also constrain that mu_x + x0 must fall within the domain, i.e. assume
% that at least one chemical lifetime occurs within the distance studied.
A = [0 1 1 0 0];
b = max(no2_x);

% Overwrite the upper and lower bounds with user input if given. Doing it
% here just lets me keep the section defining the limits organized how it
% is, which is much easier to read, and allow the user to only specify
% upper or lower limits (instead of both) if they desire.
if ~isempty(ubin)
    f_ub = ubin;
end
if ~isempty(lbin)
    f_lb = lbin;
end

if ischar(fixed_val) && strcmpi(fixed_val, 'f0')
    switch fixed_param
        case 'a'
            fixed_val = f0(1);
        case 'x0'
            fixed_val = f0(2);
        case 'mux'
            fixed_val = f0(3);
        case 'sx'
            fixed_val = f0(4);
        case 'B'
            fixed_val = f0(5);
    end
elseif ischar(fixed_val)
    E.badinput('fixed_val should be a number if not the string ''f0''')
end

switch fixed_param
    case 'a'
        if fixed_val < f_lb(1) || fixed_val > f_ub(1)
            E.callCustomError('fixed_val_out_of_range');
        end
        inds=2:5;
        nlcon = @(f) nonlin_constr([fixed_val, f]);
        fitfxn_fix = @(f) fitfxn([fixed_val, f]);
        emgfxn_fix = @(f,x) emgfxn([fixed_val, f],x);
    case 'x0'
        if fixed_val < f_lb(2) || fixed_val > f_ub(2)
            E.callCustomError('fixed_val_out_of_range');
        end
        inds=[1,3,4,5];
        nlcon = @(f) nonlin_constr([f(1), fixed_val, f(2:4)]);
        fitfxn_fix = @(f) fitfxn([f(1), fixed_val, f(2:4)]);
        emgfxn_fix = @(f,x) emgfxn([f(1), fixed_val, f(2:4)],x);
    case 'mux'
        if fixed_val < f_lb(3) || fixed_val > f_ub(3)
            E.callCustomError('fixed_val_out_of_range');
        end
        inds=[1,2,4,5];
        nlcon = @(f) nonlin_constr([f(1:2), fixed_val, f(3:4)]);
        fitfxn_fix = @(f) fitfxn([f(1:2), fixed_val, f(3:4)]);
        emgfxn_fix = @(f,x) emgfxn([f(1:2), fixed_val, f(3:4)],x);
    case 'sx'
        if fixed_val < f_lb(4) || fixed_val > f_ub(4)
            E.callCustomError('fixed_val_out_of_range');
        end
        inds=[1,2,3,5];
        nlcon = @(f) nonlin_constr([f(1:3), fixed_val, f(4)]);
        fitfxn_fix = @(f) fitfxn([f(1:3), fixed_val, f(4)]);
        emgfxn_fix = @(f,x) emgfxn([f(1:3), fixed_val, f(4)],x);
    case 'B'
        if fixed_val < f_lb(5) || fixed_val > f_ub(5)
            E.callCustomError('fixed_val_out_of_range');
        end
        inds=1:4;
        nlcon = @(f) nonlin_constr([f, fixed_val]);
        fitfxn_fix = @(f) fitfxn([f, fixed_val]);
        emgfxn_fix = @(f,x) emgfxn([f, fixed_val],x);
    otherwise
        inds = 1:5;
        nlcon = @nonlin_constr;
        fitfxn_fix = @(f) fitfxn(f);
        emgfxn_fix = @(f,x) emgfxn(f,x);
end
f0 = f0(inds);
f_lb = f_lb(inds);
f_ub = f_ub(inds);
A = A(inds);

err_attempts = 0;
attempts = 0;
best_fval = inf;
while attempts < nattempts
    try
        % Always try with our best guess first, then let it randomize
        % nattempts-1 times to see if we find a better minimum.
        [this_fitparams, this_fitresults.fval, this_fitresults.exitFlag, this_fitresults.output, this_fitresults.lambda, this_fitresults.grad, this_fitresults.Hessian]...
            = fmincon(fitfxn_fix, f0, A, b, [], [], f_lb, f_ub, nlcon, opts);
        if this_fitresults.fval < 0 || abs(imag(this_fitresults.fval))>0
            E.callError('fval_lt0','The value of the fitting function is negative or imaginary. Double check the form of the fitting function.');
        elseif this_fitresults.fval < best_fval
            best_fval = this_fitresults.fval;
            fitparams = this_fitparams;
            fitresults = this_fitresults;
        end
    catch err
        if strcmp(err.identifier,'optim:barrier:UsrObjUndefAtX0')
            if err_attempts < 10
                fprintf('fmincon input function undefined at f0; randomizing f0 (%s) to try again\n',mat2str(f0))
                attempts = attempts - 1; % don't count this as one of the regular attempts.
                err_attempts = err_attempts+1;
            else
                E.callError('fmincon_failure','fmincon has been undefined at f0 10 times using fixed parameter = %s, fixed value = %.3g, aborting to prevent infinite loop',fixed_param,fixed_val)
            end
        else
            rethrow(err);
        end
    end
    f0 = randomize_f0(f_lb, f_ub, inds);
    attempts = attempts + 1;
end

if isinf(best_fval)
    E.callError('fmincon_failure','fmincon never found a valid fit')
end

unc_opts = optimoptions('fminunc','algorithm','quasi-newton','maxfunevals',1,'display','none'); % will switch to QN anyway to get Hessian w/o gradient; this avoids that message.
[f_unc,~,~,~,~,unc_hessian] = fminunc(fitfxn_fix, fitparams, unc_opts);

if any(abs(f_unc - fitparams) > 0.01)
    warning('fminunc found a different minimum than fmincon, check it''s Hessian before using the standard deviations')
end

emgfit = emgfxn_fix(fitparams, no2_x);
ffit = unfix_params(fitparams, fixed_param, fixed_val);

% The inverse of the Hessian is an approximation to the covariance matrix
% (Dovi 1991, Appl. Math. Lett. Vol 4, No 1, pp. 87-90; also
% http://www.mathworks.com/matlabcentral/answers/153414-estimator-standard-errors-using-fmincon-portfolio-optimization-context)
% however the fmincon Hessian can be inaccurate, hence we take it from
% fminunc that is initialized at the solution from fmincon.
param_stats.sd = sqrt(diag(inv(unc_hessian)));
param_stats.percentsd = param_stats.sd ./ abs(fitparams') * 100;

% tinv gives the t value for a one-tailed distribution, we want two-tailed
% so alpha must be halved, thus 97.5% certainty one-tailed is equivalent to
% 95% two-tailed. Also, since we're dealing with a fit, we've used up two
% degrees of freedom (slope & intercept).
student_t = tinv(0.975,numel(no2_x)-2);
param_stats.ci95 = param_stats.sd .* student_t ./ sqrt(numel(no2_x)-2);
param_stats.percent_ci95 = param_stats.ci95 ./ abs(fitparams') * 100;
fitresults.fminunc_soln = f_unc;
fitresults.fminunc_hessian = unc_hessian;

% Calculate r (http://onlinestatbook.com/2/describing_bivariate_data/calculation.html)
x = no2_ld - nanmean(no2_ld);
y = emgfit - nanmean(emgfit);
param_stats.r = nansum2(x .* y) ./ sqrt(nansum2(x.^2) .* nansum2(y.^2));
% This is really a check I did r right more than anything.
param_stats.r2 = 1 - nansum2((no2_ld - emgfit).^2) / nansum2((no2_ld - nanmean(no2_ld)).^2);



if nargout >= 7
    opts = statset('funvalcheck','off'); % This prevents nlinfit from stopping if NaNs or Infs occur. But can make the CovB matrix unusually large.
    %[N.beta, N.R, N.J, N.CovB, N.MSE, N.ErrorModelInfo] = nlinfit(no2_x,no2_ld,emgfxn_fix,f0,opts);
    [N.beta, N.R, N.J, N.CovB, N.MSE, N.ErrorModelInfo] = nlinfit(no2_x,no2_ld,emgfxn_fix,fitparams,opts);
    N.emg = emgfxn_fix(N.beta,no2_x);
    N.ffit = unfix_params(N.beta, fixed_param, fixed_val);
end

if nargout >= 8
    [L.beta,~,~,~,~,~,L.Jacobian] = lsqcurvefit(emgfxn_fix, fitparams, no2_x, no2_ld, f_lb, f_ub);
    L.Cov = inv(L.Jacobian' * L.Jacobian);
    L.emg = emgfxn_fix(L.beta, no2_x);
    L.ffit = unfix_params(L.beta, fixed_param, fixed_val);
end

    function stop = outfun(x,optimvals,state)
        stop = false;
        switch state
            case 'iter'
                history.x = [history.x; x];
        end
    end

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

function ffit = unfix_params(fitparams, fixed_param, fixed_val)

switch fixed_param
    case 'a'
        ffinal = [fixed_val, fitparams];
    case 'x0'
        ffinal = [fitparams(1), fixed_val, fitparams(2:4)];
    case 'mux'
        ffinal = [fitparams(1:2), fixed_val, fitparams(3:4)];
    case 'sx'
        ffinal = [fitparams(1:3), fixed_val, fitparams(4)];
    case 'B'
        ffinal = [fitparams, fixed_val];
    otherwise 
        ffinal = fitparams;
end
ffit.a = ffinal(1);
ffit.x_0 = ffinal(2);
ffit.mu_x = ffinal(3);
ffit.sigma_x = ffinal(4);
ffit.B = ffinal(5);
end

function f0 = randomize_f0(f_lb, f_ub, inds)
% Simplest method to randomize f0 is to allow it to take on any value
% between f_lb and f_ub for fully bounded values.  For those with an
% infinite limit, use some senisibly large but finite bound.
% Reasonable upper bound for each parameter if set to infinity:
ub_lim = [5e4, 1e3, 1e3, 1e3, 5e4];
ub_lim = ub_lim(inds);
lb_lim = [0 0 -1e3 0 0];
lb_lim = lb_lim(inds);
% Must remove the fixed parameter if one exists, hence the use of inds.

rvec = rand(size(f_ub));
uu = isinf(f_ub);
f_ub(uu) = ub_lim(uu);
ll = isinf(f_lb);
f_lb(ll) = lb_lim(ll);

f0 = f_lb + rvec .* (f_ub - f_lb);
end
