function uncert = calc_fit_param_uncert(fit_params, fit_frac_sds, varargin)
%UNCERT = CALC_FIT_PARAM_UNCERT( FIT_PARAMS, FIT_FRAC_UNCERTS )
% Given the values for the EMG fitting parameters, FIT_PARAMS, and their
% fractional standard deviations, FIT_FRAC_SDS, this will calculate the
% overall uncertainty of the parameters. Note that FIT_FRAC_SDS should
% truly be the fractional uncertainty, that is, if the uncertainty is
% one-half the value of the parameter, it should be represented by 0.5, not
% 50. The return value, UNCERT is an array the same size as FIT_PARAMS and
% FIT_FRAC_SDS that contains the absolute uncertainties for each value.
% FIT_PARAMS can be either the "ffit" structure from fit_line_density.m or
% a vector of those parameters in the same order as the uncertainties. If a
% structure is given, then the uncertainties will be returned in a column
% vector, regardless of the shape of the inputs.  It is recommended that
% FIT_FRAC_SDS be derived from the Hessian of minimizing the unexplained
% variance.
%
% The calculation of uncertainty follows the supplement in Beirle
% et. al. They assign the following uncertainties:
%
%   VCD = 30%, 25% in Lu et al. 2015 
%   NOx:NO2 ratio = 10% Fit
%   confidence interval = ~10-50%
%   SME (standard mean error?) of fit results = ~10-40%
%   Choice of b (across wind integration distance) = 10%
%   Choice of wind fields = 30%
%
% This will ignore the SME as Lu et al 2015 does, because that was for
% Beirle's average over 8 different sectors, rather than the method of
% aligning wind directions. Lu goes on to specify that NOx:NO2 ratio only
% matters for emissions (E ~ a*w/x0 where w is avg. wind speed) but not
% burdens (a) or lifetime (tau = x0 / w) and the uncertainty in VCDs does
% not matter for lifetime. However, I would argue that the uncertainty in
% VCDs DOES matter for lifetime as a spatial bias in the VCDs could
% introduce essentially a second lifetime that would convolve with the true
% lifetime. So this will assume that all the sources of uncertainty count
% for all the fit parameters, except for the NOx:NO2 ratio which is only
% use for emissions (which are not computed here, but instead in
% compute_emg_emis_tau).
%
% This also does some asymmetrical scaling of the standard deviations based
% on what I saw when I went through each parameter for Atlanta and
% Birmingham and looked at how the fit changes as each parameter varied.
% The dependence is definitely asymmetrical. This is also why you really
% should pass in the uncertainty from minimizing the unexplained variance,
% as the uncertainty from minimizing the sum of squared residuals is just
% too small (usually < 0.01%).
%
%UNCERT = CALC_FIT_PARAM_UNCERT( FIT_PARAMS, FIT_FRAC_UNCERTS, NUM_OBS )
% This will use the number of observations given in the array NUM_OBS to
% reduce the error due to VCDs by sqrt(n). Here, n is conservatively taken
% as the fewest number of observations in any single grid cell given by
% NUM_OBS (usually the output NUM_VALID_OBS from CALC_LINE_DENSITY). If not
% given, the value of 25% uncertainty in VCDs will be used unaltered.
%
%CALC_FIT_PARAM_UNCERT( ___ , 'warn', false ) will disable the warning
%about the use of VCD uncertainty of 25% and if the parameter FIT_FRAC_SDS
%is too small.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
p = inputParser;
p.addOptional('num_obs_array',[]);
p.addParameter('warn',true);

p.parse(varargin{:});
pout = p.Results;

num_obs_array = pout.num_obs_array;
warning_bool = pout.warn;

if isstruct(fit_params)
    fit_params = struct2array(fit_params);
    fit_params = fit_params(:); % ensure that it is a column vector
    fit_frac_sds = fit_frac_sds(:);
end

if ~isequal(size(fit_params), size(fit_frac_sds)) || ~isvector(fit_params) || ~isvector(fit_frac_sds)
    E.badinput('FIT_PARAM_VEC and FIT_FRAC_UNCERT_VEC must be vectors of the same size');
end

if ~isnumeric(num_obs_array)
    E.badinput('If given, NUM_OBS_ARRAY must be numeric');
end

if ~isscalar(warning_bool) || (~islogical(warning_bool) && ~isnumeric(warning_bool))
    E.badinput('The value for parameter WARN must be a scalar value that can be interpreted as a logical (either a logical or numeric value)');
end

if warning_bool && any(fit_frac_sds < 1e-3)
    warning('Fractional SDs of < 0.1% detected, it is not recommended that the SDs from minimizing sum of squared residuals are used')
end

% Make input into row vectors and replicate them so that element-wise
% products with the scaling factors (defined as rows in a 2x5 matrix) work
if ~isrow(fit_params); fit_params = fit_params'; end
fit_params = repmat(fit_params, 2, 1);
if ~isrow(fit_frac_sds); fit_frac_sds = fit_frac_sds'; end
fit_frac_sds = repmat(fit_frac_sds, 2, 1);


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define the scaling factors for the fitting uncertainties. Bottom side in
% first row, top side in second.
sd_scales = emg_sd_scaling_factors;
         

if isempty(num_obs_array);
    s_vcd = 0.25;
    if warning_bool
        warning('No number of valid observations given in the line density structure, setting VCD uncertainty to 25%')
    end
else
    % If there is a number of valid observations variable, then we
    % will reduce the uncertainty in the line density by the
    % smallest number present to be conservative.
    s_vcd = 0.25 / sqrt(min(num_obs_array(num_obs_array>0)));
end
s_b = 0.1;
s_wind = 0.3;
per_uncert = sqrt((fit_frac_sds .* sd_scales).^2 + s_vcd.^2 + s_b.^2 + s_wind.^2);
uncert = abs(fit_params) .* per_uncert;
end