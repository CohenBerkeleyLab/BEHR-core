function [ emg_lu, emg_defoy, emg_conv, emg_math ] = convolution_test( ffit, x )
%CONVOLUTION_TEST Test implementation of the MATLAB conv function for EMG functions
%   Beirle et al. 2011 actually defines their EMG function as the
%   convolution of a truncated exponential and Gaussian function:
%       e(x) = exp( -(x-mu_x)/x_0 ) for x >= mu_x
%              0                    for x < 0
%       G(x) = 1/(sqrt(2 pi) sigma_x) * exp( -x^2 / (2 sigma_x^2))
%
%       F(x) = E (e (*) G)(x) + B where (*) is the convolution operator
%       (just * is multiplication).
%
%   Matlab's convolution operator acts on vectors, therefore I need to test
%   if my implementation of it is correct for this case, where we have a
%   known convolution result function before I can use it for Liu et al.
%   2016, which uses a convolution function to convolve the truncated
%   exponential with an empirically derived function representing emissions
%   sources.
%
%   Give this a structure, ffit, with the five fitting parameters as fields
%   (a, x_0, mu_x, sigma_x, B) and the x coordinates as x and it will
%   generate EMG functions from the form given in Lu et al 2015, Eq. (4)
%   and the convolution in Eq. (1).

    function e = emgfxn_lu(x, a, x_0, mu_x, sigma_x, B)
        e = a/(2 * x_0) .* exp( mu_x / x_0 + sigma_x.^2 / (2*x_0.^2) - x/x_0 )...
        .* erfc( -1/sqrt(2) * ((x-mu_x)/sigma_x - sigma_x/x_0) ) + B;
        e(isnan(e)) = Inf;
    end

    function e = emgfxn_defoy(x, a, x_0, mu_x, sigma_x, B)
        e = a/(2 .* x_0) .* exp( (sigma_x .^ 2) ./ (2 .* x_0 .^ 2) - (x - mu_x) ./ x_0 )...
            .* (1 - erf( (sigma_x .^ 2 - x_0 .* (x - mu_x)) ./ (sqrt(2) .* sigma_x .* x_0))) + B;
        e(isnan(e)) = Inf;
    end

    function e = trunc_x(x, x_0, mu_x)
        e = exp( -(x - mu_x) ./ x_0 );
        e(x<mu_x) = 0;
    end

    function e = gaussian(x, sigma_x)
        e = 1 ./ (sqrt(2*pi) .* sigma_x) * exp( -(x .^ 2) ./ (2 * sigma_x.^2) );
    end

    function e = e_conv_g(x, a, x_0, mu_x, sigma_x, B)
        x_ext = pad_x(x, 20);
        evec = trunc_x(x_ext, x_0, mu_x);
        %evec = trunc_x(x_ext, x_0, 0);
        gvec = gaussian(x, sigma_x);
        e_tmp = a ./ x_0 .* conv(evec, gvec, 'same') + B;
        e = interp1(x_ext, e_tmp, x);
    end

    function e = mathematica_conv(x, x_0, mu_x, sigma_x)
        e = 0.5 .* exp( (sigma_x .^ 2 + 2 .* mu_x .* x_0 - 2 .* x_0 .* x) ./ (2 .* x_0 .^ 2) )...
            .* erfc( (sigma_x .^ 2 + x_0 .* (mu_x - x))./(sqrt(2) .* sigma_x .* x_0));
    end
    
    emg_lu = emgfxn_lu(x, ffit.a, ffit.x_0, ffit.mu_x, ffit.sigma_x, ffit.B);
    emg_defoy = emgfxn_defoy(x, ffit.a, ffit.x_0, ffit.mu_x, ffit.sigma_x, ffit.B);
    emg_conv = e_conv_g(x, ffit.a, ffit.x_0, ffit.mu_x, ffit.sigma_x, ffit.B);
    emg_math = mathematica_conv(x, ffit.x_0, ffit.mu_x, ffit.sigma_x);

end

function x2 = pad_x(x,n,direction)
if ~exist('direction','var')
    do_pre = true;
    do_post = true;
elseif strcmpi(direction,'left') || strcmpi(direction,'pre')
    do_pre = true;
    do_post = false;
elseif strcmpi(direction,'right') || strcmpi(direction,'post')
    do_pre = false;
    do_post = true;
else
    do_pre = true;
    do_post = true;
end
del = mean(diff(x));
if do_pre
    pre = fliplr( (min(x)-del):-del:(min(x)-n*del) );
else
    pre = [];
end

if do_post
    post = (max(x) + del):del:(max(x)+n*del);
else
    post = [];
end
x2 = [pre, x, post];
end

