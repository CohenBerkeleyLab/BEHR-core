%%integPr
%%arr 07/28/2008

%..........................................................................
% Integrates vector of mixing ratios (cm-3) above pressureSurface as function of pressure (hPa)
% to get vertical column densities (cm-2). Computes piecewise in layers between two pressures.
% Assumes exponential variation between mixing ratios that are positive or zero (zero is treated as 1.e-30).
% If one or both mixing ratios is negative, assumes constant (average=(f1+f2)/2) mixing ratio in the layer.
% The uncertainties are computed assuming constant mixing ratio in each layer (using exponential
% variation can lead to strange (large) results when f1*p1 is approx f2*p2. This sometimes happens
% when using averaging kernels).
%
% Arguments (vector indices are i = 0,1,2,...n-1):
%
%  mixingRatio(i)   = input volume mixing ratios (no units)
%  pressure(i)      = vector of input pressures from largest to smallest (hPa)
%  pressureSurface  = surface pressure (where integration starts) (hPa)
%  interpPres       = pressures to interpolate the output profile to (hPa) ( added by JLL for BEHR scattering weights )
%
% Restrictions:
%  Pressures must be greater than zero and monotonically decreasing
%
% Equation for exponential variation:
%  For a segment of the vertical column between ambient pressures p1 and p2,
%  at which the trace gas volume mixing ratios are f1 and f2, the vertical
%  column density VCD (assuming number density varies exponentially with p) is
%
%    VCD   =  (f1 * p1  -  f2 * p2) /  [(b + 1) * m * g]
%
%    where:   b = ALOG( f2/f1 ) / ALOG( p2/p1 )
%
%   Also, for interpolation of f between p1 and p2:
%
%   f  =  f1 * (p/p1)^b
%
% Note:  Can get large uncertainties when  f1*p1 = f2*p2
%
%..........................................................................

%function vcd = integPr(mixingRatio, pressure, pressureSurface, mixingRatioStd, vcdStd, corrLength) %mixingRatioStd, vcdStd, corrLength = 0 or 1
function [vcd, p_out, f_out] = integPr2(mixingRatio, pressure, pressureSurface, interpPres)

E = JLLErrors;

if nargin < 4;
    interpPres = [];
    if nargout > 1
        E.callError('nargout','Without any interpPres values, p_out and f_out will not be set');
    end
else
    % Make sure the interpolation pressure is not less than the 
    interpPres = max(interpPres,min(pressure));
end

if any(pressure<0)
    E.badinput('PRESSURE must be all >= 0')
elseif any(diff(pressure)>0);
    E.badinput('PRESSURE must be monotonically decreasing')
end

if ~isscalar(pressureSurface)
    E.badinput('PRESSURESURFACE must be a scalar')
end

%   mean molecular mass (kg)  *  g (m/s2)  *  (Pa/hPa)   *   (m2/cm2)
mg  = (28.97/6.02E23)*1E-3    *   9.8      *    1E-2     *     1E4;

fmin = 1E-30;

vcd  = 0;
f    = mixingRatio;
p    = pressure;
p0   = max(pressureSurface, min(pressure));
n    = numel(p);

dvcd     = 0;
deltaVcd = zeros(numel(f),1); % Changed to make these vectors on 9/26/2014 JLL
df       = zeros(numel(f),1);



numIP = numel(interpPres);
if numIP > 0
    if ~iscolumn(p); p_out = p'; 
    else p_out = p;
    end
    f_out = f;
    for a=1:numIP
        if all(p~=interpPres(a))
            f_i = interpolate_surface_pressure(p,f,interpPres(a));
            bottom = p_out > interpPres(a);
            top = p_out < interpPres(a);
            p_out = [p_out(bottom); interpPres(a); p_out(top)];
            f_out = [f_out(bottom); f_i; f_out(top)];
        end
    end
end

f0 = interpolate_surface_pressure(p,f,p0);

p(i0) = p0;
f(i0) = f0;


% Integrate................................................................
if isnan(pressureSurface)
    vcd = nan;
    return
end

for i = 1:n-1;
    deltaVcd(i) = (f(i) + f(i + 1)) .* (p(i) - p(i + 1)) ./ (2 * mg);  %assume const mixing ratio in each layer
    b = (log(max(f(i + 1),fmin)) - log(max(f(i),fmin))) ./ (log(p(i + 1)) - log(p(i)));
    
    if f(i) >= 0 && f(i + 1)>=0 && abs(b + 1) >= 0.01;
        deltaVcd(i) = (f(i)*p(i) - f(i + 1)*p(i + 1)) ./ (b + 1) ./ mg;    %assume exponential variation in each layer
    end
    
end

if any(isnan(deltaVcd(i0:n-1))) && ~all(isnan(deltaVcd(i0:n-1)))
    warning('NaNs detected in partial columns. They will not be added into the total column density.')
end
vcd = nansum2(deltaVcd(i0:n-1));


    function [f0] = interpolate_surface_pressure(p_in,f_in,p0_in)
        pp=find(p_in>=p0_in); if isempty(pp); pp=0; end
        i0 = min(max((max(pp)),1),(n - 1));
        f0 = interp1(log([p_in(i0),p_in(i0 + 1)]), [f_in(i0),f_in(i0 + 1)], log(p0_in),'linear','extrap');   %assume linear variation in surface layer
        
        b = (log(max(f_in(i0 + 1),fmin)) - log(max(f_in(i0),fmin))) ./ (log(p_in(i0 + 1)) - log(p_in(i0)));
        if f_in(i0) >= 0 && f_in(i0 + 1) >= 0 && abs(b + 1) >= 0.01;
            f0 = f_in(i0) * ((p0_in./p_in(i0))^b);                                    %assume exponential variation in surface layer
        end
        
    end
end
