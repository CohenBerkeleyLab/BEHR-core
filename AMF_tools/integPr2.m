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
%
% Keywords:
%
%  mixingRatioStd(i)= input vector of uncertainty (standard dev) in mixing ratio (cm-2)
%  vcdStd           = output uncertainty in integrated column (cm-2)
%  corrLength       = correlation length defined in log pressure space (no units).
%                     If variation in mixing ratios between pressures p1 and p2
%                     are correlated, then corrLength = ALOG(p1) - ALOG(p2).
%                     A typical value is corrLength = 0.1 to 0.2.  If no value is given,
%                     assume default corrLength near zero. Note: when mixing ratios
%                     are obtained from binning normally distributed data, and uncertainties
%                     are computed as the standard error in the bins, then the default
%                     value of corrLength = 0 yields an integrated column uncertainty
%                     that is independent of the bin width.
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

if nargin < 4;
    interpPres = [];
    if nargout > 1
        error('integPr2:nargout','Without any interpPres values, p_out and f_out will not be set');
    end
else
    interpPres = max(interpPres,min(pressure));
end


corrLen = 0.1;
%{
if numel(corrLength) == 1;
    corrLen = corrLength;
end
%}
presSurface = max(pressure);
if numel(pressureSurface) == 1;
    presSurface = pressureSurface;
end

%   mean molecular mass (kg)  *  g (m/s2)  *  (Pa/hPa)   *   (m2/cm2)
mg  = (28.97/6.02E23)*1E-3    *   9.8      *    1E-2     *     1E4;

fmin = 1E-30;

vcd  = 0;
f    = mixingRatio;
p    = pressure;
p0   = max(presSurface, min(pressure));
n    = numel(p);

dvcd     = 0;
deltaVcd = zeros(numel(f),1); % Changed to make these vectors on 9/26/2014 JLL
df       = zeros(numel(f),1);
%{
if numel(mixingRatioStd) == n;
    df = max(mixingRatioStd,0);
end
%}
if any(diff(pressure)>0);
    error('integPr2:pressure','Pressure must be monotonically decreasing')
end


% Set starting point at surface pressure...................................
%
% pp=find(p>=p0); if isempty(pp); pp=0; end
% i0 = min(max((max(pp)),1),(n - 1));
% f0 = interp1(log([p(i0),p(i0 + 1)]), [f(i0),f(i0 + 1)], log(p0),'linear','extrap');   %assume linear variation in surface layer
%
% b = (log(max(f(i0 + 1),fmin)) - log(max(f(i0),fmin))) ./ (log(p(i0 + 1)) - log(p(i0)));
% if f(i0) >= 0 && f(i0 + 1) >= 0 && abs(b + 1) >= 0.01;
%     f0 = f(i0) * ((p0./p(i0))^b);                                    %assume exponential variation in surface layer
% end

numIP = numel(interpPres);
if numIP > 0
    p_out = p';
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

%{
if numel(mixingRatioStd) == n;   % Do uncertainty calculation..............

 % Compute sig(i) = uncertainty in vcd due to mixing ratio i

    sig   = zeros(numel(f));
    for i=1:n;
        sig(i) = df(i) .* (p(max((i-1),1)) - p(min((i+1),(n-1)))) ./ (2 .* mg);
    end

 % Combine the sig(i) using correlation length
 
    vcdVar = 0;
    for i = i0:n; %%%%not sure about these indicies
        for j = i0:n;
        vcdVar = vcdVar + sig(i) .* sig(j) .* exp(-((log(p(i)) - log(p(j))) ./ corrLen)^2);
        end
    end

    vcdStd = sqrt(max(vcdVar,0));
end
%}


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
