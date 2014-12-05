%%omiAmfAK
%%arr 07/24/2008

%..........................................................................
% Given a set of profiles, cloud and terrain parameters, computes AMF 
% and averging kernel profile. Also integrates to get vertical column densities.
% This is a combination of old codes omiAvgK.pro and calcAMF.pro (w/out uncert calcs)
% EJB (2008-05-30)
%
% Inputs:
%  pTerr = pressure of terrain (hPa)
%  pCld = pressure of cloud top (hPa)
%  cldFrac    = geometrical cloud fraction (0.0 to 1.0)
%  cldRadFrac = cloud radiance fraction (0.0 to 1.0)
%  pressure   = pressure vector from highest to lowest pressure (hPa)
%  dAmfClr    = dAMF profile vector for clear skies 
%  dAmfCld    = dAMF profile vector for overcast skies 
%  no2Profile1= NO2 mixing ratio profile vector for AMF calculation (cm-3)
%  no2Profile2= NO2 mixing ratio profile vector for integration (cm-3)
%  temperature= temperature profile vector (K)
%
% Outputs:
%  amf    = air mass factor
%  amfCld = component of amf over cloudy part of scene (if any)
%  amfClr = component of amf over clear  part of scene (if any)
%
% Set keyword ak to compute these additional outputs:
%  avgKernel = averaging kernel vector
%  vcd = directly integrated NO2 profile (cm-2)
%  vcdAvgKernel = integrated product of NO2 profile and avg kernel (cm-2)
%
% Set keyword noghost for amf based on visible column only (otherwise, assume column to ground)
%
%    omiAmfAK, pTerr, pCld,  cldFrac,  cldRadFrac,  noGhost=noGhost,  ak=ak,         $ ;;scalar inputs
%              pressure,  dAMFclr,  dAMFcld,  temperature, no2Profile1, no2Profile2, $ ;;vector inputs
%              amf, amfCld, amfClr, avgKernel, vcd, vcdAvgKernel                       ;;outputs%
%
%..........................................................................

%function [amf, amfCld, amfClr, avgKernel, vcd, vcdAvgKernel] = omiAmfAK2(pTerr, pCld, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile1, no2Profile2, noGhost, ak)
function [amf, amfCld, amfClr] = omiAmfAK2(pTerr, pCld, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile1, no2Profile2, noGhost, ak)


% Each profile is expected to be a column in the no2Profile matrix.  Check
% for this by ensuring that the first dimension of both profile matrices
% has the same length as the pressure vector
E = JLLErrors;
if size(no2Profile1,1) ~= length(pressure) || size(no2Profile2,1) ~= length(pressure);
    error(E.callError('profile_input','Profiles must be column vectors in the input matrices.  Ensure size(no2Profile,1) == length(pressure)'));
end
if size(dAmfClr,1) ~= length(pressure) || size(dAmfCld,1) ~= length(pressure);
    error(E.callError('dAmf_input','dAMFs must be column vectors in the input matrices.  Ensure size(dAmfxxx,1) == length(pressure)'));
end
if size(temperature,1) ~= length(pressure)
    error(E.callError('temperature_input','temperature must be a column vector.  Ensure size(temperature,1) == length(pressure)'));
end

minPressure = 200;  % defined tropopause pressure (hPa)

alpha = 1 - 0.003 * (temperature - 220);   % temperature correction factor vector
alpha_i=max(alpha,0.1);
alpha = min(alpha_i,10);

%vcdCld = 0;     amfClr = 0;     amfCld = 0;
dAmfClr0 = dAmfClr;     dAmfCld0 = dAmfCld;


% Integrate to get clear and cloudy AMFs....................................
vcdGnd=zeros(size(pTerr));
vcdCld=zeros(size(pTerr));
amfClr=zeros(size(pTerr));
amfCld=zeros(size(pTerr));

for i=1:numel(pTerr)
    vcdGnd(i) = integPr2(no2Profile1(:,i), pressure, pTerr(i));             
    if cldFrac(i) ~= 0 && cldRadFrac(i) ~= 0;
        vcdCld(i) = integPr2(no2Profile1(:,i), pressure, pCld(i));
    else
        vcdCld(i)=0;
    end
    if cldFrac(i) ~= 1 && cldRadFrac(i) ~= 1;
        amfClr(i) = integPr2((no2Profile1(:,i).*dAmfClr(:,i).*alpha(:,i)), pressure, pTerr(i)) ./ vcdGnd(i);
    else
        amfClr(i)=0;
    end
    if cldFrac(i) ~= 0 && cldRadFrac(i) ~= 0;
        amfCld(i) = integPr2((no2Profile1(:,i).*dAmfCld(:,i).*alpha(:,i)), pressure, pCld(i)) ./ vcdGnd(i);
    else
        amfCld(i)=0;
    end
end
% Combine clear and cloudy parts of AMFs

amf = cldRadFrac .* amfCld + (1-cldRadFrac).*amfClr;

if numel(noGhost) == 1;
    if noGhost > 0;
        amf  = amf .* vcdGnd ./ (vcdCld.*cldFrac  +  vcdGnd.*(1.-cldFrac));
    end
end

amf = max(amf,1.e-6);   % clamp at min value (2008-06-20)

%if numel(ak) == 1;
%    if ak > 0;

% Now compute averaging kernel.............................................

% These 2 sets of lines are an approximation of what we do in the OMI NO2 algorithm
%        for i=1:numel(pTerr)
%            ii = find(pressure > pTerr(i));
%            if min(ii) >= 1;
%                dAmfClr0(i,ii)=1E-30;
%            end
%            ii = find(pressure > pCld(i));
%            if min(ii) >= 1;
%                dAmfCld0(i,ii)=1E-30;
%            end

%            avgKernel(:,i) = (cldRadFrac(i).*dAmfCld0(:,i) + (1-cldRadFrac(i)).*dAmfClr0(:,i)) .* alpha(:,i) ./ amf(i);

% Integrate NO2 profile with and without averaging kernel .................
%            ii           = find(pressure >= minPressure);
%            vcd(i)          = integPr2(max(no2Profile2(i,ii), 1E-30), pressure(ii), pTerr(i));
%            vcdAvgKernel(i) = integPr2(max(no2Profile2(i,ii).*avgKernel(i,ii), 1E-30), pressure(ii), pTerr(i));
%        end
%   end
%end



