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
%
%   JLL 13 May 2015: added output for scattering weights and averaging
%   kernel, both as the weighted average of clear and cloudy conditions.
%   The averaging kernel uses the preexisting code (just uncommented), the
%   scattering weights I added myself.
%
%   Josh Laughner <joshlaugh5@gmail.com> 

%function [amf, amfCld, amfClr, avgKernel, vcd, vcdAvgKernel] = omiAmfAK2(pTerr, pCld, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile1, no2Profile2, noGhost, ak)
function [amf, amfCld, amfClr, sc_weights, avgKernel, no2Profile3, swPlev, ghost] = omiAmfAK2(pTerr, pCld, cldFrac, cldRadFrac, pressure, dAmfClr, dAmfCld, temperature, no2Profile1, no2Profile2, noGhost, ak)


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


% JLL 18 May 2015..........................................................
% Added preinitialization of these matrices, also nP will be needed to pad
% output vectors from integPr2 to allow concatenation of scattering weights
% vectors into a matrix (integPr2 will return a shorter vector if one or
% both of the pressures to interpolate to is already in the pressure
% vector). We add two to the first dimension of these matrices to make room
% for the two interpolated pressures.
padvec = zeros(1,ndims(no2Profile1));
padvec(1) = 2;
swPlev=zeros(size(no2Profile1)+padvec);
swClr=zeros(size(no2Profile1)+padvec);
swCld=zeros(size(no2Profile1)+padvec);
no2Profile3=zeros(size(no2Profile1)+padvec);
nP = size(swPlev,1);
%..........................................................................

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

    
    % JLL 19 May 2015......................................................
    % Added these lines to interpolate to the terrain & cloud pressures and
    % output a vector - this results in better agreement between our AMF and
    % the AMF calculated from "published" scattering weights.
    [~, ~, this_no2Profile3] = integPr2(no2Profile1(:,i), pressure, pTerr(i), [pTerr(i), pCld(i)]);
    [~,this_swPlev,this_swClr] = integPr2((dAmfClr(:,i).*alpha(:,i)), pressure, pTerr(i), [pTerr(i), pCld(i)]);
    [~,~,this_swCld] = integPr2((dAmfCld(:,i).*alpha(:,i)), pressure, pCld(i), [pTerr(i), pCld(i)]);
    
    if ~iscolumn(this_swPlev)
        E.badvar('this_swPlev','Must be a column vector');
    elseif ~iscolumn(this_swClr)
        E.badvar('this_swClr', 'Must be a column vector');
    elseif ~iscolumn(this_swCld)
        E.badvar('this_swCld', 'Must be a column vector');
    elseif ~iscolumn(this_no2Profile3)
        E.badvar('this_no2Profile3', 'Must be a column vector');
    end
    
    % Pad with NaNs if there are fewer than nP (number of pressures in the
    % input pressure vector + 2 for the interpolated pressures) values.
    % integPr2 outputs vectors with nP values, unless one of the interpolated
    % pressures is already in the input pressure vector.
    this_swPlev = padarray(this_swPlev, nP - length(this_swPlev), nan, 'post');
    this_swClr = padarray(this_swClr, nP -  length(this_swClr), nan, 'post');
    this_swCld = padarray(this_swCld, nP - length(this_swCld), nan, 'post');
    this_no2Profile3 = padarray(this_no2Profile3, nP - length(this_no2Profile3), nan, 'post');
    
    swPlev(:,i) = this_swPlev;
    swClr(:,i) = this_swClr;
    swCld(:,i) = this_swCld;
    no2Profile3(:,i) = this_no2Profile3;
    %......................................................................

end
% Combine clear and cloudy parts of AMFs

amf = cldRadFrac .* amfCld + (1-cldRadFrac).*amfClr;


% JLL 17 May 2015: The ghost correction factor is now pre calculated, so it can
% be returned easier. Also, should the AMF be DIVIDED by it, rather than
% multiplied? It is the ratio of total to visible column, so
% V_total = V_vis * ghost, and V = S/A, so if V_total = S/A_total, and
% A_total = A_vis/ghost, then V_total = S/(A_vis / ghost) = V_vis * ghost.
ghost = vcdGnd ./ (vcdCld .* cldFrac + vcdGnd .* (1.-cldFrac));

if numel(noGhost) == 1;
    if noGhost > 0;
        amf  = amf ./ ghost;
    end
end

amf(~isnan(amf)) = max(amf(~isnan(amf)),1.e-6);   % clamp at min value (2008-06-20), but don't replace NaNs with the min value (2016-05-12)

if numel(ak) == 1;
   if ak > 0;
       % Preallocation added 13 May 2015 - JLL.............................
       avgKernel = nan(size(swPlev));
       sc_weights = nan(size(swPlev));
       % ..................................................................
       % Now compute averaging kernel.............................................

       % These 2 sets of lines are an approximation of what we do in the OMI NO2 algorithm
       for i=1:numel(pTerr)
           %...............................................................
           % JLL 19 May 2015 - pull out the i'th vector, this will allow us
           % to remove nans for AMF calculations where needed, and also
           % check that all vectors have NaNs in the same place.
           swPlev_i = swPlev(:,i);
           swClr_i = swClr(:,i);
           swCld_i = swCld(:,i);
           not_nans_i = ~isnan(swPlev_i) & ~isnan(swClr_i) & ~isnan(swCld_i);
           if ~all(not_nans_i == (~isnan(swPlev_i) | ~isnan(swClr_i) | ~isnan(swCld_i))) && ~all(isnan(swPlev_i)) && ~all(isnan(swClr_i)) && ~all(isnan(swCld_i))
               % Error called if there are NaNs present in some but not all
               % of these vectors AND none of the vectors is all NaNs
               % If one of the vectors is all NaNs, then the mismatch is
               % okay because the AMF will just end up being a NaN anyway - 
               E.callError('nan_mismatch','NaNs are not the same in the swPlev, swClr, and swCld vectors');
           end
           %...............................................................
           
           ii = find(swPlev_i > pTerr(i));
           if min(ii) >= 1;
               swClr_i(ii)=1E-30;
           end
           ii = find(swPlev_i > pCld(i));
           if min(ii) >= 1;
               swCld_i(ii)=1E-30;
           end

           % more temporary testing code - JLL 19 May 2015.................
           ii = find(pressure > pTerr(i));
           if min(ii) >= 1;
               dAmfClr0(ii,i) = 1E-30;
           end
           ii = find(pressure > pCld(i));
           if min(ii) >= 1;
               dAmfCld0(ii,i) = 1E-30;
           end
           %...............................................................
           % Added 14-15 May 2015 to handle outputting scattering weights
           % 
           sc_weights(:,i) = (cldRadFrac(i).*swCld_i + (1-cldRadFrac(i)).*swClr_i);
           %...............................................................
           avgKernel(:,i) = sc_weights(:,i) ./ amf(i); % JLL 19 May 2015 - changed to use the scattering weights we're already calculating.
       end
   end
end


dum=1;

% Integrate NO2 profile with and without averaging kernel .................
%            ii           = find(pressure >= minPressure);
%            vcd(i)          = integPr2(max(no2Profile2(i,ii), 1E-30), pressure(ii), pTerr(i));
%            vcdAvgKernel(i) = integPr2(max(no2Profile2(i,ii).*avgKernel(i,ii), 1E-30), pressure(ii), pTerr(i));
%        end
%   end
%end



