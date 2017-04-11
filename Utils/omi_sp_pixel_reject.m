function [ omi ] = omi_sp_pixel_reject( omi_in, cloud_type, cloud_frac, rowanomaly_mode, rows )
%omi_pixel_reject: Set areaweight to 0 for any pixels that will adversely
%affect the accuracy of OMI SP maps.
%   There are a number of criteria that need to be evaluated for an OMI
%   pixel before it can be reliably used as an NO2 measurement.  This
%   function will set the areaweight value to 0 for any pixel which fails
%   these criteria.  This function is the SP counterpart to
%   "omi_pixel_reject," and has BEHR criteria replaced with regular OMI
%   criteria:
%       - Checking for column amount > 1E17 and column amount = NaN, the
%       field has been changed from "BEHRColumnAmountNO2Trop" to the OMI
%       equivalet "ColumnAmountNO2Trop" 
%        - Cloud fractions will only check against OMI cloud fractions,
%        since obviously MODIS cloud fractions are unavailable.
%
%   Inputs:
%       omi_in: An OMI structure, the result of running BEHR_main.m.
%
%       cloud_type = 'geo' or 'rad' (geometric or radiative), 'omi' is also
%           permitted and will use geometric cloud fractions.
%
%       cloud_frac = maximum allowed gloud fraction (0-1.0)
%
%       rowanomaly_mode: Can be one of 4 values, 'AlwaysByRow',
%          'RowsByTime', 'XTrackFlags', and 'XTrackFlagsLight'.
%          'AlwaysByRow' will always reject pixels in the row affected by
%          the anomaly, regardless of whether the data occurs before or
%          after the anomaly reached that row. 'XTrackFlags' uses the row
%          anomaly flag in the OMNO2 product to reject pixels.  These are
%          the two recommended methods.
%
%       rows: An optional variable that will be used to remove certain rows
%           from the calculation, for viewing angle dependence tests. If
%           empty, no rows will be removed. (This allows the user to use
%           this function even if a Rows field is not present in the OMI
%           structure). It should be either empty or a two element vector
%           giving a min and max (inclusive) value for the row numbers,
%           which are 0 to 59.
%
%   The rejection criteria are:
%       Column amount < 0: Likely a fill value or unphysical result.
%       VCD Quality Flag is not a even number: the VCD quality flag is a
%          bit array, with the least significant bit as a summary.  If this
%          bit is 1, then some error occured during the NASA retrieval and
%          we should ignore this pixel.
%       Cloud fraction too great: Cloudy pixels will have much of the
%          tropospheric column obscured, and should not contribute to the
%          average.
%       Column amount > 1E17: This magnitude of tropospheric column is
%          known to be affected by the row anomaly.
%       Column amount is NaN: NaN indicates some failure, either averaging
%          0 values, or another mathematical mistake.
%       Row anomaly: see http://www.knmi.nl/omi/research/product/rowanomaly-background.php

if ~exist('rows','var')
    rows = [];
end
omi = omi_in;

%omi.Areaweight(omi.ColumnAmountNO2Trop<=0) = 0; %Do not average in negative tropospheric column densities

if iscell(omi.vcdQualityFlags) % The flags may be a cell array or not, depending on whether this is for Data or OMI (gridded) structure
    for a=1:numel(omi.vcdQualityFlags);
        if any(mod([omi.vcdQualityFlags{a}],2)~=0)
            omi.Areaweight(a) = 0; %If any of the vcdQualityFlags value is not even (least significant bit ~= 0), do not include this element
        end
    end
else
    omi.Areaweight(mod(omi.vcdQualityFlags,2)~=0) = 0;
end

if strcmpi(cloud_type,'rad')
    omi.Areaweight(omi.CloudRadianceFraction > cloud_frac) = 0;
    omi.Areaweight(isnan(omi.CloudRadianceFraction)) = 0; % Reject pixels with NaN for cloud fraction
    omi.Areaweight(omi.CloudRadianceFraction < 0) = 0; % Reject pixels with fill values for cloud fraction
else 
    omi.Areaweight(omi.CloudFraction > cloud_frac) = 0;
    omi.Areaweight(isnan(omi.CloudFraction)) = 0; % Reject pixels with NaN for cloud fraction
    omi.Areaweight(omi.CloudFraction < 0) = 0; % Reject pixels with fill values for cloud fraction
end


%omi.Areaweight(omi.ColumnAmountNO2Trop > 1E17) = 0; %Do not include the element if the NO2 column is too great.  These are known to be affected by the row anomaly (Bucsela 2013, Atmos. Meas. Tech. 2607)
hh=find(isnan(omi.ColumnAmountNO2Trop)); omi.ColumnAmountNO2Trop(hh)=0; omi.Areaweight(hh)=0; %Set any column NaNs to 0 and do not include them in the average

xx = omi_rowanomaly(omi,rowanomaly_mode); %Remove elements affected by the row anomaly.
omi.Areaweight(xx) = 0;

% Remove elements by row, if the user specifies a row range.
if ~isempty(rows)
    rr = omi.Row < min(rows) | omi.Row > max(rows);
    omi.Areaweight(rr) = 0;
end

end

