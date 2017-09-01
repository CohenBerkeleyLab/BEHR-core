function [ omi ] = omi_pixel_reject( omi, reject_mode, varargin )
%omi_pixel_reject: Set areaweight to 0 for any pixels that will adversely
%affect the accuracy of the BEHR NO2 map.
%   There are a number of criteria that need to be evaluated for an OMI
%   pixel before it can be reliably used as an NO2 measurement.  This
%   function will set the areaweight value to 0 for any pixel which fails
%   these criteria.
%
%   Inputs:
%       omi_in: An OMI structure, the result of running BEHR_main.m.
%       cloud_type: 'modis' or 'omi'. Sets which
%          cloud product to use when rejecting pixels by cloud fraction.
%       cloud_frac: The maximum allowed cloud fraction in a pixel.  If
%          using OMI cloud product, 0.2 is recommended. With MODIS, 0 is
%          recommended.
%       rowanomaly_mode: Can be one of 4 values, 'AlwaysByRow',
%          'RowsByTime', 'XTrackFlags', and 'XTrackFlagsLight'.
%          'AlwaysByRow' will always reject pixels in the row affected by
%          the anomaly, regardless of whether the data occurs before or
%          after the anomaly reached that row. 'XTrackFlags' uses the row
%          anomaly flag in the OMNO2 product to reject pixels.  These are
%          the two recommended methods.
%       rows: An optional variable that will be used to remove certain rows
%           from the calculation, for viewing angle dependence tests. If
%           empty, no rows will be removed. (This allows the user to use
%           this function even if a Rows field is not present in the OMI
%           structure). It should be either empty or a two element vector
%           giving a min and max (inclusive) value for the row numbers,
%           which are 0 to 59.
%
%   The rejection criteria are:
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

E = JLLErrors;

p = inputParser;
p.addOptional('reject_details', struct(), @isstruct);
p.addParameter('weight_field', 'Areaweight', @ischar);

p.parse(varargin{:});
pout = p.Results;

reject_details = pout.reject_details;
weight_field = pout.weight_field;

allowed_reject_modes = {'none', 'behr', 'detailed'};
if ~ischar(reject_mode) || ~ismember(reject_mode, allowed_reject_modes)
    E.badinput('REJECT_MODE must be one of the strings %s', strjoin(allowed_reject_modes, ', '))
end

cld_fields = struct('omi', 'CloudFraction', 'modis', 'MODISCloud', 'rad', 'CloudRadianceFraction');
allowed_row_modes = {'AlwaysByRow', 'RowsByTime', 'XTrackFlags', 'XTrackFlagsLight'};
if strcmp(reject_mode, 'detailed')
    req_details_fields = {'cloud_type', 'cloud_frac', 'row_anom_mode'};
    if any(~ismember(req_details_fields, fieldnames(reject_details)))
        E.badinput('The REJECT_DETAILS struct must have the following fields for REJECT_MODE == "%s": %s', reject_mode, strjoin(req_details_fields, ', '));
    end
    
    req_details = {'cloud_type', 'cloud_frac', 'row_anom_mode', 'check_behr_amf'};
    if ~exist('reject_details', 'var') || ~isstruct(reject_details) || any(~isfield(reject_details, req_details))
        E.badinput('If REJECT_MODE == ''detailed'', the third input must be a structure with fields %s', strjoin(req_details, ', '));
    elseif ~ismember(reject_details.cloud_type, fieldnames(cld_fields))
        E.badinput('The "cloud_type" field of the REJECT_DETAIL struct must be one of the strings %s', strjoin(fieldnames(cld_fields)));
    elseif ~isnumeric(reject_details.cloud_frac) || ~isscalar(reject_details.cloud_frac)
        E.badinput('The "cloud_frac" field of the REJECT_DETAIL struct must be a scalar number');
    elseif ~ismember(reject_details.row_anom_mode, allowed_row_modes)
        E.badinput('The "row_anom_mode" field of the REJECT_DETAIL struct must be one of the strings: %s', strjoin(allowed_row_modes, ', '));
    end
end



switch reject_mode
    case 'none'
        req_fields = {};
    case 'behr'
        req_fields = {'Areaweight', 'BEHRQualityFlags'};
    case 'detailed'
        req_fields = {'Areaweight', 'XTrackQualityFlags', 'VcdQualityFlags', cld_fields.(reject_details.cloud_type)};
    otherwise
        E.notimplemented('Required fields not defined for REJECT_MODE == ''%s''', reject_mode);
end

if ~isstruct(omi) || ~isscalar(omi) || (numel(req_fields) > 0 && any(~isfield(omi, req_fields)))
    E.badinput('OMI must be a scalar structure with the fields %s', strjoin(req_fields, ', '))
end

%%%%%%%%%%%%%%%%%
% MAIN FUNCTION %
%%%%%%%%%%%%%%%%%

if strcmpi(reject_mode, 'none')
    % do nothing. this allow us to average fields that should not have any
    % bad values, or at least their bad values are unrelated to what
    % generates bad values in the VCD.
elseif strcmpi(reject_mode, 'behr')
    omi.(weight_field) = reject_by_behr_flags(omi.BEHRQualityFlags, omi.(weight_field));
elseif strcmpi(reject_mode, 'detailed')
    sel_cloud_field = cld_fields.(reject_details.cloud_type);
    omi.(weight_field) = reject_by_details(reject_details, omi.(sel_cloud_field), omi.VcdQualityFlags, omi.XTrackQualityFlags, omi.BEHRAMFTrop, omi.Row, omi.(weight_field));
end

end

function areaweight = reject_by_behr_flags(behr_flags, areaweight)
areaweight(behr_flags>0) = 0;
end

function areaweight = reject_by_details(details, cloud_frac, vcd_flags, xtrack_flags, behr_amf, rows, areaweight)

if ~isfield(details, 'rows')
    details.rows = [];
end
if ~isfield(details, 'szalim')
    % not currently implemented since I didn't grid SZA 
    details.szalim = 180;
end

areaweight(mod(vcd_flags,2) ~= 0) = 0;
areaweight(cloud_frac > details.cloud_frac) = 0;

% As of 30 Aug 2017, the PSM OMI structure does not include date. This will
% need to be rectified in the future. For now, just pass 0 as the date.
areaweight(omi_rowanomaly(xtrack_flags, rows, 0, details.row_anom_mode)) = 0;

if details.check_behr_amf
    areaweight( isnan(behr_amf) | behr_amf <= behr_min_amf_val ) = 0;
end

if ~isempty(details.rows)
    areaweight( rows < min(details.rows) | rows > max(details.rows) ) = 0;
end

%areaweight( sza > details.szalim ) = 0;

end

