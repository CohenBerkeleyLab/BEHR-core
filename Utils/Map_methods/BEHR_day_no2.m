function [ SumWeightedColumn, SumWeight, Count ] = BEHR_day_no2( OMI, varargin )
%BEHR_DAY_NO2 - handles the weighting of a day's worth of BEHR data

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

p = inputParser;
p.addParameter('mapfield','',@ischar)
p.addParameter('cloud_prod','omi',@(x) ismember(lower(x),{'omi','rad','modis'}));
p.addParameter('cloud_frac_max',0.2,@(x) (isnumeric(x) && isscalar(x)));
p.addParameter('row_anomaly','XTrackFlags', @(x) ismember(x, {'AlwaysByRow', 'RowsByTime', 'XTrackFlags', 'XTrackFlagsLight'}));
p.addParameter('rows',[],@(x) (isnumeric(x) && (numel(x) == 0 || numel(x) == 2)));
p.addParameter('sza', 180, @(x) (isnumeric(x) && isscalar(x) && x >= 0))
p.addParameter('rmserror', Inf, @(x) (isnumeric(x) && isscalar(x) && x >= 0));

p.parse(varargin{:});
pout = p.Results;
mapfield = pout.mapfield;
cloud_prod = pout.cloud_prod;
cloud_frac_max = pout.cloud_frac_max;
row_anomaly = pout.row_anomaly;
rows = pout.rows;
sza = pout.sza;
rmserror = pout.rmserror;

% Check that OMI is some sort of satellite output structure, be it BEHR or
% NASA SP
if ~isstruct(OMI)
    E.badinput('OMI must be a structure');
elseif isfield(OMI,'BEHRColumnAmountNO2Trop')
    usebehr = true;
    if isempty(mapfield)
        mapfield = 'BEHRColumnAmountNO2Trop';
    end
elseif isfield(OMI, 'ColumnAmountNO2Trop');
    usebehr = false;
    if isempty(mapfield)
        mapfield = 'ColumnAmountNO2Trop';
    end
else
    E.badinput('OMI is expected to have a tropospheric column amount NO2 field, BEHR or SP');
end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

SumWeightedColumn = zeros(size(OMI(1).(mapfield)));
SumWeight = zeros(size(OMI(1).(mapfield)));
Count = zeros(size(OMI(1).(mapfield)));

for a=1:numel(OMI)
    omi = OMI(a);
    
    if usebehr
        omi = omi_pixel_reject(omi, cloud_prod, cloud_frac_max, row_anomaly, rows, sza, rmserror); %We will set the area weight to 0 for any elements that should not contribute to the average
    else
        omi = omi_sp_pixel_reject(omi, cloud_prod, cloud_frac_max, row_anomaly, rows);
    end
    
    SumWeightedColumn = nansum2(cat(3,SumWeightedColumn, omi.(mapfield) .* omi.Areaweight),3);
    SumWeight = nansum2(cat(3, SumWeight, omi.Areaweight),3);
    if isfield(omi,'Count')
        this_count = omi.Count;
        this_count(omi.Areaweight <= 0) = 0;
        Count = nansum2(cat(3, Count, this_count),3);
    end
end

end

