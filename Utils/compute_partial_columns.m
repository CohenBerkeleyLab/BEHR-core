function [ partCol ] = compute_partial_columns( mixingRatio, pressures )
%PARCOL = COMPUTE_PARTIAL_COLUMNS( MIXINGRATIO, PRESSURES )
%   This function takes as input a vector of mixing ratios and the bin center
%   pressure that match them and integrates between each pair of adjacent bins
%   to produce a vector of partial column densities (in molec./cm^2). The
%   mixing ratios should be in parts-per-part and the pressures in hPa.
%
%   Josh Laughner <joshlaugh5@gmail.com> 4 Apr 2016

partCol = nan(numel(mixingRatio)-1,1);
for a=1:numel(mixingRatio)-1
    partCol(a) = integPr2(mixingRatio(a:a+1), pressures(a:a+1), pressures(a));
end

end

