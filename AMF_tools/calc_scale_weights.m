function [ Weights ] = calc_scale_weights(  )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;

% Load the June (2013) 100% and 50% emission sets of a priori profiles
remotedir = '/Volumes/share2/USERS/LaughnerJ/WRF/';
fullemis_file = fullfile(remotedir,'US_BEHR_FULLEMIS','monthly','WRF_BEHR_monthly_2013-06-30.nc');
halfemis_file = fullfile(remotedir,'US_BEHR_HALFEMIS','monthly','WRF_BEHR_monthly_2013-06-30.nc');

lon = ncread(fullemis_file, 'XLONG');
lat = ncread(fullemis_file, 'XLAT');
lonhalf = ncread(halfemis_file, 'XLONG');
lathalf = ncread(halfemis_file, 'XLAT');
if any(abs(lon(:) - lonhalf(:)) > 0.01) || any(abs(lat(:) - lathalf(:)) > 0.01)
    E.badgeo('The full and half emissions files have inconsistent lat/lon grids');
end

no2_full = ncread(fullemis_file, 'no2');
no2_half = ncread(halfemis_file, 'no2');

%Weights.June.w = compute_weights(no2_full, no2_half, 1, 0.5);
Weights.June.w = compute_weights(no2_half, no2_full, 0.5, 1);
Weights.lon = double(lon);
Weights.lat = double(lat);

end

function weights = compute_weights( c_year, c_base, e_year, e_base )
% Actually computes the weights for how the profiles are affected by the
% change in emissions. This form is such that parts of the profile with
% little change will have a weight of nearly 0, while those that change by
% the same proportion as the emissions have a weight of 1.
weights = ((c_year - c_base) ./ c_base) ./ ((e_year - e_base) ./ e_base);
weights = double(weights);
end