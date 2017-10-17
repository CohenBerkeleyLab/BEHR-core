function [ refl, refl_struct ] = mobley_sea_refl( sza, vza, raa, refl_struct )
%MODIS_SEA_REFL Look up sea reflectivity from Mobley 2015
%   The MCD43C1 product does not give BRDF coefficients over ocean.
%   Vasilkov et al. (2017) handled this by using two models that combined
%   direct (specular) and volumetric reflection from the ocean. However, as
%   far as I could tell, that seemed like it would require calculating the
%   normal for a surface that would yield specular reflection from the
%   solar to viewing angle, which is a difficult calculation. To do the
%   volumetric/diffuse calculation would require measurements of
%   chorolphyll concentration as well, and I was keen to avoid adding
%   another dataset dependency to BEHR.
%
%   Fortunately, I found a paper by Curtis D. Mobley, which calculated
%   ocean reflectances using a Monte Carlo method. He provided a table of
%   reflectances, which was downloadable at 
%
%   http://www.oceanopticsbook.info/view/remote_sensing/level_3/surface_reflectance_factors
%
%   These are calculated for 550 nm, and from the figure at the above link,
%   it looks like this means we can expect at least a 20% error. The usual
%   way to correct that (e.g. McLinden et al. 2014, 
%
%   REFL = MOBLEY_SEA_REFL( SZA, VZA, RAA ) returns the reflectivity for
%   given values of SZA, VZA, and RAA. If the inputs are arrays, REFL will
%   be an array of the same size. Note that SZA, VZA, and RAA are expected
%   to be defined in degrees, and RAA is expected to be defined as it is in
%   the TOMRAD table, I believe. From p. 4844 of Mobley 2015:
%       "...the azimuthal viewing direction, \phi_v, which is measured
%       relative to the Sun's azimuthal direction. That is, \phi_v = 0
%       corresponds to looking towards the Sun, \phi_v = 90 deg is looking
%       at right angles to the Sun's incident rays, and \phi_v = 180 deg is
%       looking away from the Sun."
%
%   [REFL, REFL_STRUCT = MOBLEY_SEA_REFL( ___ ) returns the structure
%   created by reading the Mobley text file so that it can be passed in
%   with the next syntax:
%
%   REFL = MOBLEY_SEA_REFL( ___, REFL_STRUCT ) passes the structure back
%   into this function so that you do not need to re-read the file, which
%   can be the bulk of the execution time.
%
%   Sources:
%
%       Curtis D. Mobley, "Polarized reflectance and transmittance
%       properties of windblown sea surfaces," Appl. Opt. 54, 4828-4849
%       (2015), doi: 10.1364/AO.54.004828
%
%       Vasilkov et al., "Accounting for the effects of surface BRDF on
%       satellite cloud and trace gas retrievals: a new approach based on
%       geometry-dependent Lambertian equivalent reflectivity applied to
%       OMI algorithms," Atmos. Meas. Tech., 10, 333-349 (2017), doi:
%       10.5194/amt-10-333-2017

if ~isequal(size(sza), size(vza)) || ~isequal(size(sza), size(raa)) || ...
        ~isnumeric(sza) || ~isnumeric(vza) || ~isnumeric(raa)
    E.badinput('SZA, VZA, and RAA must all be numeric arrays of the same size')
end

if ~exist('refl_struct', 'var')
    refl_struct = read_in_table();
end

assumed_wind = repmat(5, size(sza));

refl = interpn(refl_struct.wind_speed, refl_struct.sza, refl_struct.vza, refl_struct.raa,...
               refl_struct.refl, assumed_wind, sza, vza, raa);

end

function refl = read_in_table()
E = JLLErrors;
E.addCustomError('table_read', 'Problem reading Mobley 2015 table: %s');

table_file = fullfile(behr_repo_dir, 'Albedo', 'Mobley2015_SeaSurfRefl.txt');

fid = fopen(table_file, 'r');
table_started = false;

tline = fgetl(fid);

wind_dim = [0, 2, 4, 5, 6, 8, 10, 12, 14, 15];
sza_dim = [0:10:80, 87.5];
vza_dim = [0:10:80, 87.5];
raa_dim = 0:15:180;

curr_wind = nan;
curr_sza = nan;

refl_table = nan(length(wind_dim), length(sza_dim), length(vza_dim), length(raa_dim));

while ischar(tline)
    if ~table_started && ~isempty(strfind(tline, 'WIND SPEED'))
        table_started = true;
    end
    
    if table_started
        if ~isempty(strfind(tline, 'WIND SPEED'))
            tmpline = strsplit(tline, ';');
            tmpline{1} = strsplit(tmpline{1}, '=');
            tmpline{2} = strsplit(tmpline{2}, '=');
            curr_wind = str2double(tmpline{1}{2});
            curr_sza = str2double(tmpline{2}{2});
        else
            if isnan(curr_wind) || isnan(curr_sza)
                E.callCustomError('table_read', 'Tried to read reflectance value before wind speed and SZA read');
            end
            
            tmpline = strsplit(strtrim(tline));
            curr_vza = str2double(tmpline{1});
            curr_raa = str2double(tmpline{2});
            curr_refl = str2double(tmpline{3});
            
            inds = find_table_inds(curr_wind, curr_sza, curr_vza, curr_raa);
            refl_table(inds(1), inds(2), inds(3), inds(4)) = curr_refl;
        end
    end
    
    tline = fgetl(fid);
end

fclose(fid);

% There will be some NaNs remaining in the table because the Mobley file
% does not include output for multiple RAAs when VZA = 0 (since RAA does
% not actually mean anything in that case). We want to replicate the RAA =
% 0 values to fill those nans

% Verified that this copies the values for RAA = 0 VZA = 0 to all values of
% RAA for VZA = 0 properly with (let refl_table_old be refl_table before
% the subsitution)
%   notnans = ~isnan(refl_table_old);
%   isequal(refl_table(notnans), refl_table_old(notnans)) % true
%   sum(diff(refl_table(:,:,1,:),[],4),4) % all values are 0
vza_eq_0 = refl_table(:,:,1,1);
refl_table(:,:,1,2:end) = repmat(vza_eq_0, 1, 1, length(raa_dim)-1); 

refl = struct('refl', refl_table, 'wind_speed', wind_dim, 'sza', sza_dim, 'vza', vza_dim, 'raa', raa_dim,...
    'refl_dimensions', {{'wind_speed', 'sza', 'vza', 'raa'}});


    function found_inds = find_table_inds(wind, sza, vza, raa)
        found_inds = nan(1,4);
        vals = [wind, sza, vza, raa];
        dims = {wind_dim, sza_dim, vza_dim, raa_dim};
        for a=1:numel(vals)
            xx = vals(a) == dims{a};
            if sum(xx) ~= 1
                dim_names = {'wind', 'SZA', 'VZA', 'RAA'};
                E.callCustomError('table_read', sprintf('Could not identify indices for %s = %.1f', dim_names{a}, vals(a)));
            else
                found_inds(a) = find(xx);
            end
        end
    end

end