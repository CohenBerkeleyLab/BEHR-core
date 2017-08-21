function [ refl, refl_struct ] = coart_sea_reflectance( sza, refl_struct )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



if ~exist('refl_struct', 'var')
    refl_struct = read_coart_html();
end

refl = interp1(refl_struct.sza, refl_struct.alb, sza);

end

function refl = read_coart_html()

mydir = fileparts(mfilename('fullpath'));
coart_file = fullfile(mydir, 'coart.htm');

% Preallocate the SZA and ALB vectors - will clean up unused elements at
% the end
refl.sza = nan(1,20);
refl.alb = nan(1,20);
i_alb = 1;

% For more helpful error messages
line_num = 1;

fid = fopen(coart_file, 'r');
try
    tline = fgetl(fid);
    while ischar(tline)
        % The COART HTML output prints the albedo information in groups of
        % four lines, where each line gives the albedo information for a
        % different measurement altitude and each block is for a different
        % SZA. We will search for lines that define the reflectance at 0 km
        % altitude (so right at the surface).
        vals = strsplit(tline);
        % We need to remove any values that aren't just a number, the
        % numbers could be defined a 1.23E+2 or 4.5e-3, so we look for any
        % characters that aren't a number, a decimal, E, +, -, or NaN. This
        % will usually be extra HTML tags.
        xx = iscellcontents(regexp(vals, '[^\d\.eE\+\-(NaN)]'),'isempty');
        % We also remove any cells that only contain an empty string
        xx = xx & ~iscellcontents(vals, 'isempty');
        vals = vals(xx);
        if numel(vals) == 8
            % The measurement altitude is the third column, the SZA the second.
            % The albedo is the last column, assuming that we take the advice
            % of http://www.oceanopticsbook.info/view/radiative_transfer_theory/level_2/measures_of_reflectance
            % more specifically, the "Light and Water" book they reference in
            % that the "oceanographer's albedo" is the ratio of upwelling to
            % downwelling irradiance.
            %
            % Some other sources (of varying rigor) on the relationship between
            % radiance, irradiance, and albedo/reflectance:
            %   http://www.cesbio.ups-tlse.fr/multitemp/?p=9148
            %   http://ceeserver.cee.cornell.edu/wdp2/cee6100/6100_Labs/Lab03_Fa14_Radiance%20&%20Reflectance.pdf
            meas_alt = str2double(vals{3});
            if meas_alt == 0
                alb = str2double(vals{end});
                if ~isnan(alb)
                    refl.sza(i_alb) = str2double(vals{2});
                    refl.alb(i_alb) = alb;
                    i_alb = i_alb + 1;
                end
            end
        end
        
        line_num = line_num + 1;
        tline = fgetl(fid);
    end
catch err
    fclose(fid);
    rethrow(err);
end

fclose(fid);

refl.sza(i_alb:end) = [];
refl.alb(i_alb:end) = [];

end