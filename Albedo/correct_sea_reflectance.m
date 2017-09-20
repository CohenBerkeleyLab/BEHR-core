function [ data ] = correct_sea_reflectance( data, this_month, gome_ratio, gome_lon, gome_lat )
%CORRECT_SEA_REFLECTANCE Correct sea reflectance based on the ratio observed by GOME at two wavelengths
%   [ DATA ] = CORRECT_SEA_REFLECTANCE( DATA, THIS_MONTH, GOME_RATIO,
%   GOME_LON, GOME_LAT ) Takes in a BEHR Data structure and the numerical
%   month, as well as GOME data read by read_gome2_ler. This multiplies the
%   MODISAlbedo field in DATA by the GOME_RATIO value interpolated to each
%   pixel's lat/lon, if the AlbedoOceanFlag is set to true for that pixel.

for a=1:numel(data.Longitude)
    if data.AlbedoOceanFlag(a)
        % Must have lat and lon in this order to make a valid mesh grid for
        % interp2 to work.
        correction = interp2(gome_lat, gome_lon, gome_ratio(:,:,this_month), data.Latitude(a), data.Longitude(a));
        data.MODISAlbedo(a) = data.MODISAlbedo(a) * correction;
    end
end

end

