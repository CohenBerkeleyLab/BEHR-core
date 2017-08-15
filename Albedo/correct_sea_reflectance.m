function [ data ] = correct_sea_reflectance( data, this_month, gome_ratio, gome_lon, gome_lat )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

for a=1:numel(data.Longitude)
    if data.AlbedoOceanFlag(a)
        % Must have lat and lon in this order to make a valid mesh grid for
        % interp2 to work.
        correction = interp2(gome_lat, gome_lon, gome_ratio(:,:,this_month), data.Latitude(a), data.Longitude(a));
        data.MODISAlbedo(a) = data.MODISAlbedo(a) * correction;
    end
end

end

