function [ gome_ler_ratio, gome_lon, gome_lat ] = read_gome2_ler( ler_file, lonlim, latlim )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% The GOME-2A file was obtained from http://www.temis.nl/surface/gome2_ler.html


buffer = 25; % buffer in degrees around the lon/lat limits to ensure we get enough data to interpolate to edge pixels
gome_ler = double(h5read(ler_file, '/Minimum_LER'));
gome_lon = double(h5read(ler_file, '/Longitude'));
gome_lat = double(h5read(ler_file, '/Latitude'));
gome_wavelength = h5read(ler_file, '/Wavelength');


xx = gome_lon >= min(lonlim)-buffer & gome_lon <= max(lonlim)+buffer;
yy = gome_lat >= min(latlim)-buffer & gome_lat <= max(latlim)+buffer;

gome_lon = gome_lon(xx);
gome_lat = gome_lat(yy);

% Find the wavelength nearest to 550 (the wavelength the Mobely LUT is
% defined at) and 433 (the middle of the NO2 fitting window, c.f. Krotkov
% et al. 2017, ACPD, doi:10.5194/amt-2017-44).
[~,ind_550] = min(abs(gome_wavelength - 550));
[~,ind_433] = min(abs(gome_wavelength - 433));

ler_at_550 = gome_ler(:, ind_550, xx, yy);
ler_at_433 = gome_ler(:, ind_433, xx, yy);

gome_ler_ratio = permute(squeeze(ler_at_433 ./ ler_at_550), [2 3 1]);
[gome_lat, gome_lon] = meshgrid(gome_lat, gome_lon);

end

