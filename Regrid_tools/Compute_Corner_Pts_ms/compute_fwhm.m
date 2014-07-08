%%compute_fwhm
%%arr 12/10/2007

function [fwhm] = compute_fwhm(lat, lon, satlat, satlon, satalt)

n = length(lat); 

fwhm  = zeros(n,1);
distance = zeros(n,1);

for i = 1:n;
    distance(i) = compute_distance(lat(i), lon(i), satlat, satlon, satalt); %This is the distance from the satellite to the lat/lon specified
    if isnan(distance(i));
        fwhm(i) = NaN;
    else
        fwhm(i) = weight_distance(distance(i));
    end
end

