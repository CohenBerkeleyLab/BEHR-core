function [ in ] = is_in_pixel( lon, lat, pix_loncorn, pix_latcorn )
%IN = IS_IN_PIXEL( LON, LAT, PIX_LONCORN, PIX_LATCORN )
%   Returns a logical array IN the same size as LON and LAT that is true
%   for values of LON and LAT that are inside the polygon specified by
%   PIX_LONCORN and PIX_LATCORN.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
if ndims(lon) ~= ndims(lat) || any(size(lon) ~= size(lat))
    E.badinput('LON and LAT must be the same size')
end
if numel(pix_loncorn) ~= 4 || numel(pix_latcorn) ~= 4
    E.badinput('PIX_LONCORN and PIXLATCORN must 4 elements')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% ensure corners do not criss cross (i.e. the first is not catty-corner to
% 4 or 2.
[pix_loncorn, pix_latcorn] = uncross_pix_corners(pix_loncorn, pix_latcorn);

in = inpolygon(lon, lat, pix_loncorn, pix_latcorn);

end

