function [ data ] = avg_modis_alb_to_pixels( band3data, coart_lut, ocean_mask, data, varargin )
%AVG_MODIS_ALB_TO_PIXELS Calculate surface reflectivity from MODIS BRDFs
%   DATA = AVG_MODIS_ALB_TO_PIXELS( BAND3DATA, COART_LUT, OCEAN_MASK, DATA
%   ) Handles calculating surface reflectivity from MODIS BRDF kernels and
%   averaging the result to OMI pixels. BAND3DATA must be the structure
%   returned from READ_MODIS_ALBEDO for the correct day. COART_LUT is the
%   look up table of ocean reflectivity returned as the second output of
%   COART_SEA_REFLECTANCE. OCEAN_MASK must be a structure containing the
%   fields "mask", "lon", and "lat" that define a boolean mask for ocean
%   (true for ocean). DATA is the structure containing native OMI pixel
%   data. This function returns that structure with fields MODISAlbedo,
%   MODISAlbedoQuality, MODISAlbedoFillFlag, MODISAlbedoFile, and
%   AlbedoOceanFlag added.
%
%   Additional parameters:
%
%   'DEBUG_LEVEL' - scalar number that controls level of output to
%   terminal. 0 (default) is none, high values give more output.
%
%   'LoncornField' - string that indicates which field of Data to use for
%   longitude corners of pixels. Default is 'FoV75CornerLongitude'.
%
%   'LatcornField' - string that indicates which field of Data to use for
%   latitude corners of pixels. Default is 'FoV75CornerLatitude'.
%
%   'QualityLimit' - controls which MODIS BRDF grid cells will be used by
%   filtering for quality. In MCD43D31, 0 = best quality and 3 = low
%   quality. Only data with quality <= QualityLimit will be used. Default
%   is Inf, i.e. all data is used.

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
p = inputParser;
p.addParameter('DEBUG_LEVEL', 0);
p.addParameter('LoncornField', 'FoV75CornerLongitude');
p.addParameter('LatcornField', 'FoV75CornerLatitude');
p.addParameter('QualityLimit', Inf);

p.parse(varargin{:});
pout = p.Results;

DEBUG_LEVEL = pout.DEBUG_LEVEL;
loncorn_field = pout.LoncornField;
latcorn_field = pout.LatcornField;
max_qual_flag = pout.QualityLimit;


s=size(data.Latitude);
c=numel(data.Latitude);
MODISAlbedo = nan(s);
MODISAlbedoQuality = nan(s);
MODISAlbedoFillFlag = false(s);
ocean_flag = false(s);

%Now actually average the MODIS albedo for each OMI pixel
if DEBUG_LEVEL > 0; disp(' Averaging MODIS albedo to OMI pixels'); end

for k=1:c;
    if DEBUG_LEVEL > 3; t_total=tic; end
    
    xall=[data.(loncorn_field)(:,k); data.(loncorn_field)(1,k)];
    yall=[data.(latcorn_field)(:,k); data.(latcorn_field)(1,k)];
    
    % If there is an invalid corner coordinate, skip because we cannot
    % be sure the correct polygon will be used.
    if any(isnan(xall)) || any(isnan(yall))
        continue
    end
    
    % Next, check if we are over ocean using the ocean mask. If the mask
    % indicates that more than 50% of the pixel is ocean, then we will
    % insert a value from the look up table and move on.
    if DEBUG_LEVEL > 4; t_cut = tic; end
    
    xx_mask_lon = ocean_mask.lon >= min(xall) & ocean_mask.lon <= max(xall);
    xx_mask_lat = ocean_mask.lat >= min(yall) & ocean_mask.lat <= max(yall);
    
    if DEBUG_LEVEL > 4; fprintf('    Time to cut down ocean mask = %f\n', toc(t_cut)); end
    
    if DEBUG_LEVEL > 4; t_mask = tic; end
    [om_longrid, om_latgrid] = latlon_vec2grid(ocean_mask.lon, ocean_mask.lat, xx_mask_lon, xx_mask_lat);
    xx_ocean_mask = inpolygon(om_longrid, om_latgrid, xall, yall);
    if DEBUG_LEVEL > 4; fprintf('    Time to apply inpolygon to mask = %f\n', toc(t_mask)); end
    
    if DEBUG_LEVEL > 4; t_avg_mask = tic; end
    sub_mask = ocean_mask.mask(xx_mask_lat, xx_mask_lon);
    avg_mask = nanmean(sub_mask(xx_ocean_mask));
    
    if DEBUG_LEVEL > 4; fprintf('    Time to average ocean mask = %f\n', toc(t_avg_mask)); end
    
    if avg_mask > 0.5
        if DEBUG_LEVEL > 4; t_ocean = tic; end
        MODISAlbedo(k) = coart_sea_reflectance(data.SolarZenithAngle(k), coart_lut);
        ocean_flag(k) = true;
        if DEBUG_LEVEL > 4; fprintf('    Time to look up ocean reflectance = %f\n', toc(t_ocean)); end
        if DEBUG_LEVEL > 3; telap = toc(t_total); fprintf(' Time for MODIS alb --> pixel %u/%u = %g sec \n',k,c,telap); end
        continue
    end
    
    if DEBUG_LEVEL > 4; t_polygon = tic; end
    
    % If we're here, we're over a land pixel.
    % should be able to speed this up by first restricting based on a
    % single lat and lon vector
    xx = band3data.lons >= min(xall) & band3data.lons <= max(xall);
    yy = band3data.lats >= min(yall) & band3data.lats <= max(yall);
    
    % Checking here became necessary after we switched to using lon/lat
    % vectors instead of arrays because if only one of xx or yy is all
    % false, then brdf_quality_k will be empty but have non-zero "length"
    % in one dimension while xx_inpoly will just be a regular empty array.
    % This causes the calculation of xx_alb to fail, since xx_inpoly and
    % brdf_quality_k have different dimensions (but both are empty).
    if all(~xx) || all(~yy)
        MODISAlbedo(k) = nan;
        if DEBUG_LEVEL > 3; telap = toc(t_total); fprintf(' Time for MODIS alb --> pixel %u/%u = %g sec \n',k,c,telap); end
        continue
    end
    
    band3_iso_k = band3data.iso(yy,xx);
    band3_geo_k = band3data.geo(yy,xx);
    band3_vol_k = band3data.vol(yy,xx);
    brdf_quality_k = band3data.quality(yy,xx);

    [lon_grid, lat_grid] = latlon_vec2grid(band3data.lons, band3data.lats, xx, yy);
    xx_inpoly = inpolygon(lon_grid,lat_grid,xall,yall);
    
    % Also remove data that has too low a quality. The quality values are
    % described in the "Description" attribute for the "BRDF_Quality" SDS.
    % Lower values for the quality flag are better. Inequality comparisons
    % with NaNs are always false, but this explicitly rejects BRDF values
    % for which the quality value is a NaN (i.e. fill value).
    xx_alb = xx_inpoly & brdf_quality_k <= max_qual_flag & ~isnan(brdf_quality_k);
    
    if sum(xx_alb) == 0
        MODISAlbedo(k) = nan;
        if DEBUG_LEVEL > 3; telap = toc(t_total); fprintf(' Time for MODIS alb --> pixel %u/%u = %g sec \n',k,c,telap); end
        continue
    end
    
    if DEBUG_LEVEL > 4; fprintf('    Time to identify MODIS albedo in OMI pixel = %f\n', toc(t_polygon)); end
    
    % The 180-RAA should flip the RAA back to the standard definition (i.e.
    % the supplemental angle of what's in the data product). See the help
    % text for modis_brdf_kernels for why that matters.
    if DEBUG_LEVEL > 4; t_kernels = tic; end
    band3_vals = modis_brdf_alb(band3_iso_k(xx_alb), band3_vol_k(xx_alb), band3_geo_k(xx_alb), data.SolarZenithAngle(k), data.ViewingZenithAngle(k), 180-data.RelativeAzimuthAngle(k));
    if DEBUG_LEVEL > 4; fprintf('    Time to calculate BRDF albedo = %f\n', toc(t_kernels)); end
    
    
    % According to the MOD43 TBD
    % (https://modis.gsfc.nasa.gov/data/atbd/atbd_mod09.pdf, p. 32) the
    % Ross-Li kernel occasionally produces slightly negative albedos. In
    % practice, I have seen negative values around 1e-4 for individual
    % elements of band3_vals. Since this is apparently expected, I will
    % keep the negative values in for the average (which both avoids any
    % potential problem with biasing the albedos high and shouldn't change
    % the albedo much on average) but if the average itself is negative, we
    % reject it and insert NaN as a fill value, which should prevent
    % retrieving that pixel.
    band3_avg = nanmean(band3_vals(:));
    if band3_avg < 0
        warning('Negative average albedo detected. Setting albedo to NaN.');
        % Although we initialized these as NaNs, forcing this to NaN here
        % ensures that no changes to the initialization mess this up
        MODISAlbedo(k) = NaN;
        MODISAlbedoQuality(k) = NaN;
    else
        MODISAlbedo(k) = band3_avg;
        MODISAlbedoQuality(k) = nanmean(brdf_quality_k(xx_alb));
    end
    
    % If more than 50% of the quality values are fills, set the fill
    % warning flag. This will be used in behr_quality_flags to warn of low
    % quality MODIS data.
    if sum(isnan(brdf_quality_k(xx_inpoly)))/sum(xx_inpoly(:)) > 0.5
        MODISAlbedoFillFlag(k) = true;
    end
    
    if DEBUG_LEVEL > 3; telap = toc(t_total); fprintf(' Time for MODIS alb --> pixel %u/%u = %g sec \n',k,c,telap); end
end

data.MODISAlbedo = MODISAlbedo;
data.MODISAlbedoQuality = MODISAlbedoQuality;
data.MODISAlbedoFillFlag = MODISAlbedoFillFlag;
data.MODISAlbedoFile = band3data.files;
data.AlbedoOceanFlag = ocean_flag;

end

function [longrid, latgrid] = latlon_vec2grid(lonvec, latvec, xx_lon, xx_lat)
[longrid, latgrid] = meshgrid(lonvec(xx_lon), latvec(xx_lat));
end