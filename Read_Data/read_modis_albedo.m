function [ data ] = read_modis_albedo( modis_directory, coart_lut, ocean_mask, date_in, data, varargin )
%READ_MODIS_ALBEDO Reads MODIS MCD43C1 BRDF albedo
%   DATA = READ_MODIS_ALBEDO( MODIS_DIR, COART_LUT, OCEAN_MASK, DATE_IN, DATA ) Reads
%   MODIS MCD43C1 data from MODIS_DIR (which must be the path to the root
%   MCD43C1 directory, containing each year in a subfolder). It identifies
%   the proper file to read for the DATE_IN (a date string automatically
%   understood by Matlab or date number), reads in the BRDF kernel
%   coefficients, calculates the kernel values, and combines them to get
%   the surface reflectivity. If a pixel is over water (which means that
%   more than 50% of the coincident MCD43C1 data are fill values), the
%   surface reflectance is deduced from the COART_LUT, the look up table
%   struct returned by COART_SEA_REFLECTANCE. This table must be given as
%   loading the HTML file is problematic in a parallel loop.
%
%   Parameters:
%
%       'DEBUG_LEVEL' - increase the verbosity. Default is 0, higher
%       numbers print more information.
%
%       'LoncornField', 'LatcornField' - change which fields in DATA are
%       used as the definition of the pixel corners
%
% Important references for MODIS BRDF v006 product:
%   V006 User Guide: https://www.umb.edu/spectralmass/terra_aqua_modis/v006
%

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

if ~ischar(modis_directory)
    E.badinput('MODIS_DIRECTORY must be a string')
elseif ~exist(modis_directory, 'dir')
    E.badinput('MODIS_DIRECTORY is not a directory')
end


if isnumeric(date_in)
    if ~isscalar(date_in)
        E.badinput('If given as a number, DATE_IN must be scalar')
    end
elseif ischar(date_in)
    try
        datenum(date_in);
    catch err
        if strcmp(err.identifier, 'MATLAB:datenum:ConvertDateString')
            E.badinput('DATE_IN could not be recognized as a valid format for a date string')
        else
            rethrow(err)
        end
    end
else
    E.badinput('DATE_IN must be a string or number')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

this_year = sprintf('%04d', year(date_in));
alb_dir = fullfile(modis_directory, this_year);
julian_day = modis_date_to_day(date_in);

% As of version 6 of MCD43C1, a 16-day average is produced every day, so
% unlike version 5 of MCD43C3 where we had to look forward and back in time
% from the current date, we should be able to just pick the file for this
% day.

mcd_filename = sprintf('MCD43C1.A%04d%03d*.hdf', year(date_in), julian_day);
alb_filename = fullfile(alb_dir, mcd_filename);
alb_files = dir(alb_filename);
if numel(alb_files) < 1
    E.filenotfound('MODIS BRDF file matching pattern %s.', alb_filename);
elseif numel(alb_files) > 1
    E.toomanyfiles('Multiple MODIS BRDF files found matching pattern %s.', alb_filename);
end

mcd43_info = hdfinfo(fullfile(alb_dir,alb_files(1).name));
band3_iso = hdfreadmodis(mcd43_info.Filename, hdfdsetname(mcd43_info,1,1,'BRDF_Albedo_Parameter1_Band3'));
band3_vol = hdfreadmodis(mcd43_info.Filename, hdfdsetname(mcd43_info,1,1,'BRDF_Albedo_Parameter2_Band3'));
band3_geo = hdfreadmodis(mcd43_info.Filename, hdfdsetname(mcd43_info,1,1,'BRDF_Albedo_Parameter3_Band3'));
brdf_quality = hdfreadmodis(mcd43_info.Filename, hdfdsetname(mcd43_info,1,1,'BRDF_Quality'));

%MODIS albedo is given in 0.05 degree cells and a single file covers the
%full globe, so figure out the lat/lon of the middle of the grid cells as:
[band3_lons, band3_lats] = modis_cmg_latlon(0.05, [-180, 180], [-90, 90], 'grid');

%To speed up processing, restrict the MODIS albedo data to
%only the area we need to worry about.  This will
%significantly speed up the search for albedo values within
%each pixel.
loncorn = data.(loncorn_field)(:);
latcorn = data.(latcorn_field)(:);
loncorn(loncorn==0) = [];
latcorn(latcorn==0) = [];
lon_min = floor(min(loncorn));
lon_max = ceil(max(loncorn));
lat_min = floor(min(latcorn));
lat_max = ceil(max(latcorn));

in_lats = find(band3_lats(:,1)>=lat_min & band3_lats(:,1)<=lat_max);
in_lons = find(band3_lons(1,:)>=lon_min & band3_lons(1,:)<=lon_max);

band3_iso = band3_iso(in_lats,in_lons);
band3_vol = band3_vol(in_lats,in_lons);
band3_geo = band3_geo(in_lats,in_lons);
brdf_quality = brdf_quality(in_lats,in_lons);

band3_lats=band3_lats(in_lats,in_lons);
band3_lons=band3_lons(in_lats,in_lons);

s=size(data.Latitude);
c=numel(data.Latitude);
MODISAlbedo = nan(s);
MODISAlbedoQuality = nan(s);
MODISAlbedoFillFlag = nan(s);
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
    
    xx_mask_lon = ocean_mask.lon(1,:) >= min(xall) & ocean_mask.lon(1,:) <= max(xall);
    xx_mask_lat = ocean_mask.lat(:,1) >= min(yall) & ocean_mask.lat(:,1) <= max(yall);
    
    if DEBUG_LEVEL > 4; fprintf('    Time to cut down ocean mask = %f\n', toc(t_cut)); end
    
    if DEBUG_LEVEL > 4; t_mask = tic; end
    xx_ocean_mask = inpolygon(ocean_mask.lon(xx_mask_lat, xx_mask_lon), ocean_mask.lat(xx_mask_lat, xx_mask_lon), xall, yall);
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
    xx = band3_lons(1,:) >= min(xall) & band3_lons(1,:) <= max(xall);
    yy = band3_lats(:,1) >= min(yall) & band3_lats(:,1) <= max(yall);
    
    band3_iso_k = band3_iso(yy,xx);
    band3_geo_k = band3_geo(yy,xx);
    band3_vol_k = band3_vol(yy,xx);
    brdf_quality_k = brdf_quality(yy,xx);
    
    xx_inpoly = inpolygon(band3_lons(yy,xx),band3_lats(yy,xx),xall,yall);
    
    % Also remove data that has too low a quality. The quality values are
    % described in the "Description" attribute for the "BRDF_Quality" SDS.
    % Lower values for the quality flag are better.
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
data.MODISAlbedoFile = mcd43_info.Filename;
data.AlbedoOceanFlag = ocean_flag;

end

