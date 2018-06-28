function [ band3data ] = read_modis_albedo( modis_directory, date_in, lonlim, latlim, varargin )
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
%       'band3data' - supply the days MCD43D Band 3 parameters so that they
%       don't have to be read in again.
%
% Important references for MODIS BRDF v006 product:
%   V006 User Guide: https://www.umb.edu/spectralmass/terra_aqua_modis/v006
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
p = inputParser;
p.addParameter('DEBUG_LEVEL', 2);

p.parse(varargin{:});
pout = p.Results;

DEBUG_LEVEL = pout.DEBUG_LEVEL;

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

%To speed up processing, restrict the MODIS albedo data to only the area we
%need to worry about.  This will significantly speed up the search for
%albedo values within each pixel, plus speed up loading, since fewer IO
%reads are necessary.
[band3_lons, band3_lats, in_lons, in_lats] = modis_cmg_latlon(1/120, lonlim, latlim);
% Make band3_lats a column vector since the first dim of the modis arrays
% is latitude
band3_lats = band3_lats';

% As of version 6 of MCD43, a 16-day average is produced every day, so
% unlike version 5 of MCD43 where we had to look forward and back in time
% from the current date, we should be able to just pick the file for this
% day.

if DEBUG_LEVEL > 0
    fprintf('Reading MODIS BRDF data\n');
end

% Store the four MCD43D files used for traceability
modis_files = cell(1,4);

if DEBUG_LEVEL > 2; fprinf('    Reading band 3 f_iso\n'); end
mcd_filename = sprintf('MCD43D07.A%04d%03d*.hdf', year(date_in), julian_day);
[band3_iso, modis_files{1}] = read_band_parameter(alb_dir, mcd_filename, {in_lats, in_lons});

if DEBUG_LEVEL > 2; fprinf('    Reading band 3 f_vol\n'); end
mcd_filename = sprintf('MCD43D08.A%04d%03d*.hdf', year(date_in), julian_day);
[band3_vol, modis_files{2}] = read_band_parameter(alb_dir, mcd_filename, {in_lats, in_lons});

if DEBUG_LEVEL > 2; fprinf('    Reading band 3 f_geo\n'); end
mcd_filename = sprintf('MCD43D09.A%04d%03d*.hdf', year(date_in), julian_day);
[band3_geo, modis_files{3}] = read_band_parameter(alb_dir, mcd_filename, {in_lats, in_lons});

% Unlike the parameter files, the quality file has all seven bands in one
% file, so we need to handle it differently
mcd_filename = sprintf('MCD43D31.A%04d%03d*.hdf', year(date_in), julian_day);
alb_filename = fullfile(alb_dir, mcd_filename);
alb_files = dir(alb_filename);
if numel(alb_files) < 1
    E.filenotfound('MODIS BRDF file matching pattern %s.', alb_filename);
elseif numel(alb_files) > 1
    E.toomanyfiles('Multiple MODIS BRDF files found matching pattern %s.', alb_filename);
end
modis_files{4} = fullfile(alb_dir, alb_files(1).name);
mcd43_info = hdfinfo(modis_files{4});
brdf_quality = hdfreadmodis(modis_files{4}, hdfdsetname(mcd43_info,4,1,'BRDF_Albedo_Band_Quality_Band3'), 'log_index', {in_lats, in_lons});

% Verify that fill are the same in the 3 parameters and the quality flags.
% This assumption is used in avg_modis_alb_to_pixels in order to remove
% fill value BRDF coefficients and flag OMI pixels where >50% of the MODIS
% data is fill values.
qual_nans = isnan(brdf_quality(:));
if any(xor(qual_nans, isnan(band3_iso(:)))) || any(xor(qual_nans, isnan(band3_geo(:)))) || any(xor(qual_nans, isnan(band3_vol(:))))
    E.callError('inconsistent_fills', 'Fill values are not the same in the quality flags and one or more of the BRDF parameters');
end

band3data.lons = band3_lons;
band3data.lats = band3_lats;
band3data.iso = band3_iso;
band3data.geo = band3_geo;
band3data.vol = band3_vol;
band3data.quality = brdf_quality;
band3data.files = modis_files;
end

function [band3_param, mcd_filename] = read_band_parameter(file_dir, file_pattern, logical_indices)
E = JLLErrors;

alb_filename = fullfile(file_dir, file_pattern);
alb_files = dir(alb_filename);
if numel(alb_files) < 1
    E.filenotfound('MODIS BRDF file matching pattern %s.', alb_filename);
elseif numel(alb_files) > 1
    E.toomanyfiles('Multiple MODIS BRDF files found matching pattern %s.', alb_filename);
end

% Each MCD43D file has only a single SDS representing a single BRDF
% parameter in a single band.
mcd43_info = hdfinfo(fullfile(file_dir,alb_files(1).name));
if numel(mcd43_info.Vgroup(1).Vgroup(1).SDS) ~= 1
    E.callError('mcd43d format', 'READ_BAND_PARAMETER assumes there is only a single SDS in the file; that is not true in %s', mcd43_info.Filename);
end
band3_param = hdfreadmodis(mcd43_info.Filename, hdfdsetname(mcd43_info,1,1,1), 'log_index', logical_indices);
mcd_filename = mcd43_info.Filename;
end

