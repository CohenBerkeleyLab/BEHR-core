function [ data ] = read_modis_albedo( modis_directory, date_in, data, varargin )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% This function will likely change significantly when we merge the BRDF
% albedo branch into the update branch

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
p = inputParser;
p.addParameter('DEBUG_LEVEL', 0);
p.addParameter('LoncornField', 'FoV75CornerLongitude');
p.addParameter('LatcornField', 'FoV75CornerLatitude');
p.parse(varargin{:});
pout = p.Results;

DEBUG_LEVEL = pout.DEBUG_LEVEL;
loncorn_field = pout.LoncornField;
latcorn_field = pout.LatcornField;

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

%Find the closest MCD file
in=[0 1 -1 2 -2 3 -3 4 -4 5 -5 6 -6 7 -7 8 -8 9 -9 10 -10 11 -11 12 -12 13 -13 14 -14 15 -15 16 -16 17 -17 18 -18 19 -19 20 -20 21 -21];
for a=1:length(in);
    mcd_filename = sprintf('MCD43C3.A%04d%03d*.hdf', year(date_in), julian_day + in(a));
    alb_filename = fullfile(alb_dir, mcd_filename);
    alb_files = dir(alb_filename);
    if DEBUG_LEVEL > 1; fprintf('Looking for %s \n', alb_filename); end
    %if exist('alb_filename','file') == 2
    if ~isempty(alb_files)
        if DEBUG_LEVEL > 0; fprintf(' Found mcd43 file %s \n', alb_filename); end
        alb_file_full = fullfile(alb_dir, alb_files(1).name);
        break
    elseif a==length(in)
        error(E.filenotfound('MCD43C3 (Albedo) file within 21 days'));
    end
end

mcd43_info = hdfinfo(fullfile(alb_dir,alb_files(1).name));
band3 = hdfreadmodis(alb_file_full, hdfdsetname(mcd43_info, 1, 1, 'Albedo_BSA_Band3'));
band3 = flipud(band3);

%MODIS albedo is given in 0.05 degree cells and a single file covers the
%full globe, so figure out the lat/lon of the middle of the grid cells as:
band3_lat=-90+0.05/2:0.05:90-0.05/2; band3_lats=band3_lat'; band3_lats=repmat(band3_lats,1,7200);
band3_lon=-180+0.05/2:0.05:180-0.05/2; band3_lons=repmat(band3_lon,3600,1);

%To speed up processing, restrict the MODIS albedo data to
%only the area we need to worry about.  This will
%significantly speed up the search for albedo values within
%each pixel.
loncorn = data.(loncorn_field);
latcorn = data.(latcorn_field);
loncorn(loncorn==0) = [];
latcorn(latcorn==0) = [];
lon_min = floor(min(loncorn));
lon_max = ceil(max(loncorn));
lat_min = floor(min(latcorn));
lat_max = ceil(max(latcorn));

in_lats = find(band3_lat>=lat_min & band3_lat<=lat_max);
in_lons = find(band3_lon>=lon_min & band3_lon<=lon_max);
band3=band3(in_lats,in_lons);
band3_lats=band3_lats(in_lats,in_lons);
band3_lons=band3_lons(in_lats,in_lons);
s=size(data.Latitude);
c=numel(data.Latitude);
MODISAlbedo=zeros(s);

%Now actually average the MODIS albedo for each OMI pixel
if DEBUG_LEVEL > 0; disp(' Averaging MODIS albedo to OMI pixels'); end
for k=1:c;
    if DEBUG_LEVEL > 2; tic; end
    
    xall=[data.(loncorn_field)(:,k); data.(loncorn_field)(1,k)];
    yall=[data.(latcorn_field)(:,k); data.(latcorn_field)(1,k)];
    
    % should be able to speed this up by first restricting based on a
    % single lat and lon vector
    xx_alb = inpolygon(band3_lats,band3_lons,yall,xall);
    
    band3_vals=band3(xx_alb);  band3_zeros=band3_vals==0;
    band3_vals(band3_zeros)=NaN; band3_vals(isnan(band3_vals))=[];
    band3_avg=mean(band3_vals);
    
    %put in ocean surface albedo from LUT
    if isnan(band3_avg)==1;
        sza_vec = [5 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 89];
        alb_vec = [0.038 0.038 0.039 0.039 0.040 0.042 0.044 0.046 0.051 0.058 0.068 0.082 0.101 0.125 0.149 0.158 0.123 0.073];
        alb = interp1(sza_vec, alb_vec, data.SolarZenithAngle(k));
        band3_avg = alb;
    end
    
    MODISAlbedo(k) = band3_avg;
    if DEBUG_LEVEL > 2; telap = toc; fprintf(' Time for MODIS alb --> pixel %u/%u = %g sec \n',k,c,telap); end
end

data.MODISAlbedo = MODISAlbedo;
data.MODISAlbedoFile = alb_file_full;

end

