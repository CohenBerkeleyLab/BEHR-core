%modis_cloud_map
%Makes maps of MODIS cloud fraction

date_start = '06/21/2012';
date_end = '06/21/2012';

lonbdy = [-98 -86];
latbdy = [36 43];

DEBUG_LEVEL = 2;

modis_dir = '/Volumes/share-sat/SAT/MODIS/MOD06_L2';
modis_prefix = 'MOD06_L2.A';
modis_field = 'Cloud_Fraction';

%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Check input %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

if lonbdy(1) > lonbdy(2) || latbdy(1) > latbdy(2)
    error('modis_cld_map:latlonbdy','Lat or lon minimum greater than maximum');
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Main Program %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

days = datenum(date_start):datenum(date_end);

lats = latbdy(1):0.01:latbdy(2);
lons = lonbdy(1):0.01:lonbdy(2);
[longrid, latgrid] = meshgrid(lons, lats);
s = size(longrid);
cloudgrid = zeros(s(1), s(2), 10*numel(days));
areaweight = zeros(s(1), s(2), 10*numel(days));
count = zeros(s(1), s(2), 10*numel(days));

D=0;
for d=1:numel(days)
    curr_date = datestr(days(d),29);
    year = curr_date(1:4);
    month = curr_date(6:7);
    day = curr_date(9:10);
    modis_day = num2str(modis_date_to_day(curr_date));
    
    % Find all MODIS files for the current date
    modis_filename=sprintf('%s%s%s*.hdf',modis_prefix,year,modis_day);
    modis_files = dir(fullfile(modis_dir,year,modis_filename));
    n=numel(modis_files);
    % Loop through all files
    for a=1:n
        if DEBUG_LEVEL > 1; fprintf('Checking %s\n',modis_files(a).name); end
        hdfi = hdfinfo(fullfile(modis_dir,year,modis_files(a).name));
        mlat = hdfread(fullfile(modis_dir,year,modis_files(a).name),hdfdsetname(hdfi,1,1,'Latitude'));
        mlon = hdfread(fullfile(modis_dir,year,modis_files(a).name),hdfdsetname(hdfi,1,1,'Longitude'));
        if all(mlat(:) < latbdy(1)) || all(mlat(:) > latbdy(2)) || all(mlon(:) < lonbdy(1)) || all(mlon(:) > lonbdy(2))
            continue % If it falls outside our lat/lon boundaries, skip the file. 
        else
            D=D+1;
            % Calculate corner points based on the geometry of the pixels
            corners = geometric_corner_calc(mlat, mlon);
            latcorn = corners(:,:,2,1:4);
            loncorn = corners(:,:,1,1:4);
            
            % Load cloud fraction
            cldfrac = hdfread(fullfile(modis_dir,year,modis_files(a).name),hdfdsetname(hdfi,1,2,modis_field));
            cldfrac = double(cldfrac);
            
            for b=1:size(mlat,1)
                for c=1:size(mlat,2)
                    % Find all oversampled pixels from the grid that fall within the boundaries
                    % of this modis pixel
                    latcorn_bc = squeeze(latcorn(b,c,:)); loncorn_bc = squeeze(loncorn(b,c,:));
                    yy = find(lats > min(latcorn(b,c,:)) & lats < max(latcorn(b,c,:)));
                    xx = find(lons > min(loncorn(b,c,:)) & lons < max(loncorn(b,c,:)));
                    sub_latgrid = latgrid(yy,xx); sub_longrid = longrid(yy,xx);
                    
                    latvert = [latcorn_bc; latcorn_bc(1)]; lonvert = [loncorn_bc; loncorn_bc(1)];
                    IN = inpolygon(sub_longrid, sub_latgrid, lonvert, latvert);
                    [yyg,xxg] = meshgrid(yy,xx);
                    
                    % Calculate the areaweight for this modis pixel
                    area = polyarea(lonvert,latvert);
                    
                    % Assign cloud fraction and area weight values to the
                    % correct index of the oversampled grid
                    
                    cloudgrid(yyg(IN),xxg(IN),D) = (cloudgrid(yyg(IN),xxg(IN),D) .* count(yyg(IN),xxg(IN),D) + cldfrac(b,c))./(count(yyg(IN),xxg(IN),D)+1);
                    areaweight(yyg(IN),xxg(IN),D) = areaweight(yyg(IN),xxg(IN),D) + 1/area;
                    count(yyg(IN),xxg(IN),D) = count(yyg(IN),xxg(IN),D) + 1;
                end
            end
        end
    end
    
    % Remove the extra entried in the 3rd dimension for the cloud fraction
    % and areaweight matrices
    cloudgrid = cloudgrid(:,:,1:D);
    areaweight = areaweight(:,:,1:D);
    weighted_cloudfrac = sum(cloudgrid .* areaweight,3);
    total_weight = sum(areaweight,3);
    
    weighted_cloudfrac = weighted_cloudfrac ./ total_weight;
    weighted_cloudfrac(total_weight==0) = NaN;
    
    % Map it
    figure;
    m_proj('Albers Equal-Area Conic','lon',lonbdy,'lat',latbdy);
    m_pcolor(lons,lats,weighted_cloudfrac); shading flat
    m_states('w');
    m_grid('linestyle','none');
    cb = colorbar;
end