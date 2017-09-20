function [ DataHDF, DataTXT ] = match_hdf_txt( hdffile, txtfile, fields )
%MATCH_HDF_TEXT Matches up data in a text file to the array shape in an HDF file
%   [ DATAHDF, DATATXT ] = MATCH_HDF_TEXT( HDFFILE, TXTFILE, FIELDS )
%   will read data from the files at paths HDFFILE and TXTFILE and return
%   structures DATAHDF and DATATXT with fields FIELDS, matching up the data
%   in the text file with that from the HDF file so that the arrays have
%   the same shape.

fields = [{'Longitude', 'Latitude'}, fields(:)'];

hdf_data = prod_test_load_hdf(hdffile, fields);
hdf_data = cat_hdf(hdf_data);
DataTXT = prod_test_load_txt(txtfile, fields);

def_val = nan(size(DataTXT.Longitude));
DataHDF = make_empty_struct_from_cell(fields, def_val);

for i = 1:numel(DataTXT.Longitude)
    if mod(i,100) == 0
        fprintf('Now on %d of %d\n',i,numel(DataTXT.Longitude))
    end
    hdf_vals = match_pix_lat_lon(DataTXT.Longitude(i), DataTXT.Latitude(i));
    for j = 1:numel(fields)
        DataHDF.(fields{j})(i) = hdf_vals.(fields{j});
    end
end

    function vals = match_pix_lat_lon(lon, lat)
        vals = make_empty_struct_from_cell(fields, nan);
        for a = 1:numel(hdf_data)
            for b = 1:numel(hdf_data(a).Longitude)
                hlon = hdf_data(a).Longitude(b);
                hlat = hdf_data(a).Latitude(b);
                if abs(hlon - lon) < 0.0001 && abs(hlat - lat) < 0.0001
                    for f = 1:numel(fields)
                        vals.(fields{f}) = hdf_data(a).(fields{f})(b);
                        hdf_data.(fields{f})(b) = [];
                    end
                    return
                end
            end
        end
        
        % If we get here, we didn't find a pixel
        warning('Pixel at lon = %f, lat = %f; could not find matching pixel in HDF file', lon, lat);
    end

end

function hdf_data_out = cat_hdf(hdf_data)
fns = fieldnames(hdf_data);
hdf_data_out = make_empty_struct_from_cell(fns);
for f=1:numel(fns)
    tmp = cat_sat_data(hdf_data, fns{f});
    hdf_data_out.(fns{f}) = tmp(:);
end
end