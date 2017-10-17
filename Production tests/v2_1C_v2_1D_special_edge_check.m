function [ edge_vals, fns ] = v2_1C_v2_1D_special_edge_check(  )
%v2_1C_v2_1D_special_edge_check Quick check of edge pixel error in v2-1C

vers_D_dir = fullfile(behr_repo_dir, 'Workspaces', 'Production tests', 'v2-1D');
vers_D_files = dir(fullfile(vers_D_dir, 'OMI_BEHR_v2-1D*.mat'));

fns = {'MODISAlbedo','GLOBETerpres','BEHRAMFTrop','BEHRColumnAmountNO2Trop','BEHRAMFTropVisOnly','BEHRColumnAmountNO2TropVisOnly'};
edge_vals.jan_c = cell(size(fns));
edge_vals.jan_d = cell(size(fns));
edge_vals.june_c = cell(size(fns));
edge_vals.june_d = cell(size(fns));

for a=1:numel(vers_D_files)
    
    D = load(fullfile(vers_D_dir, vers_D_files(a).name),'Data');
    cname = strrep(vers_D_files(a).name, 'v2-1D', 'v2-1C');
    C = load(fullfile(BEHR_paths('behr_mat_dir'),cname),'Data');
    
    fdate = datenum(regexp(cname,'\d\d\d\d\d\d\d\d','once','match'),'yyyymmdd');
    
    for b=1:numel(fns)
        data_c = cat_sat_data(C.Data,fns{b});
        data_c = data_c(:,[1 60]);
        data_d = cat_sat_data(D.Data,fns{b});
        data_d = data_d(:,[1 60]);
        
        if month(fdate) == 1
            edge_vals.jan_c{b} = cat(1, edge_vals.jan_c{b}, data_c);
            edge_vals.jan_d{b} = cat(1, edge_vals.jan_d{b}, data_d);
        elseif month(fdate) == 6
            edge_vals.june_c{b} = cat(1, edge_vals.june_c{b}, data_c);
            edge_vals.june_d{b} = cat(1, edge_vals.june_d{b}, data_d);
        else
        end
            
    end
    
end

fprintf('January:\n')
print_results(edge_vals.jan_c, edge_vals.jan_d);

fprintf('\nJune:\n')
print_results(edge_vals.june_c, edge_vals.june_d);

    function print_results(edge_vals_c, edge_vals_d)
        for f=1:numel(fns)
            perdiff = reldiff(edge_vals_d{f}, edge_vals_c{f}, true)*100;
            fprintf('%s %% diff: mean = %f, std = %f\n', fns{f}, nanmean(perdiff), nanstd(perdiff)); 
        end
    end

end

