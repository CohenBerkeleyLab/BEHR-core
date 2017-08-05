function [  ] = find_and_replace_cloud_fill( start_date, end_date )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
sp_path = behr_paths.sp_mat_dir;
save_path = fullfile(sp_path,'Staging');
for d=datenum(start_date):datenum(end_date)
    fname = sprintf('OMI_SP_%s.mat',datestr(d,'yyyymmdd'));
    if ~exist(fullfile(sp_path,fname),'file')
        continue
    end
    fprintf('%s:\n',fname);
    D = load(fullfile(sp_path,fname),'Data');
    Data = D.Data;
    first_swath = true;
    for s=1:numel(Data)
        if any(Data(s).CloudFraction(:)<0) || any(Data(s).CloudRadianceFraction(:)<0)
            if first_swath
                first_swath = false;
                fprintf('\tFixing %s cloud fraction\n',fname);
            end
            Data(s).CloudFraction(Data(s).CloudFraction==-32.767) = -32767;
            Data(s).CloudRadianceFraction(Data(s).CloudRadianceFraction==-32.767) = -32767;
        end
        Data(s).TerrainReflectivity(Data(s).TerrainReflectivity>-1) = Data(s).TerrainReflectivity(Data(s).TerrainReflectivity>-1)/1000;
    end
    
    save(fullfile(save_path,fname),'Data');
end

end

