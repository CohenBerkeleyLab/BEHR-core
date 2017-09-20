function [  ] = find_and_replace_cloud_fill( start_date, end_date )
%FIND_AND_REPLACE_CLOUD_FILL Fixes mistakenly scaled cloud fill values
%   In the SP data, the CloudFraction and CloudRadianceFraction fields are
%   given as integers that are to be scaled by 1000 in order to give the
%   true value between 0 and 1. In an early version of BEHR, these fields
%   were scaled before the fill values were treated, mistakenly scaling the
%   fill values down by 1000 to. While not a fatal error, I wrote this
%   function to replace the incorrectly scaled fill values with the
%   unscaled ones.
%
%   FIND_AND_REPLACE_CLOUD_FILL( START_DATE, END_DATE ) loads OMI_SP .mat
%   files from the sp_mat_dir specified by behr_paths, corrects the fill
%   values, and saves them in a subdirectory of the sp_mat_dir, "Staging".
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

