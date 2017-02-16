function [  ] = find_vcd_fill( start_date, end_date )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
sp_path = BEHR_paths('behr_mat_dir');
%save_path = fullfile(sp_path,'Staging');
for d=datenum(start_date):datenum(end_date)
    fname = sprintf('OMI_BEHR_v2-1Arev1_%s.mat',datestr(d,'yyyymmdd'));
    if ~exist(fullfile(sp_path,fname),'file')
        continue
    end
    D = load(fullfile(sp_path,fname),'Data');
    Data = D.Data;
    dobreak = false;
    fprintf('%s:\n',fname);
    for s=1:numel(Data)
%         if any(Data(s).CloudPressure(:)>1013)
%             fprintf('%s cld P > 1013\n',fname);
%         end
%         if any(Data(s).BEHRColumnAmountNO2Trop(:)<-1e29)
%             %fprintf('\t Data(%d).BEHRColumnAmountNO2Trop < -1e29\n',s);
%         end
        if any(Data(s).BEHRAMFTrop(:)<1e-5)
            fprintf('\t Data(%d).BEHRAMFTrop < 1e-5',s);
            xx = Data(s).BEHRAMFTrop < 1e-5;
            cc = Data(s).CloudPressure > 1013;
            if all((xx(:)&cc(:)) == xx(:))
                fprintf(', all instances due to cld P > 1013 hPa\n');
            else
                fprintf('\n');
            end
        end
        if any(Data(s).ColumnAmountNO2Trop(:)<-1e29)
            fprintf('\t Data(%d).ColumnAmountNO2Trop < -1e29\n',s);
        end
        if any(Data(s).AMFTrop(:)<-30000)
            fprintf('\t Data(%d).AMFTrop < -30000\n',s);
        end
        if any(Data(s).CloudFraction(:)<0)
            fprintf('\t Data(%d).CloudFraction < 0\n',s);
        end
    end
end

end

