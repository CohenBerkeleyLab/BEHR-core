function [  ] = find_days_to_fix( start_date, end_date )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
data_path = BEHR_paths('behr_mat_dir');
out_path = '.';
%dnums=1:29;
dnums = datenum(start_date):datenum(end_date);
nd=numel(dnums);
spmd
    ndays_per_wk = ceil(nd/numlabs);
    i = (labindex-1)*ndays_per_wk + 1;
    j = min(labindex*ndays_per_wk, nd);
    
    outname = sprintf('behrrun-%d.txt',labindex);
    fid=fopen(fullfile(out_path,outname),'w');
    
    block_on = false;
    block_start = 0;
    
    for a=i:j
        data_name = sprintf('OMI_BEHR_%s.mat',datestr(dnums(a),'yyyymmdd'));
        D = load(fullfile(data_path,data_name),'Data');
        Data = D.Data;
        found = false;
        for s=1:numel(Data)
            if any(Data(s).BEHRAMFTrop(:) < 1e-5)
                if ~block_on
                    block_on = true;
                    block_start = dnums(a);
                end
                found = true;
                break
            end
        end
        if block_on && ~found
            fprintf(fid,'BEHR_main(''%s'',''%s'');\n',datestr(block_start,'yyyy-mm-dd'),datestr(dnums(a-1),'yyyy-mm-dd'));
            block_on = false;
        end
    end
    if block_on
        fprintf(fid,'BEHR_main(''%s'',''%s'');\n',datestr(block_start,'yyyy-mm-dd'),datestr(dnums(j),'yyyy-mm-dd'));
    end
    fclose(fid);
end

end

