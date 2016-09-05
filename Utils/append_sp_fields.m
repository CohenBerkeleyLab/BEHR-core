function [ ] = append_sp_fields( sp_dir, behr_dir, out_dir )
%APPEND_SP_FIELDS Appends new SP fields to BEHR mat files
%   APPEND_SP_FILEDS( SP_DIR, BEHR_DIR, OUT_DIR ) This will load the SP
%   .mat files located in SP_DIR, find the corresponding BEHR file located
%   in BEHR_DIR, append any new fields in the Data structure in the SP
%   files to the Data structure in the BEHR file and save the resulting
%   file in OUT_DIR.
%
%   Currently this will simply do all files located in SP_DIR, it has no
%   mechanism for creating a temporal subset of those files.

if strcmp(behr_dir, out_dir)
    error('BEHR_DIR and OUT_DIR should not be the same')
end


SP = dir(fullfile(sp_dir, '.mat'));

for a=1:numel(SP)
    S = load(fullfile(sp_dir, SP(a).name));
    dstr = regexp(SP(a).name, '\d\d\d\d\d\d\d\d', 'match', 'once');
    behr_name = sprintf('OMI_BEHR_%s_%s.mat', BEHR_version, dstr);
    B = load(fullfile(behr_dir, behr_name));
    
    if numel(S.Data) ~= numel(B.Data)
        error('SP and BEHR data structures for %s have inconsistent numbers of swaths', dstr)
    end
    
    fns = fieldnames(S.Data);
    xx = ~ismember(fns, fieldnames(B.Data));
    new_fields = fns(xx);
    for b=1:numel(S.Data)
        for c=1:numel(new_fields)
            B.Data.(new_fields{c}) = S.Data.(new_fields{c});
        end
    end
    
    Data = B.Data;
    OMI = B.OMI;
    
    save(fullfile(out_dir, behr_name), 'Data', 'OMI');
end


end

