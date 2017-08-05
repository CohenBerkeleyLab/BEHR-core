function [ varargout ] = load_behr_file( file_date )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


behr_file = fullfile(behr_paths.behr_mat_dir, behr_filename(file_date));
D = load(behr_file);
if nargout == 0
    Data = D.Data;
    OMI = D.OMI;
    putvar(Data,OMI);
else
    varargout{1} = D.Data;
    varargout{2} = D.OMI;
end

end

