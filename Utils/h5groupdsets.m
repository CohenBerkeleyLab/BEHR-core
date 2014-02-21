function [ dsets ] = h5groupdsets( varargin )
%h5groupdsets(info,group1,...,groupN): a shortcut function for reading the names of h5 datasets,
%assuming that the hierarchy is info -> Group(a) -> Group(b) ... ->
%Group(z).Datasets.  Takes a minimum of 2 arguments: 
%    info: an object returned by the 'h5info' function
%    group1...groupN: The group indicies.  
%For example, to list the datasets in hinfo.Group(1).Group(2).Group(1)
%enter 'h5groupdsets(hinfo, 1, 2, 1)

if length(varargin) < 2
    disp('h5groupdsets error: needs at least 2 input arguments')
else
    info = varargin{1};
    n = length(varargin) - 1;
    groups = cell2mat(varargin(2:end));
    for i=1:n
        info = info.Groups(groups(i));
    end
    m = length(info.Datasets);
    dsets = cell(n,1);
    for i = 1:m
        entry = [num2str(i), ' : ', info.Datasets(i).Name];
        dsets(i) = {entry};
    end        
end

end

