function [ dsets ] = h5dsets( varargin )
%h5dsets will iterate through the names of each data set in the group
%passed to the function and return a cell structure with those names

if nargin == 1 %If there's only one argument, assume that it is the group we want to enumerate the datasets for
    h5group = varargin{1};
    
else
    h5group = varargin{1}; %Save the h5 object (info or group)
    n = length(varargin) - 1;
    groups = cell2mat(varargin(2:end));
    for i=1:n
        h5group = h5group.Groups(groups(i));
    end
end

n = length(h5group.Datasets);
dsets = cell(n,1);
for i = 1:n
    info = [num2str(i), ' : ', h5group.Datasets(i).Name];
    dsets(i) = {info};
end

end 