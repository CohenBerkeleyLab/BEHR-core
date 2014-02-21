function [ dsets ] = h5dsets( h5group )
%h5dsets will iterate through the names of each data set in the group
%passed to the function and return a cell structure with those names

n = length(h5group.Datasets);
dsets = cell(n,1);
for i = 1:n
    info = [num2str(i), ' : ', h5group.Datasets(i).Name];
    dsets(i) = {info};
end

end

