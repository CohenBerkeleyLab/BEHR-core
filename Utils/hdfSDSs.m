function [ dsets ] = hdfSDSs( varargin )
%hdfSDSs will iterate through the names of each SDS in the Vgroup
%passed to the function and return a cell structure with those names

if nargin == 1 %If there's only one argument, assume that it is the group we want to enumerate the datasets for
    hdfi = varargin{1};
    
else
    hdfi = varargin{1}; %Save the h5 object (info or group)
    n = length(varargin) - 1;
    groups = cell2mat(varargin(2:end));
    for i=1:n
        hdfi = hdfi.Vgroup(groups(i));
    end
end

n = length(hdfi.SDS);
dsets = cell(n,1);
for i = 1:n
    info = [num2str(i), ' : ', hdfi.SDS(i).Name];
    dsets(i) = {info};
end

end 