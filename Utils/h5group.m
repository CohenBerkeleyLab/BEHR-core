function [ info ] = h5group( varargin )
%h5group(info,group1,...,groupN): a shortcut function for returning an h5group,
%assuming that the hierarchy is info -> Group(a) -> Group(b) ... ->
%Group(z)  Takes a minimum of 2 arguments: 
%    info: an object returned by the 'h5info' function
%    group1...groupN: The group indicies.  
%For example, to get hinfo.Group(1).Group(2).Group(1)
%enter "h5group(hinfo, 1, 2, 1)"

if length(varargin) < 2
    disp('h5group error: needs at least 2 input arguments')
else
    info = varargin{1};
    n = length(varargin) - 1;
    groups = cell2mat(varargin(2:end));
    for i=1:n
        info = info.Groups(groups(i));
    end
end

end

