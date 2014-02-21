function [ dsetname ] = h5dsetname( varargin )
%h5group(info,group1,...,groupN,dataset): a shortcut function for returning an h5group,
%assuming that the hierarchy is info -> Group(a) -> Group(b) ... ->
%Group(z)  Takes a minimum of 3 arguments: 
%    info: an object returned by the 'h5info' or an h5 group
%    group1...groupN: The group indicies. Optional if the group is already defined.
%    dataset: The index or name of the dataset of interest
%For example, to list the datasets in hinfo.Group(1).Group(2).Group(1)
%enter 'h5groupdsets(hinfo, 1, 2, 1)
%
% --Josh Laughner <joshlaugh5@gmail.com> 20 Feb 2014--

if length(varargin) < 2
    disp('h5group error: needs at least 2 input arguments')
else
    info = varargin{1}; %Save the h5 object (info or group)
    if length(varargin) > 2 %If there are 3 or more arguments, use all but the first and last to access the child groups
        n = length(varargin) - 2;
        groups = cell2mat(varargin(2:end-1));
        for i=1:n
            info = info.Groups(groups(i));
        end
    end
    
    tmp = varargin{end}; 
    if ischar(tmp) %If a name string was passed, take that name and check if it is a dataset name.  Does not stop execution if the name is not found.
        dset = tmp;
        for a = 1:length(info.Datasets)
            if strcmp(tmp, info.Datasets(a).Name); 
                foundit = 1;
                break; 
            else
                foundit = 0;
            end
        end
        if foundit == 0; disp('Could not verify dataset name'); end
        
    else %Otherwise, assume the input was a dataset index and find the corresponding name
        dset = info.Datasets(tmp).Name;
    end
    dsetname = [info.Name, '/', dset]; %Format the full dataset path + name in the style needed for h5read.
        
end

end

