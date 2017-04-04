function [ dsetname ] = hdfdsetname( varargin )
%hdfdsetname(info,group1,...,groupN,dataset): a shortcut function for
%returning the full name of an HDF SDS, assuming that the hierarchy is 
%info -> Vgroup(a) -> Vgroup(b) ... -> Vgroup(z)  Takes a minimum of 2 arguments: 
%    info: an object returned by the 'h5info' or an h5 group
%    group1...groupN: The group indicies. Optional if the group is already defined.
%    dataset: The index or name of the dataset of interest
%For example, to list the datasets in hinfo.Group(1).Group(2).Group(1)
%enter 'h5groupdsets(hinfo, 1, 2, 1)
%
% --Josh Laughner <joshlaugh5@gmail.com> 3 Mar 2014--

dset_path = '/';

if length(varargin) < 2
    disp('hdfsetname error: needs at least 2 input arguments')
else
    info = varargin{1}; %Save the h5 object (info or group)
    if length(varargin) > 2 %If there are 3 or more arguments, use all but the first and last to access the child groups
        n = length(varargin) - 2;
        groups = cell2mat(varargin(2:end-1));
        for i=1:n
            info = info.Vgroup(groups(i));
            dset_path = strcat(dset_path, info.Name, '/');
        end
    end
    
    tmp = varargin{end}; 
    if ischar(tmp) %If a name string was passed, take that name and check if it is a dataset name.  Does not stop execution if the name is not found.
        dset = tmp;
        for a = 1:length(info.SDS)
            if strcmp(tmp, info.SDS(a).Name); 
                foundit = 1;
                break; 
            else
                foundit = 0;
            end
        end
        if foundit == 0; disp('Could not verify dataset name'); end
        
    else %Otherwise, assume the input was a dataset index and find the corresponding name
        dset = info.SDS(tmp).Name;
    end
    try
        dsetname = strcat(dset_path, dset); %Format the full dataset path + name in the style needed for h5read.
    catch err
        % If the SDSs are not in a Vgroup and are in the top level, there
        % will be no Name field to read, so format the SDS name to indicate
        % it is in the top level.
        if strcmp(err.identifier,'MATLAB:nonExistentField')
            dsetname = ['/',dset];
        else 
            rethrow(err)
        end
    end
        
end

end

