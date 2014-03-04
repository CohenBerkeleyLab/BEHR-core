function [ dsetID ] = hdf_dsetID( varargin )
%hdf_dsetID(hdfinfo, Vgroup1, Vgroup2,...,VgroupN, SDS) is a function that
%returns an SDS (Scientific Data Set) from an HDF file that can then be
%passed to hdfread().  The SDS can be identified by its index or name, so
%for example if Longitude can be found in the file held in the hdfinfo
%object 'hdfi', as hdfi.Vgroup(1).Vgroup(1).SDS(2), then both
%"hdf_dsetID(hdfi,1,1,2)" or "hdf_dsetID(hdfi,1,1,'Latitude')" will return
%the proper SDS.
%
% --Josh Laughner <joshlaugh5@gmail.com> 3 Mar 2014--

if length(varargin) < 2
    disp('hdfsetname error: needs at least 2 input arguments')
else
    info = varargin{1}; %Save the h5 object (info or group)
    if length(varargin) > 2 %If there are 3 or more arguments, use all but the first and last to access the child groups
        n = length(varargin) - 2;
        groups = cell2mat(varargin(2:end-1));
        for i=1:n
            info = info.Vgroup(groups(i));
        end
    end
    
    tmp = varargin{end}; 
    if ischar(tmp) %If a name string was passed, take that name and check if it is a dataset name.  Does not stop execution if the name is not found.
        dset = tmp;
        for a = 1:length(info.SDS)
            if strcmp(tmp, info.SDS(a).Name); 
                SDSnumber = a;
                foundit = 1;
                break; 
            else
                foundit = 0;
            end
        end
        if foundit == 0; disp('Could not verify dataset name'); end
        
    else %Otherwise, assume the input was an SDS index
        SDSnumber = tmp;
    end
    dsetID = info.SDS(SDSnumber); %Format the full dataset path + name in the style needed for h5read.
        
end

end

