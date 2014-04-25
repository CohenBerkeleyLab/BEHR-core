function [ Data ] = behrhdf2struct( hdfi )
%Parses a BEHR HDF file (actually an HDF5 or he5 file) into a Matlab data structure
%   Useful for comparing old values in BEHR HDF files to current Matlab
%   data structures. Especially useful in conjunction with
%   "Compare_Data_Fields", which allows you to compare an aribitrary field
%   in two data structures
%
%   Josh Laughner <joshlaugh5@gmail.com> 16 Apr 2014

Data=struct([]);
n = length(hdfi.Groups(1).Groups);
for a=1:n
    
    k = strfind(hdfi.Groups(1).Groups(a).Name,'Swath');
    swath = str2double(hdfi.Groups(1).Groups(a).Name(k+5:end));
    Data(a).Swath = swath;
    
    for b=1:length(hdfi.Groups(1).Groups(a).Datasets)
        dset = h5read(hdfi.Filename, h5dsetname(hdfi,1,a,b));
        eval(sprintf('Data(a).%s = dset;',hdfi.Groups(1).Groups(a).Datasets(b).Name));
    end
    
end

end

