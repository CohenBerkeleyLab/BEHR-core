function [ data, there_are_points ] = read_omi_sp( sp_file, sp_group_path, sp_vars, data, lonlim, latlim, varargin )
%READ_OMI_SP Reads in an OMI Standard Product data file
%   DATA = READ_OMI_SP( SP_FILE, SP_VARS, DATA ) Reads in a NASA OMI .he5
%   (HDF version 5) file at the path SP_FILE. It will read in the variables
%   specified in the cell array SP_VARS and store them in the structure
%   DATA. DATA must be scalar and must have each variable name given in
%   SP_VARS as a field, plus the field Row.
%
%   Parameters:
%       dim_order - A string indicating how the dimensions should be
%       ordered, default is 'olx'. 'x' indicates across track, 'l' along
%       track, and 'o' other dimensions.
%
%       match_data - will match new data to the existing DATA structure
%       based on the longitude and latitude points (specifically, uses
%       find_submatrix2 to identify the subset of read in points that match
%       to within a tolerance of the existing Longitude and Latitude fields
%       in DATA

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
p = inputParser;
p.addParameter('dim_order', 'olx');
p.addParameter('match_data', false);
p.parse(varargin{:});
pout = p.Results;

dim_order = pout.dim_order;
match_data = pout.match_data;

if ~ischar(sp_file)
    E.badinput('SP_FILE must be a string')
elseif ~exist(sp_file, 'file')
    E.filenotfound('Could not find %s', sp_file)
end

if ~iscellstr(sp_vars)
    E.badinput('SP_VARS must be a cell array of strings')
end

if ~ischar(dim_order) || any(~ismember('olx',dim_order)) || length(dim_order) ~= 3
    E.badinput('The parameter DIM_ORDER must be a string consisting of the characters o, l, and x only')
end

if ~isscalar(match_data) || (~isnumeric(match_data) && ~islogical(match_data))
    E.badinput('The parameter MATCH_DATA must be a scalar boolean or number');
elseif match_data && ((isscalar(data.Longitude) && data.Longitude == 0) || (isscalar(data.Latitude) && data.Latitude ==0))
    warning('MATCH_DATA requested but latitude/longitude fields in DATA do not appear to be fully filled (are scalars and == 0)');
end

% By requiring that all fields are already present in the data structure,
% we ensure that the order will remain the same as we add data to it. If we
% added the fields as the variables came up, that might not be true.
if ~isstruct(data) || ~isscalar(data)
    E.badinput('DATA must be a scalar structure')
end
xx = ~isfield(data, sp_vars);
if any(xx)
    E.badinput('All variable names in SP_VARS must exist as fields in DATA (missing: %s)', strjoin(sp_vars(xx)', ', '));
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

hgrp_info = h5info(sp_file, sp_group_path);

%Read in the full latitude data set; this will be used to determine
%which pixels to read in later.
lat = double(h5readomi(sp_file, find_dset(hgrp_info, 'Latitude')));
lat = lat';
lon = double(h5readomi(sp_file, find_dset(hgrp_info, 'Longitude')));
lon = lon';

%Restrict data to that which falls within the bounds specified by the lon
%and lat limits
if ~match_data
    xx = lon > lonlim(1) & lon < lonlim(2);
    yy = lat > latlim(1) & lat < latlim(2);
else
    % NaNs mess up find_submatrix2. So what we need to do is subset out the
    % non-NaN components and find the submatrix just in there, then
    % translate that back into the full matrix. To do this, we're going to
    % assume that if any coordinates are NaNs, that an entire across-track
    % row is NaNs, which happens when doing zoom-mode days. If this
    % assumption is violated, or if there are differing numbers of NaNs,
    % abort.
    data_nans = all(isnan(data.Longitude),1);
    if any(xor(isnan(data.Longitude(:)), isnan(data.Latitude(:))))
        E.callError('nan_mismatch', 'Using match_data = true: data.Longitude and data.Latitude have different NaNs')
    elseif ~isequal(any(isnan(data.Longitude),1), data_nans)
        E.notimplemented('Using match_data = true: data.Longitude has a row that is only partially NaNs')
    end
    
    new_nans = all(isnan(lon),1);
    if any(xor(isnan(lon(:)), isnan(lat(:))))
        E.callError('nan_mismatch', 'Using match_data = true: lon and lat have different NaNs')
    end
    
    if sum(data_nans(:)) ~= sum(new_nans(:))
        E.notimplemented('Using match_data = true: the number of NaNs in the longitude in the DATA structure and the file being read are different')
    end
    
    [tmp_xx, tmp_yy] = find_submatrix2(data.Longitude(:,~data_nans), data.Latitude(:,~data_nans), lon(:,~new_nans), lat(:,~new_nans));
    if isempty(tmp_xx) || isempty(tmp_yy)
        E.callError('data_match_failure', 'Failed to find the existing data.Longitude/data.Latitude in the new files'' longitude/latitude.')
    end
    xx_sub = false(size(lon(:,~new_nans)));
    yy_sub = false(size(lon(:,~new_nans)));
    xx_sub(tmp_xx,:) = true;
    yy_sub(:,tmp_yy) = true;
    
    xx = false(size(lon));
    yy = false(size(lon));
    xx(:,~new_nans) = xx_sub;
    yy(:,~new_nans) = yy_sub;
end

%cut_alongtrack = any(xx & yy, 2);
cut_alongtrack = any(xx,2) & any(yy,2);
cut_acrosstrack = true(1,60); % keep all elements in the across track direction for now
lat = lat(cut_alongtrack, cut_acrosstrack);
lon = lon(cut_alongtrack, cut_acrosstrack); %#ok<NASGU> lon not used, will read in Latitude and Longitude like normal variables if requested. This line here just to be clear that lon, if used, should be cut down.

there_are_points = numel(lat) > 0;
if ~there_are_points
    return
end

% Row will keep track of the pixel's location in the across-track
% direction. These indices are 0 based by NASA convention.
if isfield(data,'Row')
    Row = find(cut_acrosstrack)-1;
    data.Row = repmat(Row, size(lat,1), 1);
end

for a=1:numel(sp_vars)
    dset_name = find_dset(hgrp_info, sp_vars{a});
    dset_vals = read_hdf5_dataset(sp_file, dset_name, cut_acrosstrack, cut_alongtrack, dim_order);
    %dset_vals = read_hdf5_simple(sp_file, dset_name, cut_acrosstrack, cut_alongtrack, dim_order);
    data.(sp_vars{a}) = dset_vals;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SUBFUNCTIONS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%


function dsetname = find_dset(search_group, varname)
% Recursive function that searches top-down through all groups subordinate
% the the input group to find the variable with name VARNAME
E = JLLErrors;

if ~isempty(search_group.Datasets)
    dsetname = ''; % default return value
    dsets = {search_group.Datasets.Name};
    xx = strcmp(dsets, varname);
    if any(xx)
        dsetname = strcat(search_group.Name, '/', varname);
        return
    end
else
    for a=1:numel(search_group.Groups)
        dsetname = find_dset(search_group.Groups(a), varname);
        if ~isempty(dsetname)
            return
        end
    end
    % We should only get here if we've exhausted all the groups. 
    E.callError('variable_not_found','Could not find variable %s in %s',varname,filename);
end
end


function pvec = order_to_permvec(order, xtrack_ind, ltrack_ind, n_array_dims)
% Turns a dimension order where 'x' means across track, 'l' along track,
% and 'o' any other dimensions into a vector that can be given to permute
% to reorder an array.
E = JLLErrors;

% If two dimensions in the array, then they should be along and across
% track.
if n_array_dims == 2
    order = strrep(order, 'o', '');
elseif n_array_dims > 3
    E.notimplemented('n_array_dims > 3');
end

o_ind = strfind(order, 'o');
l_ind = strfind(order, 'l');
x_ind = strfind(order, 'x');

% If one of ltrack_ind and xtrack_ind are empty, that means that we're
% probably dealing with a 1D dataset. If both are empty, that's weird.
pvec = 1:n_array_dims;
if xor(isempty(ltrack_ind), isempty(xtrack_ind))
    tmp = 1:n_array_dims;
    if isempty(ltrack_ind)
        ltrack_ind = tmp(tmp ~= xtrack_ind);
    else
        xtrack_ind = tmp(tmp ~= ltrack_ind);
    end
elseif isempty(ltrack_ind) && isempty(xtrack_ind)
    E.notimplemented('ltrack_ind and xtrack_ind are empty')
end

pvec(l_ind) = ltrack_ind;
pvec(x_ind) = xtrack_ind;
if ~isempty(o_ind)
    tmp = 1:n_array_dims;
    pvec(o_ind) = tmp( tmp ~= ltrack_ind & tmp ~= xtrack_ind );
end

end

% Benchmarking results for reading a single file: using the low level HDF5
% functions offers no appreciable performance increase over the high level
% functions in terms of run time. It does reduce memory usage by 6x in one
% test, but memory usage within both subfunctions was < 1 MB. Will use low
% level functions in case the memory reduction is valuable later (for
% parallelization perhaps)

function vals = read_hdf5_dataset(filename, dsetname, xtrack_cut, ltrack_cut, dim_order)
% Read in individual variables, apply the offset, scale factor, and fill
% value, and convert to double. This uses low-level HDF functions, I should
% benchmark this versus the simple (and more transparent in my opinion)
% approach in read_hdf5_simple

E = JLLErrors;
if ~islogical(xtrack_cut) || ~isvector(xtrack_cut)
    E.badinput('XTRACK_CUT must be a logical vector')
end
if ~islogical(ltrack_cut) || ~isvector(ltrack_cut)
    E.badinput('LTRACK_CUT must be a logical vector')
end

fileID = H5F.open(filename, 'H5F_ACC_RDONLY', 'H5P_DEFAULT');
datasetID = H5D.open(fileID, dsetname); 
dataspaceID = H5D.get_space(datasetID); 

% Figure out how many dimensions the dataset has and use that to set up the
% offset and slab size
[~, slabsize] = H5S.get_simple_extent_dims(dataspaceID);

% Figure out which dimension is along track and which across track. Note:
% these indices represent C-style ordering, which is the reverse of how the
% matrices are ordered in Matlab

xtrack_ind_c = slabsize == length(xtrack_cut);
if sum(xtrack_ind_c > 1)
    E.callError('swath_index','Multiple dimensions had the length expected for the across track dimension for %s in %s', dsetname, filename);
end
ltrack_ind_c = slabsize == length(ltrack_cut);
if sum(ltrack_ind_c > 1)
    E.callError('swath_index','Multiple dimensions had the length expected for the along track dimension for %s in %s', dsetname, filename);
end

stride = []; % an empty array tells the H5 functions to assume 1's in all dimensions
blocksize = [];
slabsize(xtrack_ind_c) = sum(xtrack_cut);
slabsize(ltrack_ind_c) = sum(ltrack_cut);
offset = zeros(size(slabsize));
offset(xtrack_ind_c) = find(xtrack_cut,1,'first')-1; % these indices are zero based, but Matlab's are 1 based
offset(ltrack_ind_c) = find(ltrack_cut,1,'first')-1;
memspaceID = H5S.create_simple(length(slabsize), slabsize, slabsize);

H5S.select_hyperslab(dataspaceID, 'H5S_SELECT_SET', offset, stride, slabsize, blocksize); 
vals = H5D.read(datasetID, 'H5ML_DEFAULT', memspaceID, dataspaceID, 'H5P_DEFAULT'); 
H5F.close(fileID);

vals=double(vals);
fillval = double(h5readatt(filename, dsetname, '_FillValue'));
scalefac = double(h5readatt(filename, dsetname, 'ScaleFactor'));
offset = double(h5readatt(filename, dsetname, 'Offset'));

fills = abs((vals - fillval)/fillval) < 1e-3;
vals(fills) = nan;
vals = (vals * scalefac) + offset;

% Matlab arrays use Fortran-style ordering, so the indices are flipped
% compared to what we had before
xtrack_ind_f = length(slabsize) + 1 - find(xtrack_ind_c);
ltrack_ind_f = length(slabsize) + 1 - find(ltrack_ind_c);
permvec = order_to_permvec(dim_order, xtrack_ind_f, ltrack_ind_f, ndims(vals));
vals = permute(vals, permvec);

end

function vals = read_hdf5_simple(filename, dsetname, xtrack_cut, ltrack_cut, dim_order)
E = JLLErrors;
if ~islogical(xtrack_cut) || ~isvector(xtrack_cut)
    E.badinput('XTRACK_CUT must be a logical vector')
end
if ~islogical(ltrack_cut) || ~isvector(ltrack_cut)
    E.badinput('LTRACK_CUT must be a logical vector')
end

vals = h5readomi(filename, dsetname);
sz = size(vals);

if isvector(vals) && numel(vals) == numel(xtrack_cut)
    vals = vals(xtrack_cut);
elseif isvector(vals) && numel(vals) == numel(ltrack_cut)
    vals = vals(ltrack_cut);
elseif ~isvector(vals)
    vals = vals(xtrack_cut, ltrack_cut, :);
    sz(1) = sum(xtrack_cut);
    sz(2) = sum(ltrack_cut);
    vals = reshape(vals, sz);
    permvec = order_to_permvec(dim_order, 1, 2, ndims(vals));
    vals = permute(vals, permvec);
else
    E.notimplemented('1D variable with size not equal to across track or along track dimension')
end
end