function [ dataset, fill_value, scale_factor, offset ] = h5readomi( filename, datasetname, varargin )
%H5READOMI Read a dataset in, removing fill values and applying scale and offset
%   DATASET = H5READOMI( FILENAME, DATASETNAME ) returns the dataset in
%   FILENAME specified by DATASETNAME with fill values converted into NaNs,
%   and the scale factor and offset applied as:
%       dataset = (dataset - offset) * scale_factor
%
%   DATASET = H5READOMI( FILENAME, DATASETNAME, 'fill_crit', FILL_CRIT )
%   allows you to specify the relative difference allowed between a value
%   in the data set and the fill value for the value to be considered a
%   fill value. Defaults to 0.001, i.e. (val - fill_val)/val < fill_crit
%   must be met for VAL to be considered a fill value.
%
%   [ DATASET, FILL_VALUE, SCALE_FACTOR, OFFSET ] = H5READOMI( __ ) returns
%   the other attributes as well as the dataset.

E = JLLErrors;
p = inputParser;
p.addParameter('fill_crit',0.001);
p.parse(varargin{:});
pout = p.Results;

fill_crit = pout.fill_crit;
if ~isnumeric(fill_crit) || ~isscalar(fill_crit) || fill_crit <= 0
    E.badinput('FILL_CRIT must be a positive, scalar number')
end

dataset = h5read(filename, datasetname);
fill_value = h5readatt(filename, datasetname, '_FillValue');
scale_factor = h5readatt(filename, datasetname, 'ScaleFactor');
offset = h5readatt(filename, datasetname, 'Offset');

fills = abs(dataset - fill_value) ./ dataset < fill_crit;
dataset(fills) = nan;
dataset = (double(dataset) - offset)*scale_factor;

end

