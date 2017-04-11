function [ D_new, fill_vals ] = prod_test_load_hdf( newfile, fields_to_check )
%[ D_NEW, D_OLD ] = PROD_TEST_LOAD_HDF( NEWFILE, FIELDS_TO_CHECK )
%   This function will load a BEHR HDF file and return the data in
%   structures like those used in the .mat files. NEWFILE and OLDFILE must
%   be strings pointing to the new and old files respectively.
%   FIELDS_TO_CHECK must be a cell array of fields in the HDF file to load.
%   This returns two structures D_NEW and D_OLD which will have each swath
%   as a top-level index and each field in that swath. FILL_VALS is a
%   vector of fill values read from the HDF attributes. It will have the
%   same order as FIELDS_TO_CHECK along the first dimension, the second
%   dimension will be for new and old. E.g., if you request the fields
%   AMFTrop, AMFStrat, and BEHRAMFTrop then FILL_VALS will have the form:
%
%       [ Fill AMFTrop (new),       Fill AMFTrop (old);
%         Fill AMFStrat (new),      Fill AMFStrat (old);
%         Fill BEHRAMFTrop (new),   Fill BEHRAMFTrop (old) ]


%%%%% INPUT CHECKING %%%%%
E = JLLErrors;
if ~ischar(newfile)
    E.badinput('NEWFILE must be a path given as a string')
elseif ~exist(newfile, 'file')
    E.badinput('NEWFILE (%s) does not exist.',newfile);
end

if ~iscellstr(fields_to_check)
    E.badinput('FIELDS_TO_CHECK must be a cell array of strings')
end

%%%%% MAIN FUNCTION %%%%%

hi_new = h5info(newfile);
[D_new, fill_vals] = read_hdf(hi_new, fields_to_check);

end

function [D, fills] = read_hdf(hi, fields_to_check)
% Utility function that actually does the reading.
E=JLLErrors;
nswaths = numel(hi.Groups(1).Groups);
D = make_empty_struct_from_cell(fields_to_check);
D = repmat(D, nswaths, 1);
for a=1:nswaths
    % Read in each field in each swath. If the read throws an error, see if
    % it's because the field doesn't exist. If so, rephrase the error more
    % succinctly.
    for b=1:numel(fields_to_check)
        try
            D(a).(fields_to_check{b}) = h5read(hi.Filename, h5dsetname(hi,1,a,fields_to_check{b}));
        catch err
            if strcmp(err.identifier,'MATLAB:imagesci:h5read:libraryError')
                dsets = {hi.Groups(1).Groups(a).Datasets.Name};
                if ~ismember(fields_to_check{b}, dsets)
                    E.badinput('The field %s is not present in the file %s',fields_to_check{b},hi.Filename);
                else
                    rethrow(err);
                end
            else
                rethrow(err)
            end
        end
    end
end

% Assume that the fill value is the same in every swath
fills = nan(numel(fields_to_check),1);
for a=1:numel(fields_to_check)
    dsets = {hi.Groups(1).Groups(1).Datasets.Name};
    xx = strcmp(dsets, fields_to_check{a});
    if sum(xx) < 1
        E.callError('find_fill','Could not find the field %s in %s',fields_to_check{a},hi.Filename);
    end
    fills(a) = hi.Groups(1).Groups(1).Datasets(xx).FillValue;
end
end

