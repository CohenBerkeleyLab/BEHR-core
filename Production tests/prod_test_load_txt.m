function [ D_new, fillvals ] = prod_test_load_txt( newfile, fields_to_check )
%[ D_NEW, D_OLD, FILLVALS ] = PROD_TEST_LOAD_TXT( NEWFILE, FIELDS_TO_CHECK )
%   This function will load a BEHR .txt file and return the data in
%   structures like those used in the .mat files. However, because of how
%   the .txt files are formatted, each field will just contain a vector
%   rather than a matrix. NEWFILE and OLDFILE must be strings pointing to
%   the new and old files respectively. This returns two structures D_NEW
%   and D_OLD which will have each swath as a top-level index and each
%   field in that swath. FILL_VALS is a vector of fill values read from the
%   HDF attributes. It will have the same order as the fields in DATA along
%   the first dimension, the second dimension will be for new and old.
%   E.g., if Data had only the fields AMFTrop, AMFStrat, and BEHRAMFTrop
%   then FILL_VALS will have the form:
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


[D_new_tmp, fillvals_new] = read_file(newfile);

% Finally, cut down the structures to just the requested fields; this
% mimics the behavior of prod_test_load_hdf.m most closely. 
fillvals = nan(numel(fields_to_check), 1);
fns_new = fieldnames(D_new_tmp);
xx = ismember(fns_new, fields_to_check);
if sum(xx) < numel(fields_to_check)
    yy = ~ismember(fields_to_check, fns_new);
    E.callError('unknown_field','Fields %s are not present in the files.',strjoin(fields_to_check(yy), ', '));
else
    fillvals = fillvals_new(xx);
    for a=1:numel(fields_to_check)
        D_new.(fields_to_check{a}) = D_new_tmp.(fields_to_check{a});
    end
end



end

function [Data, fillvals] = read_file(filename)
% Open the file and get the variable names first, then read in the rest.
fid = fopen(filename);
cleanobj = onCleanup(@() fclose(fid));
possible_fills = [-32.77, -32770, -9e9, -1.268e30, -3.402e38];

line = fgetl(fid);
varnames = strsplit(line,',');
format_spec = repmat('%f,',1,numel(varnames)-1);
format_spec = [format_spec, '%f'];
M = fscanf(fid,format_spec); % All values should be floats. Since we've already read in the first line, it'll start on the second.
M = reshape(M,numel(varnames),[])'; % It imports as one long vector, get it back into the shape it is in the file.

Data = make_empty_struct_from_cell(varnames);
fillvals = nan(size(varnames))';
for a=1:numel(varnames);
    Data.(varnames{a}) = M(:,a);
    % Guess the fill value from the possibilities. First try direct
    % equality
    foundit = false;
    for b=1:numel(possible_fills)
        if any(M(:,a) == possible_fills(b))
            fillvals(a) = possible_fills(b);
            foundit = true;
            break
        end
    end
    % Then allow some floating point error
    if ~foundit
        for b=1:numel(possible_fills)
            if any(abs(M(:,a)-possible_fills(b)) < 0.1)
                fillvals(a) = possible_fills(b);
                break
            end
        end
    end
end
end