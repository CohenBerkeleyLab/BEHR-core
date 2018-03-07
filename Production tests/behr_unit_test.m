function [ success ] = behr_unit_test( new, old, DEBUG_LEVEL, fid, fields_to_ignore )
%BEHR_UNIT_TEST Compare old and new BEHR data 
%   SUCCESS = BEHR_UNIT_TEST( NEW, OLD )
%   Takes two Data or OMI structures (NEW and OLD) and compares the values
%   of each field in the structures. If everything matches, SUCCESS will be
%   true, if not, it'll be false.
%
%   SUCCESS = BEHR_UNIT_TEST( NEW, OLD, DEBUG_LEVEL ) Allows you to control
%   the verbosity of the function. By default it prints out each individual
%   comparison (DEBUG_LEVEL == 2). This is useful if you need the detail,
%   but prints a lot of information. Passing 1 as DEBUG_LEVEL will only
%   print failed tests, passing 0 will turn it off completely.
%
%   SUCCESS = BEHR_UNIT_TEST( NEW, OLD, DEBUG_LEVEL, FID ) Redirects the
%   output printing from the terminal to the file with file ID FID (from
%   FOPEN).
%
%   SUCCESS = BEHR_UNIT_TEST( NEW, OLD, DEBUG_LEVEL, FID, FIELDS_TO_IGNORE )
%   FIELDS_TO_IGNORE is a cell array of strings that specifies fields that
%   should not be compared, usually because you know they will fail for a
%   good reason. For instance, the GitHead fields will almost always fail
%   because the new and old data were produced with different versions of
%   the code.

E = JLLErrors;

if ~exist('DEBUG_LEVEL', 'var')
    DEBUG_LEVEL = 2;
end

if ~exist('fid', 'var')
    % An fid of 1 will make fprint print to the command window as if no fid
    % was given
    fid = 1;
end

if ~exist('fields_to_ignore', 'var')
    fields_to_ignore = {};
elseif ~iscellstr(fields_to_ignore)
    E.badinput('FIELDS_TO_IGNORE must be a cell array of strings');
end

tol = 1e-4;

[test_field_names, new_old_fields_mapping] = compare_fields_present(new, old, DEBUG_LEVEL, fid, fields_to_ignore);
test_field_vales = compare_field_values(new, old, new_old_fields_mapping, tol, DEBUG_LEVEL, fid);

success = test_field_names && test_field_vales;
end

function [success, field_mapping] = compare_fields_present(new, old, DEBUG_LEVEL, fid, ignore_fields)
new_fields = fieldnames(new);
old_fields = fieldnames(old);

if length(ignore_fields) > 0
    rr = regmember(new_fields, ignore_fields);
    new_fields(rr) = [];
    rr = regmember(old_fields, ignore_fields);
    old_fields(rr) = [];

    fprintf(fid, '  Ignoring fields:\n\t%s\n', strjoin(ignore_fields, '\n\t'));
end

in_new_only = ~ismember(new_fields, old_fields);
new_only_fields = new_fields(in_new_only);
in_old_only = ~ismember(old_fields, new_fields);
old_only_fields = old_fields(in_old_only);

% Allow for fields to have changed capitalization
field_mapping = {};
b = 1;
for a=1:numel(new_only_fields)
    xx = strcmpi(new_only_fields{a}, old_only_fields);
    if sum(xx) == 1
        field_mapping(b,:) = [new_only_fields(a), old_only_fields(xx)];
        in_new_only(strcmp(new_only_fields{a}, new_fields)) = false;
        in_old_only(strcmp(old_only_fields{xx}, old_fields)) = false;
        b = b+1;
    elseif sum(xx) > 1
        fprintf(fid, '  The field %s in the new structure has multiple case-insensitive matches in the old file (%s)\n', new_only_fields{a}, strjoin(old_only_fields(xx), ', '));
    end
end

if DEBUG_LEVEL > 0
    if sum(in_new_only) > 0
        fprintf(fid, '  The following fields are only present in the new structure:\n');
        fprintf(fid, '    %s\n', strjoin(new_fields(in_new_only),'\n    '));
    end
    if sum(in_old_only) > 0
        fprintf(fid, '  The following fields are only present in the old structure:\n');
        fprintf(fid, '    %s\n', strjoin(old_fields(in_old_only),'\n    '));
    end
end

% Add the fields that exactly match to the mapping
xx = ismember(new_fields, old_fields);
xx = xx(:);
common_fields_mapping = repmat(new_fields(xx), 1, 2);

field_mapping = cat(1, field_mapping, common_fields_mapping);

success = sum(in_new_only) == 0 && sum(in_old_only) == 0;
end

function success = compare_field_values(new, old, mapping, tolerance, DEBUG_LEVEL, fid)
% new - the new Data or OMI structure
% old - the old Data or OMI structure
% mapping - an n-by-2 cell array with the new fields to check in the first
%   column and the corresponding old field in the second column.
% tolerance - the maximum absolute difference allowed between two values
%   for them to be considered equal.
% DEBUG_LEVEL - how verbose to be.

success = true;
if numel(new) ~= numel(old)
    if DEBUG_LEVEL > 0
        fprintf(fid, '  Number of swaths in new (%d) and old(%d) unequal', numel(new), numel(old));
    end
    success = false;
    return
end

for a=1:numel(new)
    if DEBUG_LEVEL > 0
        fprintf(fid, '\n  Checking swath %d\n', a);
    end
    for b=1:size(mapping,1)
        if DEBUG_LEVEL > 1
            fprintf(fid, '    Checking field %s -> %s: ', mapping{b,1}, mapping{b,2});
        end
        [eq, reason] = test_field_equality(new(a).(mapping{b,1}), old(a).(mapping{b,2}), tolerance);
        success = success && eq > 0;
        
        if DEBUG_LEVEL == 1 && eq <= 0
            % If DEBUG_LEVEL > 1, we've already printed out the field
            % name. If == 1, we only print when there's a problem
            fprintf(fid, '    Field %s/%s: ', mapping{b,1}, mapping{b,2});
        end
        if DEBUG_LEVEL > 0
            if eq > 0
                if DEBUG_LEVEL > 1
                    fprintf(fid, 'PASS\n');
                end
            else
                fprintf(fid, 'FAILED (%s)\n', reason);
            end
        end
        
    end
end

end

function [eq, reason] = test_field_equality(new_val, old_val, tolerance_scale)
reason = 'reason unspecified';
% Scale up the tolerance based on the smallest magnitude of the values
% present, if needed. This will prevent false positives where, e.g. the
% BEHR VCDs are different because of a floating point difference in the AMF
% that gets scaled up. In the old gridding code, the flag fields are stored
% as cell arrays in the grid; this causes this calculation to error (and
% isn't important because the flags are around 0-256 in value) so skip it.
% One other field is a structure, so that needs skipped too.
if isnumeric(new_val) && isnumeric(old_val)
    if isinteger(new_val) && isinteger(old_val)
        eq = isequal(new_val, old_val);
        if ~eq
            reason = 'comparison of field values with "isequal" returned false';
        end
        return
    elseif xor(isinteger(new_val), isinteger(old_val))
        reason = 'fields are not the same type; one is an integer, one is a float';
        eq = false;
        return
    else
        tolerance = round(tolerance_scale * min(abs([new_val(:); old_val(:)])), 2, 'significant');
        tolerance = max(tolerance, tolerance_scale);
        
        if ~isequal(size(new_val), size(old_val))
            eq = 0;
            reason = sprintf('size of numeric field differs - %s vs %s', mat2str(size(new_val)), mat2str(size(old_val)));
            return
        end
        newnans = isnan(new_val(:));
        oldnans = isnan(old_val(:));
        if any(xor(newnans, oldnans))
            eq = 0;
            reason = 'NaNs are different';
            return
        end
        del = new_val(~newnans) - old_val(~oldnans);
        eq = all(abs(del) <= tolerance);
        if ~eq
            reason = sprintf('at least one absolute difference exceeds tolerance of %g; min/max diff = %g/%g', tolerance, min(del), max(del));
        end
    end
else
    eq = isequal(new_val, old_val);
    if ~eq
        reason = 'comparison of field values with "isequal" returned false';
    end
end
end
