function [ diff_table, diff_stats ] = compare_data_fields( OldData, NewData, field_name )
%compare_data_fields Compare values of pixels in BEHR Data structure
%   A testing function that compares all swaths of two BEHR Data
%   structures. Outputs a table with the swath #, pixel indices, old value,
%   new value, difference, and percent difference. If the new pixel does
%   not exist in the old swath, a value of -Inf is given to it.
%
%   Josh Laughner <joshlaugh5@gmail.com> 26 May 2015

E = JLLErrors;
DEBUG_LEVEL = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isstruct(OldData) || ~isstruct(NewData)
    E.badinput('The first two arguments must be Data structures')
elseif ~isfield(OldData,'Latitude') || ~isfield(OldData, 'Longitude') || ~isfield(NewData,'Latitude') || ~isfield(NewData, 'Longitude')
    E.badinput('Both input structures must have latitude and longitude data')
elseif ~ischar(field_name)
    E.badinput('The field name must be input as a string')
elseif ~isfield(OldData,field_name) || ~isfield(NewData,field_name)
    E.badinput('The field %s is not present in one of the input structures',field_name)
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

diff_struct = make_empty_struct_from_cell({'Swath','OldIndices','NewIndices','OldValue','NewValue','AbsDifference','PercentDifference'});

% Loop through each swath, line up pixels by their lat/lon, and compare the
% values

new_swaths = nan(1,numel(NewData));
for a=1:numel(NewData)
    new_swaths(a) = mode(NewData(a).Swath(NewData(a).Swath>0));
end

for a=1:numel(OldData)
    old_swath = mode(OldData(a).Swath(OldData(a).Swath>0));
    if isnan(old_swath)
        E.callError('null_swath','Swath value is a NaN')
    elseif ~ismember(old_swath, new_swaths)
        warning('Swath %d is not present in NewData',old_swath)
        continue
    end
    
    b = find(new_swaths == old_swath);
    
    sz = size(OldData(a).Latitude);
    for p=1:numel(OldData(a).Latitude)
        lat = OldData(a).Latitude(p);
        lon = OldData(a).Longitude(p);
        old_val = OldData(a).(field_name)(p);
        [old_ind(1), old_ind(2)] = ind2sub(sz,p);
        
        yy = lat == NewData(b).Latitude;
        xx = lon == NewData(b).Longitude;
        
        pp = xx(:) & yy(:);
        
        if sum(pp) == 0
            new_val = -Inf;
            new_ind = nan(1,2);
        elseif sum(pp) == 1
            new_val = NewData(b).(field_name)(pp);
            [new_ind(1), new_ind(2)] = ind2sub(size(NewData(b).(field_name)), find(pp));
        else
            E.callError('multiple_pixels','Too many pixels found that match');
        end
        
        abs_diff = new_val - old_val;
        rel_diff = abs_diff / old_val * 100;
        
        diff_struct.Swath = cat(1,diff_struct.Swath,old_swath);
        diff_struct.OldIndices = cat(1,diff_struct.OldIndices,old_ind);
        diff_struct.NewIndices = cat(1,diff_struct.NewIndices,new_ind);
        diff_struct.OldValue = cat(1,diff_struct.OldValue,old_val);
        diff_struct.NewValue = cat(1,diff_struct.NewValue,new_val);
        diff_struct.AbsDifference = cat(1,diff_struct.AbsDifference,abs_diff);
        diff_struct.PercentDifference = cat(1,diff_struct.PercentDifference,rel_diff);
        
    end
end

diff_table = struct2table(diff_struct);

if nargout > 1
    diff_stats.NumberOfPixels = numel(diff_struct.OldValue);
    diff_stats.NumberMismatches = sum(diff_struct.PercentDifference > 0.0001);
    diff_stats.NumberMissingInNewData = sum(isinf(diff_struct.NewValue));
    
    mad = ~isnan(diff_struct.AbsDifference) & ~isinf(diff_struct.AbsDifference);
    
    diff_stats.MeanAbsDifference = nanmean(diff_struct.AbsDifference(mad));

    mpd = ~isnan(diff_struct.PercentDifference) & ~isinf(diff_struct.PercentDifference);
    
    diff_stats.MeanPercentDifference = nanmean(diff_struct.PercentDifference(mpd));
    
    xor_mat = false(numel(diff_struct.OldValue),1);
    for a=1:numel(diff_struct.OldValue)
        xor_mat(a) = xor(isnan(diff_struct.OldValue(a)),isnan(diff_struct.NewValue(a)));
    end
    diff_stats.NaNMismatches = sum(xor_mat & ~isinf(diff_struct.NewValue));
    
end

