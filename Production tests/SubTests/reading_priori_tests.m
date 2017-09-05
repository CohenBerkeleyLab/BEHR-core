function [ success ] = reading_priori_tests( data, DEBUG_LEVEL, fid )
%READING_PRIORI_TESTS Sanity check of data imported from NASA SP 2, etc.
%   Detailed explanation goes here

if ~exist('DEBUG_LEVEL', 'var')
    DEBUG_LEVEL = 2;
end

if ~exist('fid', 'var')
    % An fid of 1 will make fprint print to the command window as if no fid
    % was given
    fid = 1;
end

success = true;
success = check_pixel_corners(data) && success;

    function passfail = check_pixel_corners(data)
        % Verify that:
        %   1) If a pixel center is not NaN, all it's corners are not NaN
        %   as well
        %   2) That the pixel corners are not arranged in such a way that a
        %   border crosses (i.e. are clockwise or counter-clockwise)
        %   3) That a pixel center falls within the border defined by the
        %   corners
        
        if DEBUG_LEVEL > 0
            fprintf(fid, 'Checking validity of corner fields...\n');
        end
        
        passfail_swaths = true(size(data));
        corner_fields = {'FoV75CornerLongitude', 'FoV75CornerLatitude';...
                         'TiledCornerLongitude', 'TiledCornerLatitude'};
                     
        for a=1:numel(data)
            this_data = data(a);
            
            for b=1:size(corner_fields,1)
                lon = this_data.Longitude;
                lat = this_data.Latitude;
                loncorn = this_data.(corner_fields{b,1});
                latcorn = this_data.(corner_fields{b,2});
                
                failed_nans = false(size(lon));
                failed_crossed = false(size(lon));
                failed_contained = false(size(lon));
                failed_other = false(size(lon));
                
                zero_flag = false;
                
                for c=1:numel(lon);
                    x = lon(c);
                    y = lat(c);
                    xall = loncorn(:,c);
                    yall = latcorn(:,c);
                    
                    % Check NaNs
                    if (any(isnan([x;y])) && ~all(isnan([x;y]))) || ...
                            (any(isnan([x;y])) && ~all(isnan([xall; yall])))
                    	% NaNs are considered inconsistent under two cases:
                        %   (i) If one but not both of the center
                        %   coordinates are NaNs
                        %   (ii) If the center coordinates are NaNs but the
                        %   corners are not. The reverse is not considered
                        %   inconsistent (as of 10 May 2017) b/c there are
                        %   days when the corner product (for whatever
                        %   reason) is missing values but the main product
                        %   is not.
                        failed_nans(c) = true;
                        continue
                    elseif all(isnan([xall;yall]))
                        % If all the corners are NaNs, then we cannot (and
                        % do not need to) do any further tests.
                        continue
                    end  
                    
                    % Found some orbits where multiple pixel corners are 0
                    % for certain rows. For now, assume that these
                    % correspond to some spacecraft maneuver (e.g. orbits
                    % 29479-29482, http://projects.knmi.nl/omi/research/calibration/instrument_status_v3/events.html)
                    % As long as it coincides with a NaN in the total NO2
                    % column, we will not consider this a failure since NO2
                    % was not retrieved for that pixel.
                    if sum(xall == 0) > 1 && sum(yall == 0) > 1
                        zero_flag = true;
                        if ~isnan(this_data.ColumnAmountNO2)
                            failed_other(c) = true;
                        end
                        continue
                    end
                    
                    % Check crossed
                    [x_chk, y_chk] = uncross_pix_corners(xall, yall);
                    if ~isequal(xall, x_chk) || ~isequal(yall, y_chk)
                        failed_crossed(c) = true;
                        continue
                    end
                    
                    % Check contained
                    if ~inpolygon(x,y,xall,yall)
                        failed_contained(c) = true;
                        continue % this continue not strictly necessary, but would be if expanded the tests
                    end
                end
                % Print out the results
                field_stem = regexprep(corner_fields{b,1}, 'Longitude|Latitude', '');
                this_field_passfail = ~(any(failed_nans(:)) || any(failed_crossed(:)) || any(failed_contained(:)));
                if DEBUG_LEVEL > 0
                    fprintf(fid, '  Swath %d, %s fields: %s\n', a, field_stem, passfail_str(this_field_passfail));
                    % Print special messages before the swath
                    if zero_flag
                        fprintf(fid, ['    NOTE: at least one pixel in this swath had multiple corners set to 0.\n'...
                                      '      I''ve only seen this happen during special maneuvers of the spacecraft.\n'...
                                      '      As long as no pixels fail "for other reasons", that means that the 0 corner\n'...
                                      '      pixels have fill values for total NO2 column and so can probably be ignored.\n']);
                    end
                end
                if ~this_field_passfail
                    passfail_swaths(a) = false;
                    if DEBUG_LEVEL > 1
                        fprintf(fid, '    %d pixels had inconsistent nans.\n', sum(failed_nans(:)));
                        if DEBUG_LEVEL > 2
                            fprintf(fid, '    Indicies are:\n');
                            fprintf(fid, index_str(failed_nans, '      '));
                            fprintf(fid, '\n');
                        end
                        
                        fprintf(fid, '    %d pixels had crossed corners.\n', sum(failed_crossed(:)));
                        if DEBUG_LEVEL > 2
                            fprintf(fid, '    Indices are:\n');
                            fprintf(fid, index_str(failed_crossed, '      '));
                            fprintf(fid, '\n');
                        end
                        
                        fprintf(fid, '    %d pixels had corners and centers misaligned.\n', sum(failed_contained(:)));
                        if DEBUG_LEVEL > 2
                            fprintf(fid, '    Indices are:\n');
                            fprintf(fid, index_str(failed_contained, '      '));
                            fprintf(fid, '\n');
                        end
                        
                        fprintf(fid, '    %d pixels failed for other reasons.\n', sum(failed_other(:)));
                        if DEBUG_LEVEL > 2
                            fprintf(fid, '    Indices are:\n');
                            fprintf(fid, index_str(failed_other, '      '));
                            fprintf(fid, '\n');
                        end
                    end
                end
            end
        end
        
        passfail = all(passfail_swaths);
    end

end

function s = passfail_str(b)
if b
    s = 'PASS';
else
    s = 'FAIL';
end
end

function s = index_str(failed_bool, beginning_space)
[xi,yi] = find(failed_bool);
ci = cell(1,numel(xi));
for a=1:numel(xi)
    ci{a} = mat2str([xi(a), yi(a)]);
end
delim = ['\n', beginning_space];
s = [beginning_space, strjoin(ci, delim), '\n'];
end
