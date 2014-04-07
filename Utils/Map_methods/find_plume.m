function [ in_plume, edge_pixels ] = find_plume( matrix_in, matrix_lon, matrix_lat, threshold, center_lon, center_lat )
%Finds all pixels in the matrix_in around the center point (given by center_lon and
%center_lat) that have a value greater than the threshold value.  Returns
%the linear indices of the pixels in the plume as well as those that define
%the edge.
%   This function starts with the pixel defined by [center_lon center_lat]
%   and works outward, checking to see if each of its neighbors falls above
%   the threshold value.  Those that do are added to the plume, and their
%   neighbors in turn are checked.  This continues until the algorithm is
%   satisfied that all pixels with a concentration greater than the
%   threshold value that are continguos with the center pixel have been
%   located.
%   Josh Laughner 4 Apr 2014 <joshlaugh5@gmail.com>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Find the pixel corresponding to the center lat/lon %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take the difference between all of the lat/lon values and the center
% lat/lon, then find the smallest total difference
lat_dif = abs(matrix_lat - center_lat);
lon_dif = abs(matrix_lot - center_lon);
combined_dif = lat_dif + lon_dif;

%The linear index of the center pixel
center_pix = find(combined_dif == min(combined_dif));
if matrix_in(center_pix) < threshold;
    error('center_error:below_threshold','The center pixel value is below the stated threshold. Aborting.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initialize our starting condition and loop until the plume is located %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plume_mat = zeros(size(matrix_in));
plume_mat(center_pix) = 1;

while true
    [row, col] = find(plume_mat == 1);
    if isempty(row); break; end %If no pixels with a value of 1 are left, assume that the full plume has been found and exit.
    for a = 1:length(row)
        %Take each of the eight neighbors and test: (a) if they are already
        %part of the plume and (b) if they are higher than the threshold
        n_row = [row(a) - 1, row(a) - 1, row(a) - 1, row(a), row(a) row(a) + 1, row(a) + 1, row(a) + 1];
        n_col = [col(a) - 1, col(a), col(a) + 1, col(a) - 1, col(a) + 1, col(a) - 1, col(a), col(a) + 1];
        for b = 1:8;
            if plume_mat(n_row(b), n_col(b)) > 0 %A 0 indicates this cell has not yet been visited.  If the value is >0, skip this cell because it has already been considered
            elseif matrix_in(n_row(b), n_col(b)) >= threshold
                plume_mat(n_row(b), n_col(b)) = 1; %A value of 1 indicates that this pixel is part of the plume but its neighbors have not yet been checked.
            end                              %By setting this to 1, the loop will catch it next time and check its neighbors.
        end
        plume_mat(row(a),col(a)) = 2; %A value of 2 indicates that the neighbors of this pixel have been evaluated.
    end
end

in_plume = find(plume_mat == 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Now find the edge pixels, i.e. any with a neighbor marked 0 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[row2, col2] = find(plume_mat == 2);
edge_pixels = zeros(1,numel(row2)); edge_pixels(:)=NaN;
edge_ind = 1;
for a = 1:length(row2);
    n_row = [row2(a) - 1, row2(a) - 1, row2(a) - 1, row2(a), row2(a) row2(a) + 1, row2(a) + 1, row2(a) + 1];
    n_col = [col2(a) - 1, col2(a), col2(a) + 1, col2(a) - 1, col2(a) + 1, col2(a) - 1, col2(a), col2(a) + 1];
    
    plume_mat_sub = plume_mat(n_row,n_col);
    
    for b=1:8
        if any(plume_mat_sub(:)==0)
            lin_ind = sub2ind(size(plume_mat),row2(a),col2(a));
            edge_pixels(edge_ind) = lin_ind;
            edge_ind = edge_ind + 1;
        end
    end
end
edge_pixels(isnan(edge_pixels)) = [];

end

