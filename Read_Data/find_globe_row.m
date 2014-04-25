function [ yy, rr ] = find_globe_row( row, globe_lat_matrix, globe_lon_matrix, Data )
%Find the GLOBE terrain pressure that might belong to the current OMI row
%   Josh Laughner 26 Mar 2014 <joshlaugh5@gmail.com>

rr=find(Data.Row==row);
if isempty(rr) %Do nothing if there are no entries from that row in the lat/lon boundaries
    yy = [];
    disp('Empty row passed...');
else
    x_rightedge = reshape(Data.Loncorn([2,3],rr),1,numel(Data.Loncorn([2,3],rr)));
    x_leftedge = reshape(Data.Loncorn([1,4],rr),1,numel(Data.Loncorn([1,4],rr)));
    y_rightedge = reshape(Data.Latcorn([2,3],rr),1,numel(Data.Latcorn([2,3],rr)));
    y_leftedge = reshape(Data.Latcorn([1,4],rr),1,numel(Data.Latcorn([1,4],rr)));
    
    zz = find(x_rightedge == 0); x_rightedge(zz) = []; y_rightedge(zz) = [];
    zz = find(y_rightedge == 0); x_rightedge(zz) = []; y_rightedge(zz) = [];
    zz = find(x_leftedge == 0); x_leftedge(zz) = []; y_leftedge(zz) = [];
    zz = find(y_leftedge == 0); x_leftedge(zz) = []; y_leftedge(zz) = [];
    
    %Find the end points of the left and right row edges
    left_max = find(y_leftedge==max(y_leftedge)); left_min = find(y_leftedge==min(y_leftedge));
    right_max = find(y_rightedge==max(y_rightedge)); right_min = find(y_rightedge==min(y_rightedge));
    
    %Calculate the slope for the right side as the tangent to the middle of the
    %arc.
    right_mid = floor(length(y_rightedge)/2);
    right_slope = (y_rightedge(right_mid + 1) - y_rightedge(right_mid - 1))/(x_rightedge(right_mid + 1) - x_rightedge(right_mid - 1));
    
    %Find a box that contains the row: calculate the right side points as where
    %the tangent line to the middle intersects lines crossing through the top
    %and bottom of the left side and right side
    top_slope = (y_rightedge(right_max) - y_leftedge(left_max))/(x_rightedge(right_max) - x_leftedge(left_max));
    right_top_xint = (top_slope * x_leftedge(left_max) - right_slope * x_rightedge(right_mid) - y_leftedge(left_max) + y_rightedge(right_mid))/(top_slope - right_slope);
    right_top_yint = top_slope * (right_top_xint - x_leftedge(left_max)) + y_leftedge(left_max);
    
    bottom_slope = (y_rightedge(right_min) - y_leftedge(left_min))/(x_rightedge(right_min) - x_leftedge(left_min));
    right_bottom_xint = (bottom_slope * x_leftedge(left_min) - right_slope * x_rightedge(right_mid) - y_leftedge(left_min) + y_rightedge(right_mid))/(bottom_slope - right_slope);
    right_bottom_yint = bottom_slope * (right_bottom_xint - x_leftedge(left_min)) + y_leftedge(left_min);
    
    polyx = [x_leftedge(left_min),x_leftedge(left_max),right_top_xint,right_bottom_xint,x_leftedge(left_min)];
    polyy = [y_leftedge(left_min),y_leftedge(left_max),right_top_yint,right_bottom_yint,y_leftedge(left_min)];
    
    yy = inpolygon(globe_lon_matrix,globe_lat_matrix,polyx,polyy);
end
end

