function [  ] = cross_contours( hist_mat, centers )
%CROSS_CONTOURS(HIST_MAT, CENTERS) Compare interrelationships in EMG fitting parameters
%   Creates 10 contour plots depicting the relationship of the 5 fitting
%   parameters in the EMG. HIST_MAT is the array of bin counts output from
%   nd_binning.  CENTERS is the cell array of bin centers also output from
%   nd_binning.
%
%   The resulting set of plots looks at the cross relations between each
%   combination of variable pairs. The contours represent the total number
%   of points that have that particular combination of the two parameters
%   specified on the axes; summed over all values of the other three
%   parameters.
%
%   Josh Laughner <joshlaugh5@gmail.com> Feb 2016

figure;
varnames = {'a','x0','mu_x','sx','B'};

cmax = 0;
for a=1:4
    for b=a+1:5
       sp_ind = sub2ind([4,4],a,b-1);
       inds = 1:5;
       inds(inds==a | inds==b) = [];
       C = squeeze(sum(sum(sum(hist_mat,inds(1)),inds(2)),inds(3)));
       cmax = max(max(C(:)),cmax);
       [X,Y] = meshgrid(centers{a},centers{b});
       
       subplot(4,4,sp_ind);
       contourf(X,Y,C);
       colorbar
       xlabel(varnames{a});
       ylabel(varnames{b});
    end
end

% go back and normalize the color ranges
for a=1:4
    for b=a+1:5
        sp_ind = sub2ind([4,4],a,b-1);
        subplot(4,4,sp_ind);
        caxis([0 cmax]);
    end
end
end

