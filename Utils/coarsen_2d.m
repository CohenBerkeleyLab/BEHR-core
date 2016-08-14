function [ Xout, Yout, Zout ] = coarsen_2d( X, Y, Z, nx, ny, varargin )
%COARSEN_2D Coarsens 2D data
%   COARSEN_2D(X, Y, Z, NX, NY) will take data Z with X and Y coordinates
%   defined by X and Y and average NX by NY submatrices to generate coarser
%   output.
%
%   COARSEN_2D(X, Y, Z, NX, NY, 'corners') will output the first X and Y
%   value in each submatrix as the element for that coarser grid cell,
%   rather than an average. This should work well with pcolor, which plots
%   its X and Y input as the corners of the cells.

E = JLLErrors;

% Remove the if-present, is-true option(s) first
corner_bool = ismember('corners',varargin);
if corner_bool
    varargin(strcmp(varargin,'corners'))=[];
end
% Then do the regular parameter parsing
p = inputParser;
p.addParameter('op','mean')

p.parse(varargin{:});
pout = p.Results;

opname = pout.op;
allowed_ops = {'mean','nanmean','median','nanmedian','mode'};
if ~ismember(opname, allowed_ops)
    E.badinput('%s is not a valid averaging operation; allowed values for "op" are %s', opname, strjoin(allowed_ops, ', '));
end
opfxn = str2func(opname);

%%%%%%%%%%%%%%%%%%%%%%%
%%%% MAIN FUNCTION %%%%
%%%%%%%%%%%%%%%%%%%%%%%

outsz = [ceil(size(X,1)/nx), ceil(size(X,2)/nx)];
Xout = nan(outsz);
Yout = nan(outsz);
Zout = nan(outsz);

for a=1:outsz(1)
    for b=1:outsz(2)
        i = (a-1)*nx+1;
        j = (b-1)*ny+1;
        i2 = i+nx;
        j2 = j+ny;
        
        if i2 > size(X,1)
            i2 = size(X,1);
        end
        if j2 > size(X,2)
            j2 = size(X,2);
        end
        if corner_bool
            Xout(a,b) = X(i,j);
            Yout(a,b) = Y(i,j);
        else
            Xout(a,b) = opfxn(reshape(X(i:i2, j:j2), 1, []));
            Yout(a,b) = opfxn(reshape(Y(i:i2, j:j2), 1, []));
        end
        Zout(a,b) = opfxn(reshape(Z(i:i2, j:j2), 1, []));
    end
end


end

