%%gridIndx2
%%arr 10/22/2008

%..........................................................................
%
% Description:
%    Given a defined grid, the procedure finds the grid indices
%    corresponding to input values of latitude, and/or longitude, and/or
%    date.
%
% Arguments:
%
%   x              scalar or vector of values at which grid cell indices are to be evaluated
%   xSAVE          vector of grid cell values
%   nGridCells     number of grid cells (require nGridCells not equal to zero)
%
%
% Purpose:
% To find indices of grid cells at a given set of values in vector x.
%
%
% This file replaces Eric's gridIndx which I was having trouble translating
% from IDL.
%..........................................................................


function [gridIndex] = gridIndx2(x, xSAVE, nGridCells)

diff=abs(x-xSAVE);

a=find(diff==min(diff));

cells=1:nGridCells;

b=cells(a);

gridIndex=b;

