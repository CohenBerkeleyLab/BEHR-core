function [ fill_val ] = behr_fill_val(  )
% BEHR_FILL_VAL - returns a standard fill value to be used in BEHR calculations.

% Chosen because it is the most negative single precision floating point
% number available in matlab (-realmax('single'))
fill_val = -3.402e38;

end

