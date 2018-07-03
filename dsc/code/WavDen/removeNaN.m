function y = removeNaN(x)

% removeNaN:  This function deletes NaNs from a vector.
% Usage
%             y = removeNaN(x)
% Inputs
%   x         a given vector
% Outputs
%   y         a vector without NaN's in it
% See also
%             NaN, isnan, zeroNaN

y = x(~isnan(x)); 