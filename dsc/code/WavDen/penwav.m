function f = penwav(data,filt)

% penwav:   Least squares penalized wavelet estimation needed by the 
%           RECPENWAV procedure.
% Usage
%           f = penwav(data,filt)
% Inputs
%   data		1-d noisy signal, length(signal)= 2^J
%   filt		Quadature mirror filter for wavelet transform
% Outputs
%   f		   Vector of same length as "data" containing the estimate 
%           obtained by applying penalized least squares. 
% References	
%	         Antoniadis, A. (1996). Smoothing noisy data with tapered 
%           coiflets series. Scand. J. Statist., 23, 313-330.
% Note
%           Uses the regpengcv function.

% smoothness index
preg=3.6175; 

data_wd = FWT_PO(data, 0, filt);
data_smowd=regpengcv(data_wd,preg);
f = IWT_PO(data_smowd, 0, filt);

% Copyright (c) 1996
%
% Anestis Antoniadis
 