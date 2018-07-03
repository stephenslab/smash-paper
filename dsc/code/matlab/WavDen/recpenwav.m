function f = recpenwav(signal,h)

% recpenwav:  Smoothing using least squares penalized wavelet estimation.
%             It calls the PENWAVE procedure.
% Usage
%             f = recpenwav(signal,h)
% Inputs
%   signal	  1-d Noisy signal, length(signal)= 2^J
%   h 		  Quadrature mirror filter for wavelet transform.
%		        Optional, Default = Symmlet 8
% Outputs
%   f		     Estimate, obtained by applying linear shrinkage on the 
%             wavelet coefficients.
% References
%             1.  Antoniadis, A. (1996). Smoothing noisy data with tapered 
%                 coiflets series. Scand. J. Statist., 23, 313-330.
%             2.  Amato, U. & Vuza, D.T. (1997). Wavelet approximation of a 
%                 function from samples affected by noise. Rev. Roumanie Math. 
%                 Pure Appl., 42, 481-493.

if nargin < 2,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);

[signalnorm,coef] = NormNoise(signal,h);

aux=penwav(signalnorm,h);
f=aux./coef;

% Copyright (c) 1996
%
% Anestis Antoniadis 