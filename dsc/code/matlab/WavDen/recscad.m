function f = recscad(signal,h)

% recscad:  Smoothing using a clipped penalized wavelet estimation.
%           It calls the SCADTHRESH procedure.
% Usage
%           f = recscad(signal,h)
% Inputs
%   signal	1-d noisy signal, length(signal)= 2^J
%   h 		Quadrature mirror filter for wavelet transform.
%		      Optional, default = Symmlet 8
% Outputs
%   f		   Estimate, obtained by applying thresholding on the 
%           wavelet coefficients.
% References
%           Antoniadis, A. & Fan, J. (2001). Regularization of 
%           wavelets approximations. J. Am. Statist. Ass., 96 
%           (to appear).

if nargin < 2,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+2;
J=log2(n);

%Normalisation of the noise level to 1
[signalnorm,coef] = NormNoise(signal,h);
thr = sqrt(2*log(n) - log(1+6^2*log(n)));

%perform DWT
wcoef = FWT_PO(signalnorm,lev,h);
wthresh = wcoef;

%Extraction of the wavelet coefficients
wthresh((2^lev+1):2^J)=scadthresh(wcoef((2^lev+1):2^J),thr);

f = (1/coef)*IWT_PO(wthresh,lev,h);

% Copyright (c) 1999
%
% Anestis Antoniadis & Jianqing Fan
 