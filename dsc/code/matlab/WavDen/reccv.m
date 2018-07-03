function f = reccv(signal,type, h)

% recgcv:    Smoothing using a generalized cross validation approach.
%            It calls the CV procedure.
% Usage
%            f = reccv(signal,type,h)
% Inputs
%   signal	 1-d Noisy signal, length(signal)= 2^J
%   type	    'S' for soft thresholding, 'H' for hard thresholding
%   h 		 Quadrature Mirror Filter for Wavelet Transform.
%		       Optional, default = Symmlet 8
% Outputs
%   f		    Estimate, obtained by applying thresholding on the wavelet 
%            coefficients.
% References 	
%            Nason, G.P. (1996). Wavelet shrinkage using cross-validation.
%            J. R. Statist. Soc. B, 58, 463-479.

if nargin < 3,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;
J=log2(n);

%Normalisation of the noise level to 1
[signalnorm,coef] = NormNoise(signal,h);

%Determining the cross-validated threshold
thr=fmins('cv',sqrt(2*log(n)),[],[],signalnorm,lev,h,type);
thr=thr/sqrt((1-log(2)/log(n)));

%Extraction of the wavelet coefficients
wcoef=FWT_PO(signalnorm,lev,h);
if strcmp(type,'H'),
  wcoef(2^lev+1:2^J) = HardThresh(wcoef(2^lev+1:2^J),thr);
else
  wcoef(2^lev+1:2^J) = SoftThresh(wcoef(2^lev+1:2^J),thr);
end

reconstruct=IWT_PO(wcoef,lev,h);
f = (1/coef)*reconstruct;

% Acknowledgements
%
% This code is based on an Splus-code written by G.P. Nason
%
% Copyright (c) 2001
%
% Anestis Antoniadis, Jeremie Bigot
% Laboratoire IMAG-LMC
% University Joseph Fourier
% BP 53, 38041 Grenoble Cedex 9
% France.
%
% mailto: Anestis.Antoniadis@imag.fr
% mailto: Jeremie.Bigot@imag.fr
%
% and
%
% Theofanis Sapatinas
% Department of Mathematics and Statistics
% University of Cyprus
% P.O. Box 20537  
% CY 1678 Nicosia
% Cyprus.
%
% mailto: T.Sapatinas@ucy.ac.cy 
 
 
