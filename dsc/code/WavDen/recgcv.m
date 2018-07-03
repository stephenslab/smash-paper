function f = recgcv(signal,h)

% recgcv:    Smoothing using a generalized cross-validation approach.
%            It calls the GCV procedure.
% Usage
%            f = recgcv(signal,h)
% Inputs
%   signal	 1-d Noisy signal, length(signal)= 2^J
%   h 		 Quadrature Mirror Filter for Wavelet Transform.
%		       Optional, default = Symmlet 8
% Outputs
%   f		    Estimate, obtained by applying thresholding on the wavelet 
%            coefficients.
%  References 	
%       1. Weyrich, N. & Warhola, G.T. (1995). De-noising using wavelets and cross validation.
%          In "Approximation Theory, Wavelets and Applications", (S.P. Sing Ed.), NATO ASI, 
%          Series C, 454, pp. 523--532.
%		  2. Jansen, M., Malfait, M. & Bultheel, A. (1997). Generalized cross validation for
%          wavelet thresholding. Signal Processing, 56, 33--44.

if nargin < 2,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;
J=log2(n);

% Normalisation of the noise level to 1
[signalnorm,coef] = NormNoise(signal,h);

% Determining the cross-validated threshold
wcoef=FWT_PO(signalnorm,lev,h);
thr=mingcv(wcoef(2^lev+1:2^J));

% Applying soft thresholding to the wavelet coefficients
wcoef(2^lev+1:2^J) = SoftThresh(wcoef(2^lev+1:2^J),thr);

reconstruct=IWT_PO(wcoef,lev,h);
f = (1/coef)*reconstruct;

% Acknowledgements
%
% This code is based on a Matlab Code written and kindly provided by M. Jansen
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
 
 
