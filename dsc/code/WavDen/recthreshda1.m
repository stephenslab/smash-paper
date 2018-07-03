function f = recthreshda1(signal,h)

% recthreshda1:  Smoothing using a hypothesis testing method.
%                It calls the THRESHDA1 procedure.
% Usage
%                f = recthreshda1(signal,h)
% Inputs
%   signal	     1-d noisy signal, length(signal)= 2^J
%   h 		     Quadrature mirror filter for wavelet transform.
%		           Optional, default = Symmlet 8
% Outputs
%   f		        Estimate, obtained by applying thresholding on the 
%                wavelet coefficients.
% References
%                1. Ogden, R.T. & Parzen, E. (1996a). Change-point approach 
%                   to data analytic wavelet thresholding. Statist. Comput., 
%                   6, 93-99.
%                2. Ogden, R.T. & Parzen, E. (1996b). Data dependent wavelet 
%                   thresholding in nonparametric regression with change-point 
%                   applications. Comput. Statist. Data Anal., 22, 53-70.

if nargin < 2,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;

%Normalisation of the noise level to 1
[signalnorm,coef] = NormNoise(signal,h);
wcoef=FWT_PO(signalnorm,lev,h);

%Extraction of the wavelet coefficients
J=log2(n);
wthresh=wcoef;

for j=lev:J-1
 wthresh(2^j+1:2^(j+1))=threshda1(wcoef(2^j+1:2^(j+1)),0.05);
end
reconstruct=IWT_PO(wthresh,lev,h);

f = (1/coef)*reconstruct;

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
 
 