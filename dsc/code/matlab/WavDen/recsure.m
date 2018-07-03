function f = recsure(signal,h)

% recsure:  Smoothing using Stein's unbiased risk estimation criterion. 
%           It calls the SURESHRINK procedure.
% Usage
%           f = recsure(signal,h)
% Inputs
%   signal	1-d noisy signal, length(signal)= 2^J
%   h 		Quadrature mirror filter for wavelet transform.
%		      Optional, default = Symmlet 8
% Outputs
%   f		   Estimate, obtained by applying thresholding on the 
%           wavelet coefficients.
% References
%           Donoho, D.L. & Johnstone, I.M. (1995). Adapting to unknown smoothness 
%           via wavelet shrinkage. J. Am. Statist. Ass., 90, 1200-1224.
% Notes 
%           Uses the NormNoise and WaveShrink procedures from Wavelab.

if nargin < 2,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;

%Normalisation of the noise level to 1
[signalnorm,coef] = NormNoise(signal,h);

%Extraction of the wavelet coefficients
[reconstruct,wcoef] = WaveShrink(signalnorm,'SURE',lev,h);

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
  