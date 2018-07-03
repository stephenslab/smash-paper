function f = recbams(signal,h)

% recbams:  Smoothing using a Bayesian shrinkage method. It calls
%           the BAMSMOOTH procedure.
% Usage
%           f = recbams(signal,h)
% Inputs
%   signal	1-d Noisy signal, length(signal)= 2^J
%   h 		Quadrature Mirror Filter for Wavelet Transform
%		      Optional, Default = Symmlet 8
% Outputs
%   f		   Estimate, obtained by applying a Bayesian wavelet shrinkage
%           on the wavelet coefficients.
% References	
%           Vidakovic, B. & Ruggeri, F. (2000). BAMS method: theory and 
%           simulations. Discussion Paper, Institute of Statistics and 
%           Decision Sciences, Duke University, USA.

%perform DWT
if nargin < 2,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;

%Extraction of the wavelet coefficients
J = log2(n);
reconstruct = bamsmooth(signal,lev,h,J-1);
f = reconstruct;

% Copyright (c) 2000  
%
% B. Vidakovic & F. Ruggeri

 