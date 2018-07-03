function f = recdecompsh(signal,h)

% recdecompsh:  Smoothing method using a deterministic/stochastic 
%               approach. It calls the DECOMPSH procedure.
% Usage
%               f = recdecompsh(signal,h)
% Inputs
%   signal	    1-d Noisy signal, length(signal)= 2^J
%   h 		    Quadrature Mirror Filter for Wavelet Transform
%		          Optional, Default = Symmlet 8
% Outputs
%   f		       Estimate, obtained by applying thresholding on the 
%               wavelet coefficients.
% Reference
%               Huang, H.-C. & Cressie, N. (2000). Deterministic/Stochastic 
%               Wavelet Decomposition for Recovery of Signal from Noisy Data. 
%               Technometrics, 42, 262-276.

if nargin < 2,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;

%Extraction of the wavelet coefficients
J = log2(n);
[reconstruct,yd,ys,ns, beta] = Decompsh(signal,lev,h,J-1,0);

f = reconstruct;

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
  