function f = recvisu(signal,type,h)

% recvisu:   Smoothing using the VISUSHRINK procedure.
% Usage
%            f = recvisu(signal,type,h)
% Inputs
%   signal	 1-d Noisy signal, length(signal)= 2^J
%   type	    'S' for soft thresholding, 'H' for hard thresholding.
%		       Optional, default = hard thresholding
%   h 		 Quadrature mirror filter for wavelet transform
%		       Optional, default = Symmlet 8
% Outputs
%   f		    Estimate, obtained by applying thresholding on the 
%            wavelet coefficients.
% References
%            Donoho, D.L. & Johnstone, I.M. (1994). Ideal spatial 
%            adaptation by wavelet shrinkage. Biometrika, 81, 425-455.
% Note 
%            Uses the ThreshWaveb procedure from Wavelab.

if nargin < 3,
	h = MakeONFilter('Symmlet',8);
end

if nargin < 2,
	type = 'H';
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;

%Extraction of the wavelet coefficients
f = ThreshWaveb(signal,type,0,0,sqrt(2*log(n)),lev,h);

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
  