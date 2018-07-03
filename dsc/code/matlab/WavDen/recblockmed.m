function f = recblockmed(method,signal,l,thet0,h)

% recblockmed:  Smoothing using a Bayesian block thresholding method.
%               It calls the POSTBLOCKMED procedure.
% Usage
%               f = recblockmed(method,signal,l,thet0,h)
% Inputs
%   method	    'Augment' for augmenting the data, 'Truncate' for truncating 
%                the data
%   signal	     1-d Noisy signal, length(signal)= 2^J
%   l		        Length of the blocks. Level specific block length is enabled 
%                via the choice l=[]
%   thet0	     Parameter vector [PIE0,TAU0,RH0] where PIE0, TAU0 and RH0 are 
%                the initial values for the maximization algorithm
%   h 		     Quadrature Mirror Filter for Wavelet Transform
%		           Optional, default = Symmlet 8
% Outputs
%   f		        Estimate, obtained by applying thresholding on the wavelet 
%                coefficients.
% References
%                Abramovich, F., Besbeas, P. & Sapatinas, T. (2000). Empirical Bayes 
%                approach to block wavelet function estimation. Technical Report, 
%                Department of Mathematics and Statistics, University of Cyprus, 
%                Cyprus.

if nargin < 5,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;

%perform DWT
wcoef=FWT_PO(signal,lev,h);

%extract wavelet coefficients by level
J=log2(n);
wthresh=wcoef;

%estimate sigma
highestlev=wcoef(2^(J-1)+1:2^J);
sigma=median(abs(highestlev))/0.6745;

for j=lev:J-1
 %allow for level specific block length
 if isempty(l), l2=j; else, l2=l; end
 wthresh(2^j+1:2^(j+1))=postblockmed(method,wcoef(2^j+1:2^(j+1)),l2,sigma,thet0);
end  
reconstruct=IWT_PO(wthresh,lev,h);

f=reconstruct;

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
  