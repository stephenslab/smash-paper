function f = rechybblockmed(method,signal,l,thet0s,thet0b,h)

% rechybblockmed:  Smoothing using a hybrid Bayesian block thresholding scheme. 
%                  The POSTSINGLEMED procedure is applied on the first few 
%                  resolution levels (which is a user-choice) and the POSTBLOCKMED 
%                  procedure is applied on the remaining resolution levels.
% Usage
%                  f = rechybblockmed(method,signal,l,thet0s,thet0b,h)
% Inputs
%   method	       'Augment' for augmenting the data, 'Truncate' for truncating 
%                  the data
%   signal	       1-d Noisy signal, length(signal)= 2^J
%   l		          Length of the blocks. Level specific block length is enabled 
%                  via the choice l=[]
%   thet0s	       Parameter vector : [PIE0s,TAU0s] where PIE0s and TAU0s are the 
%                  initial values for the EM algorithm
%   thet0b	       Parameter vector : [PIE0b,TAU0b,RH0b] where PIE0b, TAU0b and RH0b 
%                  are the initial values for the maximization algorithm
%   h 		       Quadrature Mirror Filter for Wavelet Transform.
%		             Optional, default = Symmlet 8
% Outputs
%   f	         	 Estimate, obtained by applying thresholding on the wavelet 
%                  coefficients.
% References
%                  Abramovich, F., Besbeas, P. & Sapatinas, T. (2000). Empirical 
%                  Bayes approach to block wavelet function estimation. Technical 
%                  Report, Department of Mathematics and Statistics, University 
%                  of Cyprus, Cyprus.

if nargin < 6,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;

%perform DWT
wcoef=FWT_PO(signal,lev,h);

%extract wavelet coeff's by level
J=log2(n);
J1=lev+2;
wthresh=wcoef;

%estimate sigma
highestlev=wcoef(2^(J-1)+1:2^J);
sigma=median(abs(highestlev))/0.6745;
for j=lev:J1
   wthresh(2^j+1:2^(j+1))=postsinglemed(wcoef(2^j+1:2^(j+1)),sigma,thet0s);
end
for j=(J1+1):J-1
 %allow for level specific block length
 if isempty(l), l2=j; else, l2=l; end 
 wthresh(2^j+1:2^(j+1))=postblockmed(method,wcoef(2^j+1:2^(j+1)),l2,sigma,thet0b);
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
  