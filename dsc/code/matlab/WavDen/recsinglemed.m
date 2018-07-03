function f = recsinglemed(signal,thet0,h)

% recsinglemed:  Smoothing using the POSTSINGLEMED procedure.
% Usage
%                f = recsinglemed(signal,thet0,h)
% Inputs
%   signal	     1-d noisy signal, length(signal)= 2^J
%   thet0	     Parameter vector : [PIE0,TAU0] where PIE0 and TAU0 are 
%                the initial values for the EM algorithm
%   h 		     Quadrature mirror filter for wavelet transform
%		           Optional, default = Symmlet 8
% Outputs
%   f		        Estimate, obtained by applying thresholding on the the 
%                wavelet coefficients.
% References	
%                Johnstone, I.M. & Silverman, B.W. (1998). Empirical Bayes 
%                approaches to mixture problems and wavelet regression. 
%                Technical Report, School of Mathematics, University of Bristol, 
%                UK.

if nargin < 3,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;

%perform DWT
wcoef=FWT_PO(signal,lev,h);

%extract wavelet coeff's by level
J=log2(n);
wthresh=wcoef;

%estimate sigma
highestlev=wcoef(2^(J-1)+1:2^J);
sigma=median(abs(highestlev))/0.6745;
for j=lev:J-1
 wthresh(2^j+1:2^(j+1))=postsinglemed(wcoef(2^j+1:2^(j+1)),sigma,thet0);
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
  