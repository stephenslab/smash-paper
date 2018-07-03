function f = recsinglemed3(signal,thet0,h)

% recsinglemed3:  Smoothing using the POSTSINGLEMED3 procedure.
% Usage
%                 f = recsinglemed3(signal,thet0,h)
% Inputs
%   signal	      1-d noisy signal, length(signal)= 2^J
%   thet0	      Parameter vector [PIE0,TAU0,SIGMA0] where PIE0, TAU0 
%                 and SIGMA0 are the initial values for the EM algorithm
%   h 		      Quadrature mirror filter for wavelet transform
%		            Optional, default = Symmlet 8
% Outputs
%   f		         Estimate, obtained by applying thresholding on the 
%                 wavelet coefficients.
% References
%                 1. Abramovich, F., Sapatinas, T. & Silverman, B.W. (1998). 
%                    Wavelet thresholding via a Bayesian approach. 
%                    J. R. Statist. Soc. B, 60, 725-749.
%                 2. Johnstone, I.M. & Silverman, B.W. (1998). Empirical Bayes 
%                    approaches to mixture problems and wavelet regression. 
%                    Technical Report, School of Mathematics, University of 
%                    Bristol, UK.
% Notes	
%                 This is an empirical Bayes approach (using EM to estimate PIE, 
%                 TAU and SIGMA) based on the posterior median. It is actually a 
%                 variant of the approach of Abramovich, Sapatinas & Silverman 
%                 (1998) using estimated SIGMA from median absolute deviation and 
%                 Besov type estimates for PIE and TAU, and of Johnstone & Silverman 
%                 (1998) who just estimated SIGMA from median absolute deviation and 
%                 considered EM estimates for PIE and TAU. 

if nargin < 3,
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
param = em2(thet0,wcoef,lev);
for j=lev:J-1
  wthresh(2^j+1:2^(j+1))=postsinglemed3(wcoef(2^j+1:2^(j+1)),param.pie(j), ...
  param.tau(j),param.sigma);
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
  