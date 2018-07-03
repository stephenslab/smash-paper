function f = recsinglemean3(signal,thet0,h)

% recsinglemean3: Smoothing using the POSTSINGLEMEAN3 procedure.
% Usage
%                 f = recsinglemean3(signal,thet0,h)
% Inputs
%   signal	      1-d Noisy signal, length(signal)= 2^J
%   thet0	      Parameter vector [PIE0,TAU0,SIGMA0] where PIE0, TAU0 
%                 and SIGMA0 are the initial values for the EM algorithm
%   h 		      Quadrature mirror filter for wavelet transform.
%		            Optional, default = Symmlet 8
% Outputs
%   f		         Estimate, obtained by applying thresholding on the 
%                 wavelet coefficients.
% References
%                 Clyde, M. & George, E.I. (2000). Flexible empirical Bayes 
%                 estimation for wavelets. J. R. Statist. Soc. B, 62, 681-698.

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
param = em2(thet0,wcoef,lev);
for j=lev:J-1
  wthresh(2^j+1:2^(j+1))=postsinglemean3(wcoef(2^j+1:2^(j+1)),param.pie(j), ...
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
  