function f = recblockJS(method,signal,h)

% recblockJS:  Smoothing using a blockwise James-Stein rule.
%              It calls the BlockJS procedure.
% Usage
%              f = recblockJS(method,signal,h)
% Inputs
%   method	   'Augment' for augmenting the data, 'Truncate' for truncating 
%              the data
%   signal	   1-d Noisy signal, length(signal)= 2^J
%   h 		   Quadrature Mirror Filter for Wavelet Transform
%		         Optional, Default = Symmlet 8
% Outputs
%   f		      Estimate, obtained by applying thresholding on the wavelet 
%              coefficients.
% References
%              Cai, T.T. (1999). Adaptive wavelet estimation: a block 
%              thresholding and oracle inequality approach. Ann. Statist., 
%              27, 898-924.

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

%length of the blocks
L = floor(log(n));
for j=lev:J-1
 wthresh(2^j+1:2^(j+1))=blockJS(method,wcoef(2^j+1:2^(j+1)),L,sigma);
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
  
