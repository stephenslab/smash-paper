function f = recneighblock(signal,h)

% recneighblock: Smoothing using the NEIGHBLOCK procedure.
% Usage
%                f = recneighblock(signal,h)
% Inputs
%   signal	     1-d Noisy signal, length(signal)= 2^J
%   h 		     Quadrature mirror filter for wavelet transform.
%		           Optional, default = Symmlet 8
% Outputs
%   f 	        Estimate, obtained by applying thresholding on the
%                wavelet coefficients.
% References
%                Cai, T.T. & Silverman, B.W. (2001). Incorporating information
%                on neighboring coefficients into wavelet estimation. Sankhya,
%                Series A, 63 (to appear).

if nargin < 2,
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
  wthresh(2^j+1:2^(j+1))=neighblock(wcoef(2^j+1:2^(j+1)),floor(log(n)/2),sigma);
end
reconstruct=IWT_PO(wthresh,lev,h);

f=reconstruct;

% Acknowledgements
%
% This code is based on an Splus-code written by T.T. Cai & B.W. Silverman.
%
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

