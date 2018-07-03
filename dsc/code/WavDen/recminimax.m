function f = recminimax(signal,type,h)

% recminimax:  Smoothing using minimax thresholds. It calls 
%              the MINIMAX procedure.
% Usage
%              f = recminimax(signal,type,h)
% Inputs
%   signal	   1-d Noisy signal, length(signal)= 2^J, J = 7,8,9 or 10 
%   type	      'S' for soft thresholding, 'H' for hard thresholding.
%		         Optional, default= hard thresholding
%   h 		   Quadrature Mirror Filter for Wavelet Transform.
%		         Optional, default = Symmlet 8
% Outputs
%   f 		   Estimate, obtained by applying thresholding on the 
%              wavelet coefficients.
% References	
%              1. Donoho, D.L. & Johnstone, I.M. (1994). Ideal spatial 
%                 adaptation by wavelet shrinkage. Biometrika, 81, 425-455.
%              2. Bruce, A.G. & Gao, H.-Y. (1996). Understanding WaveShrink: 
%                 variance and bias estimation. Biometrika, 83, 727-745. 

if nargin < 3,
	h = MakeONFilter('Symmlet',8);
end

if nargin < 2,
	type = 'H';
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;
J=log2(n);

if isempty(intersect(J,[7 8 9 10])),
  fprintf('Warning: length(signal) must be 128, 256, 512 or 1024');
  f= NaN;
  return;
end

thrhard = [2.913 3.117 3.312 3.497];
thrsoft = [1.669 1.859 2.045 2.226];

% Choosing the MiniMax Threshold
if strcmp(type,'H'),
  thr = thrhard(J-6);
else
  thr = thrsoft(J-6);
end

%Normalisation of the noise level to 1
[signalnorm,coef] = NormNoise(signal,h);
wcoef = FWT_PO(signalnorm,lev,h);
wthresh = wcoef;

%Extraction of the wavelet coefficients
if strcmp(type,'H'),
	wthresh((2^lev+1):2^J)=HardThresh(wcoef((2^lev+1):2^J),thr);
else
	wthresh((2^lev+1):2^J)=SoftThresh(wcoef((2^lev+1):2^J),thr);
end
f = (1/coef)*IWT_PO(wthresh,lev,h);

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
  