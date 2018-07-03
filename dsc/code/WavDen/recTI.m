function f = recTI(signal,type,h)

% recTI:     Smoothing using the Translation-Invariant procedure.
% Usage
%            f = recTI(signal,type,h)
% Inputs
%   signal	 1-d Noisy signal, length(signal)= 2^J
%   type	    'S' for soft thresholding, 'H' for hard thresholding.
%		       Optional, default= hard thresholding
%   h 		 Quadrature mirror filter for wavelet transform.
%		       Optional, default = Symmlet 8
% Outputs
%   f		    Estimate, obtained by applying thresholding on the 
%            wavelet coefficients.
% References
%            Coifman, R.R. & Donoho, D.L. (1995). Translation-invariant 
%            de-noising. In "Wavelets and Statistics", Antoniadis, A. & 
%            Oppenheim, G. (Eds.), Lect. Notes Statist., 103, pp. 125-150, 
%            New York: Springer-Verlag.

if nargin < 3,
	h = MakeONFilter('Symmlet',8);
end

if nargin < 2,
	type = 'H';
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;
thr = sqrt(2* log(n*log2(n)));

%Normalisation of the noise level to 1
[signalnorm,coef] = NormNoise(signal,h);

%Extraction of the wavelet coefficients
[tiwt] = FWT_TI(signalnorm, lev,h);
[nrow,ncol]  = size(tiwt);

if strcmp(type,'H'),
	tiwt(:,2:ncol) = HardThresh(tiwt(:,2:ncol),thr);
else
	tiwt(:,2:ncol) = SoftThresh(tiwt(:,2:ncol),thr);
end
		
f = (1/coef)*IWT_TI(tiwt,h);
      
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
       
       