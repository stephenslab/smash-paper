function f = threshda1(wcoef,alpha) 

% threshda1:  Estimate the wavelet coefficients using a hypothesis testing
%             method needed by the RECTHRESHAD1 procedure.
% Usage
%             f = threshda1(wcoef,alpha)
% Inputs
%   wcoef     Empirical wavelet coefficients at a given resolution level
%   alpha	  Level of the hypothesis test
% Outputs
%   f	        Denoised wavelet coefficients.

%initialisations
n=length(wcoef); 

%Setting the threshold
wcsort = sort(wcoef.^2);
k=n;
while (k > 0) & ((1-chi2cdf(wcsort(k),1)^k) < alpha),
	k = k-1;
end
if k == 0,
	thr = 0;
else
	thr = sqrt(wcsort(k));
end

f=SoftThresh(wcoef,thr);

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
  