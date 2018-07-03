function f = recmixed(signal,h)

% recmixed:  Smoothing method using a mixed-effects model.
%            It calls the GCVME procedure.
% Usage
%            f = recmixed(signal,h)
% Inputs
%   signal	 1-d Noisy signal, length(signal)= 2^J
%   h 		 Quadrature Mirror Filter for Wavelet Transform
%		       Optional, Default = Symmlet 8
% Outputs
%   f		    Estimate, obtained by applying Bayesian wavelet shrinkage
%            for nonparametric mixed-effects models.
% Reference
%            Huang, S.Y. & Lu, H.H.-S. (2000). Bayesian wavelet shrinkage
%            for nonparametric mixed-effects models. Statist. Sinica, 10,
%            1021-1040.
% Note
%            Use WaveLab, the function gcvme and the function fmin in
%            the optimization toolbox of Matlab. The range of the
%            smoothing parameter is [0.01,5] in fmin('gcvme',0.01,5),
%            which make it different from Table 1 in Haung & Lu (2000).

if nargin < 2,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+3;
Jmax=log(n)/log(2);

% Estimation of 1/sigma
[signalnorm,coef] = NormNoise(signal,h);

% Performing the wavelet transform
wc=FWT_PO(signal,lev,h);

Ty=wc(1:2^lev);
Sy=wc(2^lev+1:2^Jmax);
minbl=fmin('gcvme',0.01,5,[],n,coef,Sy);
betamnoise=(Sy).^2- minbl*coef^2*log(n);
betanew=(abs(betamnoise)+betamnoise)/2;
betan=betanew./(Sy);
f=IWT_PO([Ty betan],lev,h);

% Acknowledgement
%
% This code is based on a code written by S.Y. Huang & H.H.-S. Lu.
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


