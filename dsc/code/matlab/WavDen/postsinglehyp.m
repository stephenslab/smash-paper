function f = postsinglehyp(wcoef,sigma,thet0)

% postsinglehyp: Estimate the wavelet coefficients using a Bayesian 
%                hypothesis testing approach, needed by the RECSINGLEHYP 
%                procedure.
% Usage
%                f = postsinglehyp(wcoef,sigma,thet0)
% Inputs
%   wcoef	     Empirical wavelet coefficients at a given resolution level
%   sigma	     Level of noise
%   thet0	     Parameter vector : [PIE0,TAU0] where PIE0 and TAU0 are the 
%                initial values for the EM algorithm
% Outputs
%   f		        Denoised wavelet coefficients.

%initialisations
n=length(wcoef); 

%parameter estimation
thet=em_s(thet0,sigma,wcoef);
pie=thet(1); tau=thet(2);

e=exp(-tau^2*wcoef.^2/(2*sigma^2*(sigma^2+tau^2)));
Omikron=(1-pie)/pie*sqrt(sigma^2+tau^2)/sigma*e;
eta=1./(1+Omikron);
threshval=(Omikron<1).*wcoef;

f=threshval;

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
 