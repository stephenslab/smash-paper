function f = postsinglemed(wcoef,sigma,thet0)

% postsinglemed: Estimate the wavelet coefficients by their posterior median
%                needed by the RECSINGLEMED and RECHYBBLOCKMED procedures.
% Usage
%                f = postsinglemed(wcoef,sigma,thet0)
% Inputs
%   wcoef	     Empirical wavelet coefficients at a given resolution level
%   sigma	     Level of noise
%   thet0	     Parameter vector [PIE0,TAU0] where PIE0 and TAU0 are the 
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
x=1/2+min(Omikron,1)/2;
y=tau^2/(sigma^2+tau^2)*abs(wcoef)-sigma*tau/sqrt(sigma^2+tau^2)*...
  norminv(x);
threshval=sign(wcoef).*max(0,y); 

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
 