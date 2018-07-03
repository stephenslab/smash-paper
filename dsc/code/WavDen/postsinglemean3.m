function f = postsinglemean3(wcoeflev,pie,tau,sigma)

% postsinglemean3: Estimate the wavelet coefficients by their posterior
%                  mean needed by the RECSINGLEMEAN3 procedure. 
% Usage
%                  f = postsinglemean3(wcoeflev,pie,tau,sigma)
% Inputs
%   wcoeflev	    Empirical wavelet coefficients at a given resolution level
%   pie		       Proportion of nonzero wavelet coefficients at a given 
%                  resolution level
%   tau		       Variance of the wavelet coefficients at a given resolution level
%   sigma	       Level of noise
% Outputs
%   f		          Denoised wavelet coefficients.

%initialisations
n=length(wcoeflev); 

e=exp(-tau^2*wcoeflev.^2/(2*sigma^2*(sigma^2+tau^2)));
Omikron=(1-pie)/pie*sqrt(sigma^2+tau^2)/sigma*e;
eta=1./(1+Omikron);
threshval=eta.*tau^2/(sigma^2+tau^2).*wcoeflev;

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
 