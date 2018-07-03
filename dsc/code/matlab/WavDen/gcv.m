function  gcv = gcv(x,t)

% gcv:    Applies soft-threshold t to all elements of x and computes the generalized 
%         cross-validation threshold needed by the MINGCV procedure. 
% Usage
%         gcv = gcv(x,t)
% Inputs
%   x		 input wavelet coefficients, length= 2^J
%   t     soft threshold
% Outputs
%   gcv	 generalized cross-validation threshold

  N = length(x);
  res = abs(x) - t;
  N0 = sum(res <= 0);
  y = sign(x).*res.*(res >= 0);

  if (N0 == 0)
     gcv = 0;
  else
     gcv = N * (norm(y-x) / N0 )^2 ;
  end

% Acknowledgements
%
% This code is based on a Matlab Code written and kindly provided by M. Jansen
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

 
