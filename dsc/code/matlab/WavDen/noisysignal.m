function f = noisysignal(funfcn,n,rsnr)

% noisysignal: Generate noisy signal from a specified function.
% Usage
%              f = noisysignal(funfcn,n,rsnr)
% Inputs
%   funfcn	   String: 'Spikes','Corner','Blip','Wave' + signals made
%	            by the function MakeSignal
%   n	 	      Desired signal length
%   rsnr	      Signal to noise ratio
% Outputs
%   f		      Artificial 1-d noisy signal.


%signal
fval=MakeSignalNewb(funfcn,n);

%noise
stdnoise=std(fval)/rsnr;
epsilon=randn(size(fval))*stdnoise;

%signal + noise
fvalnoisy=fval+epsilon; 

f=fvalnoisy;

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
 