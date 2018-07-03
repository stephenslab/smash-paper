function f = neighblock(wcoef,L0,sigma)

% neighblock: Group the empirical wavelet coefficients into disjoint blocks 
%             of length L0, and extract the wavelet coefficient needed by 
%             the RECNEIGHBLOCK procedure.
% Usage
%             f = neighblock(wcoef,L0,sigma)
% Inputs
%   wcoef	  Empirical wavelet coefficients at a given resolution level
%   L0		  Length of the blocks
%   sigma	  Level of noise
% Outputs
%   f    	  Denoised wavelet coefficients.

   L1 = max(1,floor(L0/2));
	L = L0 + 2*L1;
	thresh = 4.505241*L*sigma^2;
	m = length(wcoef);
	
	nb = floor(m/L0); % number of blocks of length L0
	yy = [wcoef((m-L1+1):m) wcoef wcoef(1:L1)];
	SS = 1:nb;
	for b = 1:nb,
		SS(b) = sum(yy(((b-1)*L0+1):((b-1)*L0+L)).^2);
	end
	factor = max(0,1-thresh./SS);
	for b = 1:nb,
		wcoef(((b-1)*L0+1):(b*L0)) = wcoef(((b-1)*L0+1):(b*L0)) * factor(b);
	end
	
	% last block
	if m > nb*L0,
		Lb = floor((L-(m-nb*L0))/2);
		SSL = sum([wcoef((m-L+Lb+1):m) wcoef(1:Lb)].^2);
		wcoef((nb*L0+1):m) = wcoef((nb*L0+1):m)*max(0,1-thresh/SSL);
	end
	
   f = wcoef;
   
% Acknowledgement
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

	 