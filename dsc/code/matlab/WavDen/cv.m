function M = cv(thr,x,L,h,type)

% gcv:    Generalized cross-validation approach needed by the RECCV procedure.
% Usage
%         M = cv(thr,x,L,h,type)
% Inputs
%   thr   Threshold
%   x		 1-d Noisy signal, length(signal)= 2^J
%   L		 Low resolution level
%   h		 Quadrature Miror Filter
%   type	 'S' for soft thresholding, 'H' for hard thresholding
% Outputs
%   M 	 Mean squared error.
% References 	
%            Nason, G.P. (1996). Wavelet shrinkage using cross-validation.
%            J. R. Statist. Soc. B, 58, 463-479.


if (L < 0)
     fprintf('Warning: primary level must be => 0\n');
     M = NaN;
     return;
end

% Initialisations
n=length(x);
n1=n/2;
J=log2(n);

% Selecting half of the data points
xodd=[];
xeven=[];
for k=1:n1,
   xodd=[xodd; x(2*k-1)];
   xeven=[xeven; x(2*k)];
end;

% Performing WT
wcodd=FWT_PO(xodd,L,h);
wceven=FWT_PO(xeven,L,h);

if strcmp(type,'H'),
  wcodd(2^L+1:n1) = HardThresh(wcodd(2^L+1:n1),thr);
  wceven(2^L+1:n1) = HardThresh(wceven(2^L+1:n1),thr);
else
  wcodd(2^L+1:n1) = SoftThresh(wcodd(2^L+1:n1),thr);
  wceven(2^L+1:n1) = SoftThresh(wceven(2^L+1:n1),thr);
end
hatxodd = IWT_PO(wcodd,L,h);
hatxeven = IWT_PO(wceven,L,h);

% Interpolating
barxodd = 0.5.*(hatxodd(1:n1-1) + hatxodd(2:n1));
barxodd(n1) = hatxodd(1);
barxeven = 0.5.*(hatxeven(1:n1-1) + hatxeven(2:n1));
barxeven(n1) = hatxeven(1);

M=norm(xodd-barxeven,'fro')^2+norm(xeven-barxodd,'fro')^2;

% Acknowledgements
%
% This code is based on an Splus-code written by G.P. Nason
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

 
