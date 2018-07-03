function f = blockJS(method,wcoef,l,sigma)

% blockJS:   Group the empirical wavelet coefficients into disjoint blocks of 
%            length l, and extract the wavelet coefficients using a James-Stein 
%            rule needed by the RECBLOCKJS procedure.
% Usage
%            f = blockJS(method,wcoef,l,sigma)
% Inputs
%   method	 'Augment' for augmenting the data, 'Truncate' for truncating the data
%   wcoef	 Empirical wavelet coefficients at a given resolution level
%   l		    Length of the blocks
%   sigma	 Level of noise
% Outputs
%   f		    Denoised wavelet coefficients.
% Notes
%            The first few empirical wavelet coefficients might be re-used to fill 
%            the last block (which is called the "Augmented" case) or the last few 
%            remaining empirical wavelet coefficients might not be used in the 
%            inference (which is called the "Truncated" case), should l not divide 
%            the length of wcoef exactly. 

%form k blocks of length l
lambda = 4.50524;
n=length(wcoef); 
%augment/truncate data if necessary
if n < l & strcmp(method,'Truncate'),
  S = sum(wcoef.^2);
  f = max(0, 1-lambda*l*sigma^2/S)*wcoef;
else
if rem(n,l)~=0,
  if strcmp(method,'Truncate')
    wcoefnew=wcoef(1:floor(n/l)*l);
  elseif strcmp(method,'Augment') 
    wcoefnew=[wcoef wcoef(1:l-rem(n,l))];
  else 
    disp(sprintf('BlockJS: I don''t recognize <<%s>>',method))
    disp('Allowable methods are:')
    disp('Truncate')
    disp('Augment')
    return
  end
  nnew=length(wcoefnew);
else
  wcoefnew=wcoef;
  nnew=n;
end

wcoefnew=wcoefnew(:)';
blocks=reshape(wcoefnew,l,nnew/l)';

%operate on each block
for i=1:nnew/l
   S = sum(blocks(i,:).^2);
   threshval(:,i) = max(0, 1-lambda*l*sigma^2/S)*blocks(i,:)';   
end

threshold=threshval(:)';
if strcmp(method,'Truncate')
  f=[threshold wcoef(floor(n/l)*l+1:end)];
elseif strcmp(method,'Augment')
  f=threshold(1:n);
end
end

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
