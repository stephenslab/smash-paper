function w = postblockmean(method,wcoef,l,sigma,thet0)

% postblockmean: Group the empirical wavelet coefficients into disjoint blocks 
%                of length l, and extract the wavelet coefficients needed by 
%                the RECBLOCKMEAN procedure.
% Usage
%                w = postblockmean(method,wcoef,l,sigma,thet0)
% Inputs
%   method	     'Augment' for augmenting the data, 'Truncate' for truncating 
%                the data
%   wcoef	     Empirical wavelet coefficients at a given resolution level
%   l		        Length of the blocks
%   sigma	     Level of noise
%   thet0	     Parameter vector : [PIE0,TAU0,RH0] where PIE0, TAU0 and RH0 
%                are the initial values for the maximization algorithm
% Outputs
%   w		        Denoised wavelet coefficients.
% Notes
%                The first few empirical wavelet coefficients might be re-used 
%                to fill the last block (which is called the "Augmented" case) 
%                or the last few remaining empirical wavelet coefficients might 
%                not be used in the posterior-based inference (which is called 
%                the "Truncated" case), should l not divide the length of wcoef 
%                exactly. 

%form k blocks of length l
n=length(wcoef);

%augment/truncate data if necessary
if rem(n,l)~=0,
  if strcmp(method,'Truncate')
    wcoefnew=wcoef(1:floor(n/l)*l);
  elseif strcmp(method,'Augment') 
    wcoefnew=[wcoef wcoef(1:l-rem(n,l))];
  else 
    disp(sprintf('MakeSignal: I don''t recognize <<%s>>',method))
    disp('Allowable methods are:')
       disp('Truncate'),
       disp('Augment'), break
  end
  nnew=length(wcoefnew);
else
  wcoefnew=wcoef;
  nnew=n;
end
wcoefnew=wcoefnew(:)';
blocks=reshape(wcoefnew,l,nnew/l)';

%parameter estimation
%uoptions=optimset('TolFun',1e-12,'TolX',1e-12);
%thet=fminsearch('mgmll',thet0,uoptions,sigma,blocks);
%thet=fminsearch('mgmll',thet,uoptions,sigma,blocks);
thet=fmins('mgmll',thet0,[],[],sigma,blocks);
thet=fmins('mgmll',thet,[],[],sigma,blocks);

pie=ilt(thet(1)); tau=abs(thet(2)); rho=atan(thet(3))*2/pi;

thet=[pie tau rho];

%initialisations
P=zeros(l);
for i=1:l
 for j=1:l
   P(i,j)=rho^abs(i-j);
 end
end
V=tau^2*P; 
A=inv(sigma^2*inv(V)+eye(l));

%operate on each block
for i=1:nnew/l 
  dtilda=A*blocks(i,:)';
  Omikron=(1-pie)/pie*sqrt(det(V)/sigma^(2*l)*det(A))*...
     exp(-blocks(i,:)*A*blocks(i,:)'/(2*sigma^2));
  threshval(:,i)=1/(1+Omikron)*dtilda;
end

threshold=threshval(:)';
if strcmp(method,'Truncate')
  w=[threshold wcoef(floor(n/l)*l+1:end)];
elseif strcmp(method,'Augment')
  w=threshold(1:n);
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
 