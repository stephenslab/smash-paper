function xreg=regpengcv(x,p)

% regpengcv:  Computes the regularization parameter for penalized least 
%             squares wavelet estimation by generalized cross-validation 
%             needed by the PENWAV procedure.
% Usage
%             xreg = regpengcv(x,p)
% Inputs
% 	x 		     Coefficients of the periodic wavelet transform of a a 1-d signal
% 	p 		     Lp norm (p=2 for least-squares)
% Outputs
%	xreg       Regularized wavelet coefficients.
% References
%             Amato, U. & Vuza, D.T. (1997). Wavelet approximation of a function 
%             from samples affected by noise. Rev. Roumanie Math. Pure Appl., 42, 
%             481-493.

n=length(x);
nn=n;
scale=0;
while nn>1,;
nn=nn/2; 
scale=scale+1;
end
ap=zeros(scale,1);
regpar=-10:0.5:+10;
regpar=10.^regpar;
istart=2 .^[0:scale-1]+1;
iend=2 .^([0:scale-1]+1);
istart(1)=1;
for i=1:scale
  ap(istart(i):iend(i))=2.^(2.*(i-1.)*p);
end
gcvmin=1.d+30;
cost=2.^(-2*scale*p);
for ir=1:length(regpar)
  rp=regpar(ir);
  if rp < cost 
    v2=ap ./ (1. + rp.*ap);
  else
    v2=1. - 1. ./(1. + rp.*ap); 
end
  v1=(x' .* v2).^2;
  v1(1)=0.;
  v2(1)=0.;
  gcv=sum(v1)/sum(v2)^2;
  if gcv <= gcvmin 
    gcvmin=gcv;
    rpgcv=rp;
  else
    break
  end
end

xreg=(x' ./ (1. + rpgcv * ap))'; 

% Copyright (c) 1997
%
% U. Amato & V.T. Vuza
