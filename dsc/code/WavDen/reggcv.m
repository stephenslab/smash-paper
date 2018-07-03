function xreg=reggcv(x,p,lmin)

% regpengcv:  Computes the regularization parameter for penalized least 
%             squares wavelet estimation by generalized cross-validation 
%             needed by the PENWAV procedure.
% Usage
%             xreg = reggcv(x,p,lmin)
% Inputs
% 	x 		     Coefficients of the periodic wavelet transform of a a 1-d signal
% 	p 		     Lp norm (p=2 for least-squares)
% 	lmin 	     lowest resolution level above which penalization is applied
% Outputs
%	xreg 	     Regularized wavelet coefficients.
% References
%             Amato, U. & Vuza, D.T. (1997). Wavelet approximation of a function 
%             from samples affected by noise. Rev. Roumanie Math. Pure Appl., 42, 
%             481-493.

n=size(x,1);
xreg=x;
ap=zeros(n,1);
nn=n;
scale=0;
while nn>1,;
nn=nn/2; 
scale=scale+1;
end
regpar=-10:0.5:+10;
regpar=10.^regpar;
istart=[2.^(linspace(1,scale,scale)-1)+1];
istart(1)=1;
iend=[2.^linspace(1,scale,scale)];
if lmin >= scale-1 return, end % GCV not effective, reg. par. = 0
if lmin > 0 ap(1:iend(lmin),1)=2^(2*max([0 lmin])*p); end
for i=lmin+1:scale
  ap(istart(i):iend(i))=2.^(2.*(i-1)*p);
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
  v1=(x .* v2).^2;
  gcv=sum(v1)/sum(v2)^2;
  if gcv <= gcvmin 
    gcvmin=gcv;
    rpgcv=rp;
  else
%    break
  end
end
xreg=x ./ (1.0 + rpgcv * ap) %%%%% NOTE NO FACTOR AT MOMENT %%%%%%

% Copyright (c) 1997
%
% U. Amato & V.T. Vuza

  