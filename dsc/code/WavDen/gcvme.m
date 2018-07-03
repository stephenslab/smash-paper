function vofd = gcvme(c,n,sigma,Sy)

% This is a generalized cross-validation procedure for mixed-effects
% models needed by the RECMIXED procedure.

format long e

betamnoise=(Sy).^2-c*sigma^2*log(n);
betanew=(abs(betamnoise)+betamnoise)/2;
betanew=betanew./((Sy).^2);
vbeta=(1.-betanew);
vvbeta=vbeta.*Sy;
vnew=sum(vvbeta.^2)*n;

down=(sum(betanew==0)+sum((betanew-1).*(betanew~=0)))^2 ;

vofd=(vnew/down)/(sigma^2) ;

% Copyright (c) 1999
%
% S.Y. Huang & H.H.-S Lu


 