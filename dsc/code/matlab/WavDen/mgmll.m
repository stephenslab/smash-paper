function f = mgmll(thet,sigma,dat)

% mgmll:  Multivariate Gaussian mixture log-likelihood.
%   MGMLL(THET,SIGMA,DAT) calculates the log-likelihood of the
%   the l-variate Gaussian mixture
%
%   X ~ pie N(0 , B) + (1-pie) N(0 , sigma^2*I)
%
%   where B = sigma^2 * I + tau^2 * R
%   with r(i,j) = rho^|i-j|, i,j=1,..l.
%
%   DAT is an m-by-l matrix of observations
%   THET are the parameters, in form, [ PIE TAU RHO ].
%
%   NB. SIGMA is known.

%initialisations
[m,l]=size(dat);
pie=ilt(thet(1)); tau=abs(thet(2));
%rho=sin(thet(3));
rho=atan(thet(3))*2/pi;

%construct matrix B
R=zeros(l);
for i=1:l
 for j=1:l
   R(i,j)=rho^abs(i-j);
 end
end
B=sigma^2*eye(l)+tau^2*R;

%var-covariance matrices
Sigma1=B; Sigma2=sigma^2*eye(l);

pdf1=-l/2*log(2*pi)-1/2*log(det(Sigma1))-1/2*dat*inv(Sigma1)*dat';
pdf1=exp(diag(pdf1));
pdf2=-l/2*log(2*pi)-1/2*log(det(Sigma2))-1/2*dat*inv(Sigma2)*dat';
pdf2=exp(diag(pdf2));
pdf=pie*pdf1+(1-pie)*pdf2;

%log-likelihood
f=-sum(log(pdf));

% Copyright (c) 2001:
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
