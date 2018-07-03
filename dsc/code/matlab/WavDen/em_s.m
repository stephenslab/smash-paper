function f = em_s(thet,sigma,dat)

% em_s:  EM algorithm for mixture of univariate normals (special case).
%        EM_S(THET,SIGMA,DAT) estimates the parameters from a mixture 
%        of 2 univariate normals, viz 
%
%        pie * N(0,sigma^2+tau^2) + (1-pie) * N(0,sigma^2)
%
%   where SIGMA is known, using the EM-algorithm.
%
%   THET are the unknown parameters, in form, [pie tau].

%initialisations
pie=thet(1); tau2=thet(2)^2;
n=length(dat);

%posterior quantities
e=exp(-tau2*dat.^2/(2*sigma^2*(sigma^2+tau2)));
Omikron=(1-pie)/pie*sqrt(sigma^2+tau2)/sigma*e;
eta=1./(1+Omikron);

%first estimates
pie_hat=sum(eta)/n;
tau2_hat=max(0,sum(eta.*dat.^2)/sum(eta)-sigma^2);

while (abs(pie_hat-pie)+abs(tau2-tau2_hat))>1e-6
   pie=pie_hat;
   tau2=tau2_hat;
   e=exp(-tau2*dat.^2/(2*sigma^2*(sigma^2+tau2)));
   Omikron=(1-pie)/pie*sqrt(sigma^2+tau2)/sigma*e;
   eta=1./(1+Omikron);
   pie_hat=sum(eta)/n;
   tau2_hat=max(0,sum(eta.*dat.^2)/sum(eta)-sigma^2);
end

%final estimates
f=[pie_hat sqrt(tau2_hat)];

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
 