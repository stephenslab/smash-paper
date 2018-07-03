function f = em2(thet,datall,lev)

% em2:  EM algorithm for mixture of univariate normals (special case).
%       EM(THET,DATALL,DATLEV) estimates the parameters from a mixture 
%       of 2 univariate normals, viz 
%
%       pie * N(0,sigma^2+tau^2) + (1-pie) * N(0,sigma^2)
%
%       using the EM-algorithm.
%
%   THET are the unknown parameters, in form, [pie tau sigma].

%initialisations
n=length(datall);
J=log2(n);
pie=1:(J-1); tau2=1:(J-1);
pie_hat=1:(J-1); tau2_hat=1:(J-1);
pie(lev:(J-1))=thet(1); tau2(lev:(J-1))=thet(2)^2; sigma2=thet(3)^2;

% overall posterior quantities
eall=1:n;
Omikronall=1:n;
for j = lev:(J-1),
  eall(2^j+1:2^(j+1))=exp(-tau2(j)*datall(2^j+1:2^(j+1)).^2 ...
  /(2*sigma2*(sigma2+tau2(j))));
  Omikronall(2^j+1:2^(j+1))=(1-pie(j))/pie(j)*sqrt(sigma2+tau2(j)) ...
  /sqrt(sigma2)*eall(2^j+1:2^(j+1));
end
etaall=1./(1+Omikronall);

%first estimates
sigma2_hat=sum((1-etaall(2^lev+1:2^J)).*(datall(2^lev+1:2^J).^2)) ...
/(2^J-2^lev-sum(etaall(2^lev+1:2^J)));
for j =lev:(J-1),
  pie_hat(j)=sum(etaall(2^j+1:2^(j+1)))/2^j;
  tau2_hat(j)=max(0,sum(etaall(2^j+1:2^(j+1)).*datall(2^j+1:2^(j+1)).^2) ...
  /sum(etaall(2^j+1:2^(j+1)))-sigma2_hat);
end

while max(abs(pie_hat-pie)) + max(abs(tau2-tau2_hat)) + abs(sigma2-sigma2_hat) >1e-3 ,
   pie=pie_hat;
   tau2=tau2_hat;
   sigma2=sigma2_hat;
   
   for j = lev:(J-1),
     eall(2^j+1:2^(j+1))=exp(-tau2(j)*datall(2^j+1:2^(j+1)).^2 ...
     /(2*sigma2*(sigma2+tau2(j))));
     Omikronall(2^j+1:2^(j+1))=(1-pie(j))/pie(j)*sqrt(sigma2+tau2(j)) ...
     /sqrt(sigma2)*eall(2^j+1:2^(j+1));
   end
   etaall=1./(1+Omikronall);
   
   sigma2_hat=sum((1-etaall(2^lev+1:2^J)).*datall(2^lev+1:2^J).^2) ...
   /(2^J-2^lev-sum(etaall(2^lev+1:2^J)));
   for j =lev:(J-1),
     pie_hat(j)=sum(etaall(2^j+1:2^(j+1)))/2^j;
     tau2_hat(j)=max(0,sum(etaall(2^j+1:2^(j+1)).*datall(2^j+1:2^(j+1)).^2) ...
     /sum(etaall(2^j+1:2^(j+1)))-sigma2_hat);
   end
end

%final estimates
f=struct('pie',pie_hat,'tau',sqrt(tau2_hat),'sigma',sqrt(sigma2_hat));

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

 