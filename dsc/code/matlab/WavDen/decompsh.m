function [y,yd,ys,ns, beta]  = decompsh(x,L,qmf,Lmax,show)

% decompsh:  Signal is decomposed by an orthogonal wavelet tranform:  w = Tx
%            The tranformed signal is modeled as:  w = beta + epsilon, where 
%            epsilon is a Gaussian noise. The expected denoised component b is 
%            decomposed as: beta = mu + eta, where mu is the deterministic 
%            component and eta the stochastic component. This procedure is called 
%            by the RECDECOMPSH procedure.
% Usage 
%            [y,yd,ys,ns, beta] = decompsh(x,L,qmf,Lmax,show)
% Inputs
%    x       1-d signal. length(x)= 2^J
%    L       Level of coarsest scale L
%            1<L<<J
%    qmf     Quadrature Mirror Filter for Wavelet Transform
%            Optional, Default = Symmlet 8
%    Lmax    Maximum level to retain in the reconstruction
%            (wcoef(L>Lmax) = 0),  L  <Lmax <J
%            Lmax(Default)=J-2:    finest details are supposed to be only noise
%    show    if 1 then show figures showing the process, if 0 don't
% Outputs 
%    y       Final (stochastic + deterministic) denoised estimate
%    yd      Deterministic structure estimate
%    ys      Stochastic component estimate
%    ns      Noise component
%    beta    Wavelet transform of estimate
% Notes
%            Variance estimation in finest scale is not done with the 
%            semivariogram method, but with traditional median absolute
%            deviation.
% REFERENCES
% Huang, H.-C. and Cressie, N. (2000). 
% Deterministic/stochastic wavelet decomposition for recovery of signal from noisy data. 
% Technometrics, 42, 262-276.
%
% ACKNOWLEDGEMENT
% This function was contributed by Camilo La Rota (Lab TIMC, Grenoble). 

x  = ShapeAsRow(x);
[n,J] = dyadlength(x) ;

  if (nargin < 2)|(isempty(L)),
      L=3;
  end
  if ((nargin < 3)|(isempty(qmf))),
      qmf = makeonfilter('Symmlet',8);
  end
  if ((nargin < 4) | (isempty(Lmax)))
      Lmax =J-2;  %J-1 level is supposed to be only noise for biological signals
  end
  if (nargin < 5)
      show=1;  %show decomposition process 
  end
%

%
%-------------------Wavelet Tranform ------------------------
  wcoef         = FWT_PO(x,L,qmf) ;
%
% noise variance estimation with high resolution wavelet coefficients
  wfine    = wcoef(dyad(J-1));
  %noisevar = varvar(wfine);  % variance estimation with Huang method is not as robust as they say ...
  noisevar=(mad2(wfine))^2;  % estimate with donoho method, 

  if show 
       fig1=figure; 
       set(fig1,'Name','Deterministic Component Estimation','Units','normalized','Position',[0.05,0.55,0.4,0.4]); 
  end          

% ------------------------mean structure computation (mu)------------------------------------------------


  mu            = wcoef;
  
  %cancel high resolution coefficients dues to noise
  for j=J-1:-1:Lmax+1,
        mu(dyad(j))       = zeros(size(dyad(j)));
       if show subplot(221); plot(1:n,x,'b',1:n,IWT_PO(mu,L,qmf),'r');title(['Aproximation up to level ' num2str(j)]);axis tight; pause;end         
  end

  for j=Lmax:-1:L,             %loop through scales
   wj= wcoef(dyad(j));

   if show   %show quantile plot
     subplot(2,2,2);qnqplot(wj);title('qnqplot and estimated threshold'); axis tight;
   end          

   nj= length(wj);
   [wjs,l]=sort(wj);
   qj= norminv((1:nj)/(nj+1),0,1);
   Tj                = median(abs(wj))/0.6745;   %slope of quantiles-normalquantiles plot
   qwj(l)            = qj;                 
   resj              = wj-(Tj*qwj);              %residual of wavelet coefficient to corresponding quantile
   for i=1:nj
      if (abs(wj(i))>=abs(Tj*qwj(i)))
           qwj(i)=0;
      end
   end
   lambdaj         = Tj*max(abs(qwj));          %estimation of threshold for scale j

   if show   %show quantile plot
     t=axis;
     t=t(1):(t(2)-t(1))/100:t(2);
     thres=lambdaj*ones(size(t));
     hold;plot(t,thres,'r:',t,thres*(-1),'r:'); hold;
   end          
   
   muj       = ((resj.^2)./(Tj^2 + resj.^2)).*wj;      %deterministic shrinkage, determinstic component estimation


    for k = 1:2^j,
     if abs(wj(k)) <= lambdaj
        muj(k)       = 0;         %adaptive hardthreshold of deterministic component
     end
    end

    mu(dyad(j))       = muj;

    if show 
       subplot(2,1,2);
       plot(1:nj,wj,'b*',1:nj,muj,'ko');title(['Wavelet coefficients shrinkage, Level:' num2str(j)]);axis tight;hold;  
       plot(1:nj,muj,'k');    %show estimated deterministic component for scale j
       t=axis;
       t=t(1):(t(2)-t(1))/100:t(2);
       thres=lambdaj*ones(size(t));
       plot(t,thres,'r:',t,thres*(-1),'r:'); hold;
    end

   clear qwj resj muj; 
   if show subplot(221); plot(1:n,x,'b',1:n,IWT_PO(mu,L,qmf),'r');title(['Aproximation up to level ' num2str(j)]);axis tight; pause; end         

 end       %end loop trough scales

  yd    = IWT_PO(mu,L,qmf);          % intermediate result, deterministic component
  yd   = ShapeLike(yd,x);
%---------------------------end deterministic estimation-----------------------------
  if show 
       fig2=figure; 
       set(fig2,'Name','Stochastic Component Estimation','Units','normalized','Position',[0.55,0.55,0.4,0.4]); 
       subplot(2,1,1); plot(1:n,x,'b',1:n,yd,'r');title('Estimated Deterministic Component'); axis tight; 
  
       fig3=figure; 
       set(fig3,'Name','Estimated components','Units','normalized','Position',[0.55,0.05,0.4,0.4]); 
       subplot(2,1,1); plot(1:n,x,'b',1:n,yd,'k');title('Estimated Deterministic Component'); axis tight; pause;
  end          

  

%  --------------Stochastic component estimation (eta),(scale independence hypothesis, decompI)----------------------
  if show 
  end          
  beta = mu;

  for j=Lmax:-1:L,        %loop trough scales
     muj = mu(dyad(j));
     nj= length(muj);
     %stochastic process covariance estimation
     Dj       = wcoef(dyad(j)) - muj;
     Dj2      = Dj*Dj';
     Dvar     =  (Dj2)/(2^j) - noisevar;
     sigma2j = max(Dvar,0);

    %  ----------Bayesian shrinkage, estimation with stochastic and mean components-----------
     dbj          = (sigma2j/(sigma2j+noisevar)) * Dj;
     beta(dyad(j)) = muj + dbj;
    
     if show 
          figure(fig2); 
          subplot(2,1,1); plot(1:n,x,'b',1:n,yd,'r',1:n,IWT_PO(beta,L,qmf),'k');title(['Deterministic+Stochastic Aproximation up to level ' num2str(j)]); axis tight; 
          

          subplot(2,1,2);plot(1:nj,wcoef(dyad(j)),'b:' ,1:nj,muj,'r:'); axis tight; hold; 
          plot(1:nj,wcoef(dyad(j)),'b*' ,1:nj,muj,'ro'); axis tight;  
          plot(1:nj,beta(dyad(j)),'k+',1:nj,beta(dyad(j)),'k' ); hold; 
          title(['Stochastic wavelet coefficient shrinkage at level:' num2str(j)]);
          pause;
     end
     
     clear muj;
     clear Dj;
  end                 %end loop trough scales


 %------------------- end stochastic component estimation ------------------------ 
   

% Inverse transformation
  y    = IWT_PO(beta,L,qmf);      
  y    = ShapeLike(y,x);
  ys   = y-yd;
  ns   = y-x;
  
 if show 
   figure(fig3);
   subplot(2,1,1);plot(1:n,x,'b',1:n,yd,'r:',1:n,y,'k');title('Final Estimation (Deterministic+Stochastic)'); axis tight
   subplot(2,1,2);plot(1:n,ns,'g:',1:n,yd,'b',1:n,ys,'k');title('Estimated components'); axis tight
   pause;
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
