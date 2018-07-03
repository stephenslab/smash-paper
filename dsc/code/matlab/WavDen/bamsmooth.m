function f = bamsmooth(data,coarsest,filt,finest_lev)

% bamsmooth:    Extract wavelet coefficients from data needed by the RECBAMS 
%               procedure. It calls the BAYESRULE procedure.
% Usage
%               f = bamsmooth(data,coarsest,filt,finest_lev)
% Inputs
%   data		    1-d noisy signal, length(signal)= 2^J
%   coarsest	 Low-frequency cutoff for shrinkage
%   filt		    Quadature mirror filter for wavelet transform
%   finest_lev	 Highest resolution level
% Outputs
%   f           Vector of same length as "data" containing the estimate obtained 
%               by applying model-induced wavelet shrinkage on the wavelet 
%               coefficients. 

 data_wd = FWT_PO (data, coarsest, filt);
 data_smowd = data_wd;

 %----------setting the parameters: mu----------------------
 niz=data_wd(dyad(finest_lev));
 nizo=sort(niz); mm=length(nizo);
 low=floor(1/4 * mm); high=floor(3/4 * mm);
 pseudos = abs(nizo(high)-nizo(low))/1.5; %tukey's book
 mu = 1/pseudos^2;

 %----------setting pars: tau and epsilon (levelwise)-------
 shape=ones(1, finest_lev);
 tauu= sqrt((std(data).^2 - 1/mu));
 for i = finest_lev:-1:coarsest,
  eps(i) = 1 - 1/(i - coarsest +1)^1.5;
  tau(i) = tauu;
  
 %----------------------------------------------------------
  a = data_wd(dyad(i));
  b = 0 .* a;
  nn=length(a);
   for  j = 1:nn,
       [zz,  b(j)] = bayesrule(a(j), mu, tau(i), eps(i));
   end
  data_smowd(dyad(i))= b; 
 end
 
 f = IWT_PO( data_smowd, coarsest, filt);
 
 % Copyright (c) 2000
 %
 % Brani Vidakovic & Fabrizio Ruggeri 