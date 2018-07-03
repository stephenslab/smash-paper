function thr = mingcv(x)

% mingcv:   Determining the generalized cross-validated threshold.
%           It calls the GCV procedure.
% Usage
%           thr = mingcv(x)
% Input
%   x	      vector containing input wavelet coefficients 
%           that will be shrunk with threshold thr
% Ouput
%   thr     generalized cross-validated threshold

%initialisatie Fibonacci - getallen
 F0 = 1;
 F1 = 1;
 F2 = 2;
 aantal = 1;
 
 %berekening van het benodigd aantal iteratiestappen
 m = 50;
 b = max(abs(x));
 a =  b / m; %(in 0 kunnen we gcv niet berekenen)
  eps = 0.0001;
  while (b - a) / F2 > eps
    F0 = F1;
    F1 = F2;
    F2 = F0 + F1;
    aantal = aantal + 1;
  end

  firstmove = 0;
  while (firstmove == 0)

       v = a + (F1/F2) * (b - a);
       fv = gcv(x,v);
       u = a + (F0/F2) * (b - a);
       fu = gcv(x,u);
     
       for k = 1:aantal - 1
        if (fu > fv)
           a = u;
           u = v;
           fu = fv;
           v = b-u+a;
           fv = gcv(x,v);
       firstmove = 1;
        else
           b = v;
           v = u;
           fv = fu;
           u = a+b-v;
           fu = gcv(x,u);
        end
       end
     
       if (b < 10^(-10))
         firstmove = 1;
       end % if
     
       a = b/m;
     
  end  %while (firstmove == 0)


  thr = b;


% Acknowledgements
%
% This code is based on a Matlab Code written and kindly provided by M. Jansen
%
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

 
