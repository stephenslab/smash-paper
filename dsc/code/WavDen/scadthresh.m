function xscad = scadthresh(x,t) 

% This function is needed by the RECSCAD procedure.

tau=ones(size(x))*t;

xscad=sign(x).*max(abs(x)-tau,0).*(abs(x)<= 2*tau)...
+ (abs(x) <= 3.7*tau) .* (abs(x) > 2*tau) .* (2.7*x-3.7.*tau.*sign(x))./1.7 ...
+  (abs(x) > 3.7*tau).*x;

% Copyright (c) 1999
%
% Anestis Antoniadis & Jianqing Fan
  