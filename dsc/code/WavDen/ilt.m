function f=ilt(thet)

% ILT:  Inverse logistic transform
%
% ILT(THET) calculates the inverse of the logistic transform, ie 
%
%                                       1
%                     ilt(thet) = --------------
%                                 1 + exp(-thet)
%

f=1./(1+exp(-thet)); 