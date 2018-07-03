% Wavelet denoising example on an electrical consumption signal.
%
% This signal has been thoroughly analyzed in "Misiti, Misiti, Oppenheim \& Poggi (1994).
% Decomposition en ondelettes et methodes comparatives: etude d'une courbe de charge electrique.
% Revue de Statistique Appliquee, 17, 57--77".
%
% It is particularly interesting because of noise introduced whenever a defect is
% present in the  monitoring equipment. The data consist of measurement of a complex, 
% highly-aggregated plant: the electrical load consumption, sampled minute by minute, 
% over a 3-day period.
% 
% This example loads the original 1-D signal, chooses a portion of it, denoises it
% with two denoising procedures and plots the results.
% 
% Load the signal  
s=file2var('eleccum.dat'); 

% Choose a portion of this signal.
x=(1:4096);  
signal=s(x);

% Plot the original signal
plot(x,signal); axis([-20 4097 100 600]); title('Electrical Signal'); 
pause

% Denoise the plotted portion of the signal using some of the procedures
% implemented in our library. Two such cases are illustrated below.

% Denoise using the Translation Invariant procedure with soft thresholding.
f = recTI(signal,'S');

% Denoise using the Neighblock procedure.
g=recneighblock(signal);

% Plot the results.
subplot(3,1,1);
plot(x,signal); axis([-20 4097 100 600]); title('Electrical Signal');
subplot(3,1,2);
plot(x,f); axis([-20 4097 100 600]); title('Soft TI Denoising')
subplot(3,1,3);
plot(x,g); axis([-20 4097 100 600]); title('Neighblock Denoising')

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

