fignum=1;
figure(fignum)
n=256;
names={'Step' 'Wave' 'Blip' 'Blocks' 'Bumps' 'HeaviSine' ...
       'Doppler' 'Angles' 'Parabolas' 'Time Shifted Sine' ...
       'Spikes' 'Corner'};
xt=(1:n)/n;
signal=zeros(12,n);
nsignal=zeros(12,n);
 axis([0 1 0 1])
for test=1:12,
%    signal(test,:) = MakeSignalNewb(names(test), 256);
     nsignal(test,:) = noisysignal(names(test), 256, 3);
     axis([0 1 0 1])
	 subplot(4,3,test);
%     plot(xt,signal(test,:))
      plot(xt,nsignal(test,:),'+')
	  axis([0 1 0 1])
	 xlabel(char(names(test)));
end
gname=sprintf('Noisyallfunctions');
laprint(fignum,gname,'keepticklabels','asonscreen')

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
  