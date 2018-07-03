function sig = MakeSignalNewb(signal_name,n)

% MakeSignalNewb:   Make artificial 1-d signal
% 	                 sig = MakeSignalNewb(signal_name,n)
% Inputs 
%  signal_name      string: 'Step','Wave','Blip','Blocks','Bumps',
%					     'HeaviSine','Doppler','Angles',
%					     'Parabolas','Time Shifted Sine','Spikes','Corner'
% 	n 			        Desired signal length
% Ouputs
%	sig			     Artificial 1-d signal.

	t = (1:n)./n;
	if strcmp(signal_name,'Step'),
		sig = 0.2 + 0.6*(t > 1/3 & t <= 0.75); 

	elseif strcmp(signal_name,'Wave'),
		sig = 0.5 + (0.2.*cos(4*pi*t)) + (0.1.*cos(24*pi*t));
 
	elseif strcmp(signal_name,'Blip'),
		sig = (0.32 + (0.6.*t) + 0.3*exp(-100*((t-0.3).^2))).*(t <= 0.8) + ...
         	  (-0.28 + (0.6.*t) + 0.3*exp(-100*((t-1.3).^2))).*(t > 0.8);
	elseif strcmp(signal_name,'Blocks'),
     	pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
     	hgt = [4 (-5) 3 (-4) 5 (-4.2) 2.1 4.3  (-3.1) 2.1 (-4.2)];
     	sig = 2*ones(size(t));
     	for j=1:length(pos)
       	 sig = sig + (1 + sign(t-pos(j))).*(hgt(j)/2) ;
     	end
     	sig = (0.6/9.2)*sig + 0.2;
  
	elseif strcmp(signal_name,'Bumps'),
    	pos = [ .1 .13 .15 .23 .25 .40 .44 .65  .76 .78 .81];
    	hgt = [ 4  5   3   4  5  4.2 2.1 4.3  3.1 5.1 4.2];
    	wth = [.005 .005 .006 .01 .01 .03 .01 .01  .005 .008 .005];
    	sig = zeros(size(t));
    	for j =1:length(pos)
       		sig = sig + hgt(j)./(( 1 + (abs(t - pos(j))./wth(j))).^4);
   		 end 
    	sig = ((0.6/5.3437952)*sig) + 0.2;

	elseif strcmp(signal_name,'HeaviSine'),
		sig = 4.*sin(4*pi.*t) - sign(t - .3) - sign(.72 - t) + 5;
    	sig = (0.6/9)*sig + 0.2;

	elseif strcmp(signal_name,'Doppler'),
    	sig = sqrt(t.*(1-t)).*sin((2*pi*1.05) ./(t+.05)) + 0.5;
    	sig = 0.6*sig + 0.2;

	elseif strcmp(signal_name,'Angles'),
     sig = ((2*t + 0.5).*(t <= 0.15)) + ...
        ((-12*(t-0.15) + 0.8).*(t > 0.15 & t <= 0.2)) + ...
        0.2*(t > 0.2 & t <= 0.5) + ...
        ((6*(t - 0.5) + 0.2).*(t > 0.5 & t <= 0.6)) + ...
        ((-10*(t - 0.6) + 0.8).*(t > 0.6 & t <= 0.65)) + ...
        ((-0.5*(t - 0.65) + 0.3).*(t > 0.65 & t <= 0.85)) + ...
        ((2*(t - 0.85) + 0.2).*(t > 0.85));
	elseif strcmp(signal_name,'Parabolas'),
     pos = [0.1 0.2 0.3 0.35 0.37 0.41 0.43 0.5 0.7 0.9];
     hgt = [(-30) 60 (-30) 500 (-1000) 1000 (-500) 7.5 (-15) 7.5];
     sig = zeros(size(t));
	 for j =1:length(pos)
       sig = sig + hgt(j).*((t-pos(j)).^2).*(t > pos(j));
     end
    	sig = sig + 0.8;

	elseif strcmp(signal_name,'Time Shifted Sine') | ...
           strcmp(signal_name,'TSh Sine'),
      u = t;
    	for j =1:4,
        u = 0.5*(1-cos(pi*u));
    	end
        sig = 0.3*sin(3*pi*(u+t)) + 0.5;
   
	elseif strcmp(signal_name,'Spikes'),
		sig = 15.6676 .* (exp(-500.*(t-0.23).^2) + ...
		2.*exp(-2000.*(t-0.33).^2) +  4.*exp(-8000.*(t-0.47).^2) + ...
		3.*exp(-16000.*(t-0.69).^2) +exp(-32000.*(t-0.83).^2) );
		sig=(0.6/range(sig)).*sig+0.2;
		
	elseif strcmp(signal_name,'Corner'),
		sig = t;
		sig(t <= 0.5) = 62.387.*10.*t(t <= 0.5).^3.*(1-4.*t(t <= 0.5).^2);
		sig(0.5 < t & t <= 0.8) = 62.387.*3.*(0.125-t(0.5 < t & t <= 0.8).^3).*t(0.5 < t & t <= 0.8).^4;
		sig(t > 0.8) = 62.387.*59.443.*(t(t > 0.8)-1).^3;
		sig=(0.6/range(sig)).*sig+0.6;
	else
	   disp('MakeSignalNewb: I don*t recognize');signal_name
       disp('Allowable Names are:')
       disp('Step'),
       disp('Wave'),
       disp('Blip'),
       disp('Blocks'),
       disp('Bumps'),
       disp('HeaviSine'),
       disp('Doppler'),
       disp('Angles'),
       disp('Parabolas'),
       disp('Time Shifted Sine'),
       disp('Spikes'),
       disp('Corner')
    end
    
% Acknowledgement
%
% This code is based on a code provided by Buckheit, Chen, Donoho,
% Johnstone & Scargle.
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
 