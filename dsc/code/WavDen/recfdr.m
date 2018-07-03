function f = recfdr(signal,type,h)

% recfdr:   Smoothing using the false discovery rate method.
% Usage
%           f = recfdr(signal,type,h)
% Inputs
%   signal	1-d Noisy signal, length(signal)= 2^J
%   type	   'S' for soft thresholding, 'H' for hard thresholding.
%		      Optional, default = hard thresholding
%   h 		Quadrature Mirror Filter for Wavelet Transform
%		      Optional, default = Symmlet 8
% Outputs
%   f		   Estimate, obtained by applying thresholding on the wavelet 
%           coefficients.
% References
%           1. Abramovich, F. & Benjamini, Y. (1995). Thresholding of 
%              wavelet coefficients as multiple hypotheses testing procedure. 
%              In "Wavelets and Statistics", Antoniadis, A. & Oppenheim, G. 
%              (Eds.), Lecture Notes in Statistics 103, pp. 5-14, New York: 
%              Springer-Verlag. 
%           2. Abramovich, F. & Benjamini, Y. (1996). Adaptive thresholding 
%              of wavelet coefficients. Comput. Statist. Data Anal., 22, 
%              351-361.

if nargin < 2,
	type = 'H';
end

if nargin < 3,
	h = MakeONFilter('Symmlet',8);
end

%initialisations
n=length(signal);
lev=floor(log2(log(n)))+1;

%perform DWT
d=FWT_PO(signal,0,h);
J=log2(n);

%estimate sigma
highestlev=d(2^(J-1)+1:2^J);
sigma=median(abs(highestlev))/0.6745;

%calculate the threshold
alpha=0.05;
minit=length(d(2:n));
dinit=d(2:n);
thinit=norminv(1-alpha/2)*sigma;
if J>12,
	ninit=3;
elseif J>10,
	ninit=2;
else
	ninit=1;
end

for k=1:ninit,
	dinit1 = dinit(abs(dinit) >= thinit);
	minit=length(dinit1);
	if (minit == 0)
	  thresh = max(abs(d))*1.0001;
	else
	  thinit = norminv(1-(alpha*minit)/(2*n))*sigma;
	  minit1 = length(dinit1(abs(dinit1) >= thinit));
	  dinit = dinit1;
	end
end

if (sigma > 0),
	m = length(d);
	minit = length(dinit);
	p = 2*(1-normcdf(abs(dinit)/sigma));
	[psort,index] = sort(p);
	j = 1:minit;
	m0 = max(j(p(index) <= (alpha*j)/m));
	if (not(isnan(m0)) & m0 <= minit),
	  thresh = abs(dinit(index(m0)));
  	elseif isnan(m0),
	  thresh = max(abs(dinit))*1.0001;
	else
	  thresh = 0;
	end
else
	thresh = 0;
end

%Extraction of the wavelet coefficients
wcoef=FWT_PO(signal,lev,h);
if strcmp(type,'H'),
  wcoef(2^lev+1:2^J) = HardThresh(wcoef(2^lev+1:2^J),thresh);
else
  wcoef(2^lev+1:2^J) = SoftThresh(wcoef(2^lev+1:2^J),thresh);
end

f=IWT_PO(wcoef,lev,h); 

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
  
