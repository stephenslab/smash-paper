function f = simulations(test,M,n,rsnr,thet0s,thet0s2,thet0b,hint,show)

% simulations:  The main routine to produce the figures and tables given in
%               the paper by Antoniadis, Bigot & Sapatinas (2001). We refer 
%               to this paper which discusses in detail wavelet shrinkage and 
%               wavelet thresholding estimators in nonparametric regression. 
%               Various estimators have been compared in an extensive simulation 
%               on a variety of sample sizes, test functions, signal-to-noise
%               ratios and wavelet filters.
% Usage
%               f = simulations(test,M,n,rsnr,thet0s,thet0s2,thet0b,hint,show)
% Inputs
%   test        Integer pointing to Signal function (see comments in MakeSignalNewb.m)
%			       'Step', 'Wave', 'Blip', 'Blocks', 'Bumps', 'HeaviSine', 'Doppler', 
%               'Angles', 'Parabolas', 'Time Shifted Sine', 'Spikes', 'Corner'
%   M           Total number of simulations per setting
%   n           Size of the sample (number of data points)
%   rsnr  	    Root of signal-to-noise ratio
%   thet0s	    Initial parameter estimates for single : thet0s = [PIE0 TAU0]
%			       when Sigma is estimated through MAD
%   thet0s2     Initial parameter estimates for single : thet0s2 = [PIE0 TAU0 SIGMA0]
%			       when Sigma is part of the EM algorithm
%   thet0b	    Initial parameter estimates for blocks
%   hint 	    Quadrature mirror filter for wavelet transform
%		          (Enter 1 if you want Symmlet 8; Enter 2 if you want Coiflet 3)
%   show	       Enter 1 if you want figures showing boxplots of the various criterias 
%               to measure 
% Outputs      
%   A structure f with the following fields:
%   amse_f      Average Mean Squared Error
%   rmse_f	    Root Mean Squared Error
%   rmsb_f      Root Mean Squared Bias
%   amxdv_f	    Average Maximum Deviation
%   band_f	    Average width of the 2 SD band
%   aL1norm_f	 Average L1-Norm
%   aCPUtime_f  Average CPU Time
% References
%               Antoniadis, A., Bigot, J. & Sapatinas, T. (2001). Wavelet estimators 
%               in nonparametric regression: description and simulative comparison. 
%               Techical Report, Laboratoire IMAG-LMC, University Joseph Fourier, 
%               France.
% Notes    
%               Wavelet regression methods used :
%
%     1 - Classical Wavelet methods
%     ------------------------------ 
% 		1  VisuShrink (Hard) 		        : recvisu.m
% 		2  VisuShrink (Soft) 		        : recvisu.m
% 		3  SureShrink (Classical)          : recsure.m
% 		4  SureShrink (Hybrid Version)     : rechybsure.m
% 		5  TI (Hard) 			              : recTI.m
% 		6  TI (Soft) 			              : recTI.m
%		7  Minimax (Hard)		              : recminimax.m
%		8  Minimax (Soft)		              : recminimax.m
% 		9  CV (Nason's CV with a Hard)     : reccv.m
% 		10 CV (Nason's CV with a Soft)     : reccv.m
% 		11 NeighBlock (Cai and Silverman)  : recneighbl.m
% 		12 BlockJS (Cai, Augment) 	        : recblockJS.m
% 		14 BlockJS (Cai, Truncate) 	     : recblockJS.m
% 		14 Threshda1 (Ogden and Parzen)	  : recthreshda1.m
%		15 FDR (Hard)			              : recfdr.m
%		16 FDR (Soft)			              : recfdr.m
%		17 Least Squares Penalized	        : recpenwav.m
%		18 SCAD				                 : recscad.m
% 
%		2 - Bayesian Methods
% 		------------------------------
%		19 Mixed-effects model (Huang and Lu)                   : recmixed.m
% 		20 Decompsh (Huang and Cressie)	                       : recdecompsh.m
% 		21 Blocking Posterior Median (Augment)	                 : recblockmed.m
%		22 Blocking Posterior Median (Truncate)	              : recblockmed.m
% 		23 Blocking Posterior Median (Hybrid Version, Augment)  : rechybblockmed.m
% 		24 Blocking Posterior Median (Hybrid Version, Truncate) : rechybblockmed.m
% 		25 Blocking Posterior Mean (Augment)                    : recblockmean.m
%		26 Blocking Posterior Mean (Truncate)                   : recblockmean.m
% 		27 Blocking Posterior Mean (Hybrid Version, Augment)    : rechybblockmean.m
% 		28 Blocking Posterior Mean (Hybrid Version, Truncate)   : rechybblockmean.m
% 		29 Single Posterior Median
% 		(Silverman and Johnstone : Sigma is estimated through MAD) 	
%						                                            : recsinglemed.m
% 		30 Single Posterior Median (Sigma is part of the EM algorithm) 	
%						                                            : recsinglemed3.m
% 		31 Single Posterior Mean
% 		(Clyde and George : Sigma is estimated through MAD) 		
%						                                            : recsinglemean.m
% 		32 Single Posterior Mean
% 		(Clyde and George : Sigma is part of the EM algorithm) 		
%						                                            : recsinglemean3.m
% 		33 Single Hypothesis Testing (Vidakovic)                : recsinglehyp.m			                                  : recsinglehyp.m
%		34 Bams (Vidakovic and Rugeri) 	                       : recbams.m
%
%     (set to 1 in the case of real data)
  
names={'Step' 'Wave' 'Blip' 'Blocks' 'Bumps' 'HeaviSine' ...
       'Doppler' 'Angles' 'Parabolas' 'Time Shifted Sine' ...
       'Spikes' 'Corner'};
	   
if exist('thet0s')==0 thet0s=[0.5 1]; end;
if exist('thet0s2')==0 thet0s2=[0.5 1 0.7]; end;	
if exist('thet0b')==0 thet0b=[0.5 1 0.7]; end;	
if exist('n') ==0 n=64; end  
if exist('M') ==0 M=1; end  
if exist('rsnr') ==0 rsnr=7; end	
if exist('test') ==0 test=1; end
if exist('hint')==0 hint=1; end
if hint==1,
	h=MakeONFilter('Symmlet',8);
else
	h=MakeONFilter('Coiflet',3);
end
if exist('show')==0 show=1; end

fprintf('\n n=%4i test=%3i h=%1i M=%3i show=%3i rsnr=%3i\n',...
   n,test,hint,M,show,rsnr);
fname=sprintf('t%in%ih%ir%i.out',test,n,hint,rsnr); %filename containing the output
fidname=fopen(fname,'w');
fprintf(fidname,'\n sample size = %5i \n', n);
fprintf(fidname,'\n test function = %5i \n', test);
fprintf(fidname,'\n Number of Simulations = %5i \n', M);
fprintf(fidname,'\n QMF = %5i \n', hint);
fprintf(fidname,'\n rsnr = %5i \n', rsnr);
fprintf(fidname,'\n\n');
% 
%
% Initializations for Simulations

% MSE
mse_visu_H = zeros(M, 1);
mse_visu_S = zeros(M, 1);
mse_sure = zeros(M, 1);
mse_hybsure = zeros(M, 1);
mse_TI_H = zeros(M, 1);
mse_TI_S = zeros(M, 1);
mse_minimax_H = zeros(M, 1);
mse_minimax_S = zeros(M, 1);
mse_fdr_H = zeros(M, 1);
mse_fdr_S = zeros(M, 1);
mse_cv_H = zeros(M, 1);
mse_cv_S = zeros(M, 1);
mse_neighbl = zeros(M, 1);
mse_blockJS_A = zeros(M, 1);
mse_blockJS_T = zeros(M, 1);
mse_penwav = zeros(M, 1);
mse_scad = zeros(M, 1);
mse_mixed = zeros(M, 1);
mse_decompsh = zeros(M, 1);
mse_thrda1 = zeros(M, 1);
mse_bams = zeros(M, 1);
mse_blmed_A = zeros(M, 1);
mse_blmean_A = zeros(M, 1);
mse_blmed_T = zeros(M, 1);
mse_blmean_T = zeros(M, 1);
mse_singlmed = zeros(M, 1);
mse_singlmed2 = zeros(M, 1);
mse_singlmean = zeros(M, 1);
mse_singlmean2 = zeros(M, 1);
mse_singlhyp = zeros(M, 1);
mse_hybmed_A = zeros(M, 1);
mse_hybmean_A = zeros(M, 1);
mse_hybmed_T = zeros(M, 1);
mse_hybmean_T = zeros(M, 1);

% MXDV

mxdv_visu_H = zeros(M, 1);
mxdv_visu_S = zeros(M, 1);
mxdv_sure = zeros(M, 1);
mxdv_hybsure = zeros(M, 1);
mxdv_TI_H = zeros(M, 1);
mxdv_TI_S = zeros(M, 1);
mxdv_minimax_H = zeros(M, 1);
mxdv_minimax_S = zeros(M, 1);
mxdv_fdr_H = zeros(M, 1);
mxdv_fdr_S = zeros(M, 1);
mxdv_cv_H = zeros(M, 1);
mxdv_cv_S = zeros(M, 1);
mxdv_neighbl = zeros(M, 1);
mxdv_blockJS_A = zeros(M, 1);
mxdv_blockJS_T = zeros(M, 1);
mxdv_penwav = zeros(M, 1);
mxdv_scad = zeros(M, 1);
mxdv_mixed = zeros(M, 1);
mxdv_decompsh = zeros(M, 1);
mxdv_thrda1 = zeros(M, 1);
mxdv_bams = zeros(M, 1);
mxdv_blmed_A = zeros(M, 1);
mxdv_blmean_A = zeros(M, 1);
mxdv_blmed_T = zeros(M, 1);
mxdv_blmean_T = zeros(M, 1);
mxdv_singlmed = zeros(M, 1);
mxdv_singlmed2 = zeros(M, 1);
mxdv_singlmean = zeros(M, 1);
mxdv_singlmean2 = zeros(M, 1);
mxdv_singlhyp = zeros(M, 1);
mxdv_hybmed_A = zeros(M, 1);
mxdv_hybmean_A = zeros(M, 1);
mxdv_hybmed_T = zeros(M, 1);
mxdv_hybmean_T = zeros(M, 1);

% L1 Norm

L1_visu_H = zeros(M, 1);
L1_visu_S = zeros(M, 1);
L1_sure = zeros(M, 1);
L1_hybsure = zeros(M, 1);
L1_TI_H = zeros(M, 1);
L1_TI_S = zeros(M, 1);
L1_minimax_H = zeros(M, 1);
L1_minimax_S = zeros(M, 1);
L1_fdr_H = zeros(M, 1);
L1_fdr_S = zeros(M, 1);
L1_cv_H = zeros(M, 1);
L1_cv_S = zeros(M, 1);
L1_neighbl = zeros(M, 1);
L1_blockJS_A = zeros(M, 1);
L1_blockJS_T = zeros(M, 1);
L1_penwav = zeros(M, 1);
L1_scad = zeros(M, 1);
L1_mixed = zeros(M, 1);
L1_decompsh = zeros(M, 1);
L1_thrda1 = zeros(M, 1);
L1_bams = zeros(M, 1);
L1_blmed_A = zeros(M, 1);
L1_blmean_A = zeros(M, 1);
L1_blmed_T = zeros(M, 1);
L1_blmean_T = zeros(M, 1);
L1_singlmed = zeros(M, 1);
L1_singlmed2 = zeros(M, 1);
L1_singlmean = zeros(M, 1);
L1_singlmean2 = zeros(M, 1);
L1_singlhyp = zeros(M, 1);
L1_hybmed_A = zeros(M, 1);
L1_hybmean_A = zeros(M, 1);
L1_hybmed_T = zeros(M, 1);
L1_hybmean_T = zeros(M, 1);

% CPU Time

CPU_visu_H = zeros(M, 1);
CPU_visu_S = zeros(M, 1);
CPU_sure = zeros(M, 1);
CPU_hybsure = zeros(M, 1);
CPU_TI_H = zeros(M, 1);
CPU_TI_S = zeros(M, 1);
CPU_minimax_H = zeros(M, 1);
CPU_minimax_S = zeros(M, 1);
CPU_fdr_H = zeros(M, 1);
CPU_fdr_S = zeros(M, 1);
CPU_cv_H = zeros(M, 1);
CPU_cv_S = zeros(M, 1);
CPU_neighbl = zeros(M, 1);
CPU_blockJS_A = zeros(M, 1);
CPU_blockJS_T = zeros(M, 1);
CPU_penwav = zeros(M, 1);
CPU_scad = zeros(M, 1);
CPU_mixed = zeros(M, 1);
CPU_decompsh = zeros(M, 1);
CPU_thrda1 = zeros(M, 1);
CPU_bams = zeros(M, 1);
CPU_blmed_A = zeros(M, 1);
CPU_blmean_A = zeros(M, 1);
CPU_blmed_T = zeros(M, 1);
CPU_blmean_T = zeros(M, 1);
CPU_singlmed = zeros(M, 1);
CPU_singlmed2 = zeros(M, 1);
CPU_singlmean = zeros(M, 1);
CPU_singlmean2 = zeros(M, 1);
CPU_singlhyp = zeros(M, 1);
CPU_hybmed_A = zeros(M, 1);
CPU_hybmean_A = zeros(M, 1);
CPU_hybmed_T = zeros(M, 1);
CPU_hybmean_T = zeros(M, 1);

% Initializations for estimates

visu_H = zeros(1, n);
visu_S = zeros(1, n);
sure = zeros(1, n);
hybsure = zeros(1, n);
TI_H = zeros(1, n);
TI_S = zeros(1, n);
minimax_H = zeros(1, n);
minimax_S = zeros(1, n);
fdr_H = zeros(1, n);
fdr_S = zeros(1, n);
cv_H = zeros(1, n);
cv_S = zeros(1, n);
neighbl = zeros(1, n);
blockJS_A = zeros(1, n);
blockJS_T = zeros(1, n);
penwav = zeros(1, n);
scad = zeros(1, n);
mixed = zeros(1, n);
decompsh = zeros(1, n);
thrda1 = zeros(1, n);
bams = zeros(1, n);
blmed_A = zeros(1, n); 
blmean_A = zeros(1, n); 
blmed_T = zeros(1, n); 
blmean_T = zeros(1, n); 
singlmed = zeros(1, n);
singlmed2 = zeros(1, n);
singlmean = zeros(1, n);
singlmean2 = zeros(1, n);
singlhyp = zeros(1, n); 
hybmed_A = zeros(1, n);
hybmean_A = zeros(1, n);
hybmed_T = zeros(1, n);
hybmean_T = zeros(1, n);

all_visu_H = zeros(1, n);
all_visu_S = zeros(1, n);
all_sure = zeros(1, n);
all_hybsure = zeros(1, n);
all_TI_H = zeros(1, n);
all_TI_S = zeros(1, n);
all_minimax_H = zeros(1, n);
all_minimax_S = zeros(1, n);
all_fdr_H = zeros(1, n);
all_fdr_S = zeros(1, n);
all_cv_H = zeros(1, n);
all_cv_S = zeros(1, n);
all_neighbl = zeros(1, n);
all_blockJS_A = zeros(1, n);
all_blockJS_T = zeros(1, n);
all_penwav = zeros(1, n);
all_scad = zeros(1, n);
all_mixed = zeros(1, n);
all_decompsh = zeros(1, n);
all_thrda1 = zeros(1, n);
all_bams = zeros(1, n);
all_blmed_A = zeros(1, n);
all_blmean_A = zeros(1, n);
all_blmed_T = zeros(1, n);
all_blmean_T = zeros(1, n);
all_singlmed = zeros(1, n);
all_singlmed2 = zeros(1, n);
all_singlmean = zeros(1, n);
all_singlmean2 = zeros(1, n);
all_singlhyp = zeros(1, n);
all_hybmed_A = zeros(1, n);
all_hybmean_A = zeros(1, n);
all_hybmed_T = zeros(1, n);
all_hybmean_T = zeros(1, n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_visu_H2 = zeros(1, n);
all_visu_S2 = zeros(1, n);
all_sure2 = zeros(1, n);
all_hybsure2 = zeros(1, n);
all_TI_H2 = zeros(1, n);
all_TI_S2 = zeros(1, n);
all_minimax_H2 = zeros(1, n);
all_minimax_S2 = zeros(1, n);
all_fdr_H2 = zeros(1, n);
all_fdr_S2 = zeros(1, n);
all_cv_H2 = zeros(1, n);
all_cv_S2 = zeros(1, n);
all_neighbl2 = zeros(1, n);
all_blockJS2_A = zeros(1, n);
all_blockJS2_T = zeros(1, n);
all_penwav2 = zeros(1, n);
all_scad2 = zeros(1, n);
all_mixed2 = zeros(1, n);
all_decompsh2 = zeros(1, n);
all_thrda12 = zeros(1, n);
all_bams2 = zeros(1, n);
all_blmed2_A = zeros(1, n);
all_blmean2_A = zeros(1, n);
all_blmed2_T = zeros(1, n);
all_blmean2_T = zeros(1, n);
all_singlmed22 = zeros(1, n);
all_singlmed222 = zeros(1, n);
all_singlmean22 = zeros(1, n);
all_singlmean222 = zeros(1, n);
all_singlhyp2 = zeros(1, n);
all_hybmed2_A = zeros(1, n);
all_hybmean2_A = zeros(1, n);
all_hybmed2_T = zeros(1, n);
all_hybmean2_T = zeros(1, n);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reseting the generator of random numbers to its initial state

rand('state',0); 
randn('state',0);

% Generating the signal
signal = MakeSignalNewb(names(test), n);

% Starting the Simulations

nsignal=zeros(1,n);

for m=1:M
fprintf(' Evaluating simulation number ');
%
% Generate data
%
nsignal(1,:) = noisysignal(names(test), n, rsnr);

% Calculating the signal smooth estimates
 
 t=cputime;
 visu_H = recvisu(nsignal(1,:),'H',h);
 CPU_visu_H(m,1) = cputime -t;
 
 t=cputime;
 visu_S = recvisu(nsignal(1,:),'S',h);
 CPU_visu_S(m,1) = cputime -t;
 
 t=cputime;
 sure = recsure(nsignal(1,:),h);
 CPU_sure(m,1) = cputime -t;
 
 t=cputime;
 hybsure = rechybsure(nsignal(1,:),h);
 CPU_hybsure(m,1) = cputime -t;
 
 t=cputime;
 TI_H = recTI(nsignal(1,:),'H',h);
 CPU_TI_H(m,1) = cputime -t;
 
 t=cputime;
 TI_S = recTI(nsignal(1,:),'S',h);
 CPU_TI_S(m,1) = cputime -t;
 
 t=cputime;
 minimax_H = recminimax(nsignal(1,:),'H',h);
 CPU_minimax_H(m,1) = cputime -t;
 
 t=cputime;
 minimax_S = recminimax(nsignal(1,:),'S',h);
 CPU_minimax_S(m,1) = cputime -t;
 
 t=cputime;
 fdr_H = recfdr(nsignal(1,:),'H',h);
 CPU_fdr_H(m,1) = cputime -t;
 
 t=cputime;
 fdr_S = recfdr(nsignal(1,:),'S',h);
 CPU_fdr_S(m,1) = cputime -t;
 
 t=cputime;
 cv_H = reccv(nsignal(1,:),'H',h);
 CPU_cv_H(m,1) = cputime -t;
 
 t=cputime;
 cv_S = reccv(nsignal(1,:),'S',h);
 CPU_cv_S(m,1) = cputime -t;
 
 t=cputime;
 neighbl = recneighblock(nsignal(1,:),h);
 CPU_neighbl(m,1) = cputime -t;
 
 t=cputime;
 blockJS_A = recblockJS('Augment',nsignal(1,:),h);
 CPU_blockJS_A(m,1) = cputime -t;
 
 t=cputime;
 blockJS_T = recblockJS('Truncate',nsignal(1,:),h);
 CPU_blockJS_T(m,1) = cputime -t;
 
 t=cputime;
 penwav = recpenwav(nsignal(1,:),h);
 CPU_penwav(m,1) = cputime -t;
 
 t=cputime;
 scad = recscad(nsignal(1,:),h);
 CPU_scad(m,1) = cputime -t;
 
 t=cputime;
 mixed = recmixed(nsignal(1,:),h);
 CPU_mixed(m,1) = cputime -t;
 
 t=cputime;
 decompsh = recdecompsh(nsignal(1,:),h);
 CPU_decompsh(m,1) = cputime -t;
 
 t=cputime;
 thrda1 = recthreshda1(nsignal(1,:),h);
 CPU_thrda1(m,1) = cputime -t;
 
 t=cputime;
 bams = recbams(nsignal(1,:),h);
 CPU_bams(m,1) = cputime -t;
 
 t=cputime;
 blmed_A = recblockmed('Augment',nsignal(1,:),[],thet0b,h);
 CPU_blmed_A(m,1) = cputime -t;
 
 t=cputime;
 blmean_A = recblockmean('Augment',nsignal(1,:),[],thet0b,h);
 CPU_blmean_A(m,1) = cputime -t;
 
 t=cputime;
 blmed_T = recblockmed('Truncate',nsignal(1,:),[],thet0b,h);
 CPU_blmed_T(m,1) = cputime -t;
 
 t=cputime;
 blmean_T = recblockmean('Truncate',nsignal(1,:),[],thet0b,h);
 CPU_blmean_T(m,1) = cputime -t;
 
 t=cputime;
 singlmed = recsinglemed(nsignal(1,:),thet0s,h);
 CPU_singlmed(m,1) = cputime -t;
 
 t=cputime;
 singlmed2 = recsinglemed3(nsignal(1,:),thet0s2,h);
 CPU_singlmed2(m,1) = cputime -t;
 
 t=cputime;
 singlmean = recsinglemean(nsignal(1,:),thet0s,h);
 CPU_singlmean(m,1) = cputime -t;
 
 t=cputime;
 singlmean2 = recsinglemean3(nsignal(1,:),thet0s2,h);
 CPU_singlmean2(m,1) = cputime -t;
 
 t=cputime;
 singlhyp = recsinglehyp(nsignal(1,:),thet0s,h);
 CPU_singlhyp(m,1) = cputime -t;
 
 t=cputime;
 hybmed_A = rechybblockmed('Augment',nsignal(1,:),[],thet0s,thet0b,h);
 CPU_hybmed_A(m,1) = cputime -t;
 
 t=cputime;
 hybmean_A = rechybblockmean('Augment',nsignal(1,:),[],thet0s,thet0b,h);
 CPU_hybmean_A(m,1) = cputime -t;
 
 t=cputime;
 hybmed_T = rechybblockmed('Truncate',nsignal(1,:),[],thet0s,thet0b,h);
 CPU_hybmed_T(m,1) = cputime -t;
 
 t=cputime;
 hybmean_T = rechybblockmean('Truncate',nsignal(1,:),[],thet0s,thet0b,h);
 CPU_hybmean_T(m,1) = cputime -t;
 
fprintf('%4i%',m); % simulation number
fprintf('\n');

 % Calculating the Mean Squared Error
 mse_visu_H(m,1) = nanmean((visu_H-signal).^2);
 mse_visu_S(m,1) = nanmean((visu_S-signal).^2);
 mse_sure(m,1) = nanmean((sure-signal).^2);
 mse_hybsure(m,1) = nanmean((hybsure-signal).^2);
 mse_TI_H(m,1) = nanmean((TI_H-signal).^2);
 mse_TI_S(m,1) = nanmean((TI_S-signal).^2);
 mse_minimax_H(m,1) = nanmean((minimax_H-signal).^2);
 mse_minimax_S(m,1) = nanmean((minimax_S-signal).^2);
 mse_fdr_H(m,1) = nanmean((fdr_H-signal).^2);
 mse_fdr_S(m,1) = nanmean((fdr_S-signal).^2);
 mse_cv_H(m,1) = nanmean((cv_H-signal).^2);
 mse_cv_S(m,1) = nanmean((cv_S-signal).^2);
 mse_neighbl(m,1) = nanmean((neighbl-signal).^2);
 mse_blockJS_A(m,1) = nanmean((blockJS_A-signal).^2);
 mse_blockJS_T(m,1) = nanmean((blockJS_T-signal).^2);
 mse_penwav(m,1) = nanmean((penwav-signal).^2);
 mse_scad(m,1) = nanmean((scad-signal).^2);
 mse_mixed(m,1) = nanmean((mixed-signal).^2);
 mse_decompsh(m,1) = nanmean((decompsh-signal).^2);
 mse_thrda1(m,1) = nanmean((thrda1-signal).^2);
 mse_bams(m,1) = nanmean((bams-signal).^2);
 mse_blmed_A(m,1) = nanmean((blmed_A-signal).^2);
 mse_blmean_A(m,1) = nanmean((blmean_A-signal).^2);
 mse_blmed_T(m,1) = nanmean((blmed_T-signal).^2);
 mse_blmean_T(m,1) = nanmean((blmean_T-signal).^2);
 mse_singlmed(m,1) = nanmean((singlmed-signal).^2);
 mse_singlmed2(m,1) = nanmean((singlmed2-signal).^2);
 mse_singlmean(m,1) = nanmean((singlmean-signal).^2);
 mse_singlmean2(m,1) = nanmean((singlmean2-signal).^2);
 mse_singlhyp(m,1) = nanmean((singlhyp-signal).^2);
 mse_hybmed_A(m,1) = nanmean((hybmed_A-signal).^2);
 mse_hybmean_A(m,1) = nanmean((hybmean_A-signal).^2);
 mse_hybmed_T(m,1) = nanmean((hybmed_T-signal).^2);
 mse_hybmean_T(m,1) = nanmean((hybmean_T-signal).^2);
 
 % Stacking the various estimates
 
 all_visu_H = all_visu_H + visu_H;
 all_visu_S = all_visu_S + visu_S;
 all_sure = all_sure + sure;
 all_hybsure =  all_hybsure+ hybsure;
 all_TI_H = all_TI_H+ TI_H;
 all_TI_S = all_TI_S+TI_S;
 all_minimax_H = all_minimax_H+ minimax_H;
 all_minimax_S =  all_minimax_S+minimax_S;
 all_fdr_H =  all_fdr_H+ fdr_H;
 all_fdr_S = all_fdr_S+fdr_S;
 all_cv_H = all_cv_H+cv_H;
 all_cv_S = all_cv_S+cv_S;
 all_neighbl = all_neighbl+neighbl;
 all_blockJS_A =  all_blockJS_A+blockJS_A;
 all_blockJS_T =  all_blockJS_T+blockJS_T;
 all_penwav = all_penwav+penwav;
 all_scad =  all_scad+scad;
 all_mixed =  all_mixed+mixed;
 all_decompsh = all_decompsh+decompsh;
 all_thrda1 =  all_thrda1+thrda1;
 all_bams =  all_bams+bams;
 all_blmed_A = all_blmed_A+blmed_A;
 all_blmean_A = all_blmean_A+blmean_A;
 all_blmed_T = all_blmed_T+blmed_T;
 all_blmean_T = all_blmean_T+blmean_T;
 all_singlmed =  all_singlmed+singlmed;
 all_singlmed2 = all_singlmed2+singlmed2;
 all_singlmean = all_singlmean+singlmean;
 all_singlmean2 = all_singlmean2+singlmean2;
 all_singlhyp =  all_singlhyp+singlhyp;
 all_hybmed_A = all_hybmed_A+hybmed_A;
 all_hybmean_A =  all_hybmean_A+hybmean_A;
 all_hybmed_T = all_hybmed_T+hybmed_T;
 all_hybmean_T =  all_hybmean_T+hybmean_T;


 all_visu_H2 = all_visu_H2 + visu_H.^2;
 all_visu_S2 = all_visu_S2+ visu_S.^2;
 all_sure2 = all_sure2 + sure.^2;
 all_hybsure2 =  all_hybsure2+ hybsure.^2;
 all_TI_H2 = all_TI_H2 + TI_H.^2;
 all_TI_S2 = all_TI_S2 + TI_S.^2;
 all_minimax_H2 = all_minimax_H2 + minimax_H.^2;
 all_minimax_S2 =  all_minimax_S2+minimax_S.^2;
 all_fdr_H2 =  all_fdr_H2+ fdr_H.^2;
 all_fdr_S2 = all_fdr_S2+fdr_S.^2;
 all_cv_H2 = all_cv_H2+cv_H.^2;
 all_cv_S2 = all_cv_S2+cv_S.^2;
 all_neighbl2 = all_neighbl2+neighbl.^2;
 all_blockJS2_A =  all_blockJS2_A+blockJS_A.^2;
 all_blockJS2_T =  all_blockJS2_T+blockJS_T.^2;
 all_penwav2 = all_penwav2+penwav.^2;
 all_scad2 =  all_scad2+scad.^2;
 all_mixed2 =  all_mixed2+mixed.^2;
 all_decompsh2 = all_decompsh2+decompsh.^2;
 all_thrda12 =  all_thrda12+thrda1.^2;
 all_bams2 =  all_bams2+bams.^2;
 all_blmed2_A = all_blmed2_A+blmed_A.^2;
 all_blmean2_A = all_blmean2_A+blmean_A.^2;
 all_blmed2_T = all_blmed2_T+blmed_T.^2;
 all_blmean2_T = all_blmean2_T+blmean_T.^2;
 all_singlmed22 =  all_singlmed22+singlmed.^2;
 all_singlmed222 = all_singlmed222+singlmed2.^2;
 all_singlmean22 = all_singlmean2+singlmean.^2;
 all_singlmean222 = all_singlmean222+singlmean2.^2;
 all_singlhyp2 =  all_singlhyp2+singlhyp.^2;
 all_hybmed2_A = all_hybmed2_A+hybmed_A.^2;
 all_hybmean2_A =  all_hybmean2_A+hybmean_A.^2;
 all_hybmed2_T = all_hybmed2_T+hybmed_T.^2;
 all_hybmean2_T =  all_hybmean2_T+hybmean_T.^2;

 % Calculating the Maximun Deviation
 
 mxdv_visu_H(m,1) = nanmax(abs(visu_H-signal));
 mxdv_visu_S(m,1) = nanmax(abs((visu_S-signal)));
 mxdv_sure(m,1) = nanmax(abs((sure-signal)));
 mxdv_hybsure(m,1) = nanmax(abs((hybsure-signal)));
 mxdv_TI_H(m,1) = nanmax(abs((TI_H-signal)));
 mxdv_TI_S(m,1) = nanmax(abs((TI_S-signal)));
 mxdv_minimax_H(m,1) = nanmax(abs((minimax_H-signal)));
 mxdv_minimax_S(m,1) = nanmax(abs((minimax_S-signal)));
 mxdv_fdr_H(m,1) = nanmax(abs((fdr_H-signal)));
 mxdv_fdr_S(m,1) = nanmax(abs((fdr_S-signal)));
 mxdv_cv_H(m,1) = nanmax(abs((cv_H-signal)));
 mxdv_cv_S(m,1) = nanmax(abs((cv_S-signal)));
 mxdv_neighbl(m,1) = nanmax(abs((neighbl-signal)));
 mxdv_blockJS_A(m,1) = nanmax(abs((blockJS_A-signal)));
 mxdv_blockJS_T(m,1) = nanmax(abs((blockJS_T-signal)));
 mxdv_penwav(m,1) = nanmax(abs((penwav-signal)));
 mxdv_scad(m,1) = nanmax(abs((scad-signal)));
 mxdv_mixed(m,1) = nanmax(abs((mixed-signal)));
 mxdv_decompsh(m,1) = nanmax(abs((decompsh-signal)));
 mxdv_thrda1(m,1) = nanmax(abs((thrda1-signal)));
 mxdv_bams(m,1) = nanmax(abs((bams-signal)));
 mxdv_blmed_A(m,1) = nanmax(abs((blmed_A-signal)));
 mxdv_blmean_A(m,1) = nanmax(abs((blmean_A-signal)));
 mxdv_blmed_T(m,1) = nanmax(abs((blmed_T-signal)));
 mxdv_blmean_T(m,1) = nanmax(abs((blmean_T-signal)));
 mxdv_singlmed(m,1) = nanmax(abs((singlmed-signal)));
 mxdv_singlmed2(m,1) = nanmax(abs(singlmed2-signal));
 mxdv_singlmean(m,1) = nanmax(abs((singlmean-signal)));
 mxdv_singlmean2(m,1) = nanmax(abs((singlmean2-signal)));
 mxdv_singlhyp(m,1) = nanmax(abs((singlhyp-signal)));
 mxdv_hybmed_A(m,1) = nanmax(abs((hybmed_A-signal)));
 mxdv_hybmean_A(m,1) = nanmax(abs((hybmean_A-signal)));
 mxdv_hybmed_T(m,1) = nanmax(abs((hybmed_T-signal)));
 mxdv_hybmean_T(m,1) = nanmax(abs((hybmean_T-signal)));
 
 % Calculating the L1 Norm
 
 L1_visu_H(m,1) = nansum((abs(visu_H-signal)));
 L1_visu_S(m,1) = nansum((abs(visu_S-signal)));
 L1_sure(m,1) = nansum((abs(sure-signal)));
 L1_hybsure(m,1) = nansum(abs((hybsure-signal)));
 L1_TI_H(m,1) = nansum(abs((TI_H-signal)));
 L1_TI_S(m,1) = nansum(abs((TI_S-signal)));
 L1_minimax_H(m,1) = nansum(abs((minimax_H-signal)));
 L1_minimax_S(m,1) = nansum(abs((minimax_S-signal)));
 L1_fdr_H(m,1) = nansum(abs((fdr_H-signal)));
 L1_fdr_S(m,1) = nansum(abs((fdr_S-signal)));
 L1_cv_H(m,1) = nansum(abs((cv_H-signal)));
 L1_cv_S(m,1) = nansum(abs((cv_S-signal)));
 L1_neighbl(m,1) = nansum(abs((neighbl-signal)));
 L1_blockJS_A(m,1) = nansum(abs((blockJS_A-signal)));
 L1_blockJS_T(m,1) = nansum(abs((blockJS_T-signal)));
 L1_penwav(m,1) = nansum(abs((penwav-signal)));
 L1_scad(m,1) = nansum(abs((scad-signal)));
 L1_mixed(m,1) = nansum(abs((mixed-signal)));
 L1_decompsh(m,1) = nansum(abs((decompsh-signal)));
 L1_thrda1(m,1) = nansum(abs((thrda1-signal)));
 L1_bams(m,1) = nansum(abs((bams-signal)));
 L1_blmed_A(m,1) = nansum(abs((blmed_A-signal)));
 L1_blmean_A(m,1) = nansum(abs((blmean_A-signal)));
 L1_blmed_T(m,1) = nansum(abs((blmed_T-signal)));
 L1_blmean_T(m,1) = nansum(abs((blmean_T-signal)));
 L1_singlmed(m,1) = nansum(abs((singlmed-signal)));
 L1_singlmed2(m,1) = nansum(abs((singlmed2-signal)));
 L1_singlmean(m,1) = nansum(abs((singlmean-signal)));
 L1_singlmean2(m,1) = nansum(abs((singlmean2-signal)));
 L1_singlhyp(m,1) = nansum(abs((singlhyp-signal)));
 L1_hybmed_A(m,1) = nansum(abs((hybmed_A-signal)));
 L1_hybmean_A(m,1) = nansum(abs((hybmean_A-signal)));
 L1_hybmed_T(m,1) = nansum(abs((hybmed_T-signal)));
 L1_hybmean_T(m,1) = nansum(abs((hybmean_T-signal)));

end

% Calculating the Root Mean Squared Error

rmse_visu_H = sqrt(nanmean(mse_visu_H));
rmse_visu_S = sqrt(nanmean(mse_visu_S));
rmse_sure = sqrt(nanmean(mse_sure));
rmse_hybsure = sqrt(nanmean(mse_hybsure));
rmse_TI_H = sqrt(nanmean(mse_TI_H));
rmse_TI_S = sqrt(nanmean(mse_TI_S));
rmse_minimax_H = sqrt(nanmean(mse_minimax_H));
rmse_minimax_S = sqrt(nanmean(mse_minimax_S));
rmse_fdr_H = sqrt(nanmean(mse_fdr_H));
rmse_fdr_S = sqrt(nanmean(mse_fdr_S));
rmse_cv_H = sqrt(nanmean(mse_cv_H));
rmse_cv_S = sqrt(nanmean(mse_cv_S));
rmse_neighbl = sqrt(nanmean(mse_neighbl));
rmse_blockJS_A = sqrt(nanmean(mse_blockJS_A));
rmse_blockJS_T = sqrt(nanmean(mse_blockJS_T));
rmse_penwav = sqrt(nanmean(mse_penwav));
rmse_scad = sqrt(nanmean(mse_scad));
rmse_mixed = sqrt(nanmean(mse_mixed));
rmse_decompsh = sqrt(nanmean(mse_decompsh));
rmse_thrda1 = sqrt(nanmean(mse_thrda1));
rmse_bams = sqrt(nanmean(mse_bams));
rmse_blmed_A = sqrt(nanmean(mse_blmed_A));
rmse_blmean_A = sqrt(nanmean(mse_blmean_A));
rmse_blmed_T = sqrt(nanmean(mse_blmed_T));
rmse_blmean_T = sqrt(nanmean(mse_blmean_T));
rmse_singlmed = sqrt(nanmean(mse_singlmed));
rmse_singlmed2 = sqrt(nanmean(mse_singlmed2));
rmse_singlmean = sqrt(nanmean(mse_singlmean));
rmse_singlmean2 = sqrt(nanmean(mse_singlmean2));
rmse_singlhyp = sqrt(nanmean(mse_singlhyp));
rmse_hybmed_A = sqrt(nanmean(mse_hybmed_A));
rmse_hybmean_A = sqrt(nanmean(mse_hybmean_A));
rmse_hybmed_T = sqrt(nanmean(mse_hybmed_T));
rmse_hybmean_T = sqrt(nanmean(mse_hybmean_T));

% Calculating the Root Mean Squared Bias

m_visu_H = all_visu_H/M;
m_visu_S = all_visu_S/M;
m_sure = all_sure/M;
m_hybsure = all_hybsure/M;
m_TI_H = all_TI_H/M;
m_TI_S = all_TI_S/M;
m_minimax_H = all_minimax_H/M;
m_minimax_S = all_minimax_S/M;
m_fdr_H = all_fdr_H/M;
m_fdr_S = all_fdr_S/M;
m_cv_H = all_cv_H/M;
m_cv_S = all_cv_S/M;
m_neighbl = all_neighbl/M;
m_blockJS_A = all_blockJS_A/M;
m_blockJS_T = all_blockJS_T/M;
m_penwav = all_penwav/M;
m_scad = all_scad/M;
m_mixed = all_mixed/M;
m_decompsh = all_decompsh/M;
m_thrda1 = all_thrda1/M;
m_bams = all_bams/M;
m_blmed_A = all_blmed_A/M; 
m_blmean_A = all_blmean_A/M;
m_blmed_T = all_blmed_T/M; 
m_blmean_T = all_blmean_T/M; 
m_singlmed = all_singlmed/M;
m_singlmed2 = all_singlmed2/M;
m_singlmean = all_singlmean/M;
m_singlmean2 = all_singlmean2/M;
m_singlhyp = all_singlhyp/M; 
m_hybmed_A = all_hybmed_A/M;
m_hybmean_A = all_hybmean_A/M;
m_hybmed_T = all_hybmed_T/M;
m_hybmean_T = all_hybmean_T/M;
 
rmsb_visu_H = sqrt(nanmean((m_visu_H-signal).^2));
rmsb_visu_S = sqrt(nanmean((m_visu_S-signal).^2));
rmsb_sure = sqrt(nanmean((m_sure-signal).^2));
rmsb_hybsure = sqrt(nanmean((m_hybsure-signal).^2));
rmsb_TI_H = sqrt(nanmean((m_TI_H-signal).^2));
rmsb_TI_S = sqrt(nanmean((m_TI_S-signal).^2));
rmsb_minimax_H = sqrt(nanmean((m_minimax_H-signal).^2));
rmsb_minimax_S = sqrt(nanmean((m_minimax_S-signal).^2));
rmsb_fdr_H = sqrt(nanmean((m_fdr_H-signal).^2));
rmsb_fdr_S = sqrt(nanmean((m_fdr_S-signal).^2));
rmsb_cv_H = sqrt(nanmean((m_cv_H-signal).^2));
rmsb_cv_S = sqrt(nanmean((m_cv_S-signal).^2));
rmsb_neighbl = sqrt(nanmean((m_neighbl-signal).^2));
rmsb_blockJS_A = sqrt(nanmean((m_blockJS_A-signal).^2));
rmsb_blockJS_T = sqrt(nanmean((m_blockJS_T-signal).^2));
rmsb_penwav = sqrt(nanmean((m_penwav-signal).^2));
rmsb_scad = sqrt(nanmean((m_scad-signal).^2));
rmsb_mixed = sqrt(nanmean((m_mixed-signal).^2));
rmsb_decompsh = sqrt(nanmean((m_decompsh-signal).^2));
rmsb_thrda1 = sqrt(nanmean((m_thrda1-signal).^2));
rmsb_bams = sqrt(nanmean((m_bams-signal).^2));
rmsb_blmed_A = sqrt(nanmean((m_blmed_A-signal).^2));
rmsb_blmean_A = sqrt(nanmean((m_blmean_A-signal).^2));
rmsb_blmed_T = sqrt(nanmean((m_blmed_T-signal).^2));
rmsb_blmean_T = sqrt(nanmean((m_blmean_T-signal).^2));
rmsb_singlmed = sqrt(nanmean((m_singlmed-signal).^2));
rmsb_singlmed2 = sqrt(nanmean((m_singlmed2-signal).^2));
rmsb_singlmean = sqrt(nanmean((m_singlmean-signal).^2));
rmsb_singlmean2 = sqrt(nanmean((m_singlmean2-signal).^2));
rmsb_singlhyp = sqrt(nanmean((m_singlhyp-signal).^2));
rmsb_hybmed_A = sqrt(nanmean((m_hybmed_A-signal).^2));
rmsb_hybmean_A = sqrt(nanmean((m_hybmean_A-signal).^2));
rmsb_hybmed_T = sqrt(nanmean((m_hybmed_T-signal).^2));
rmsb_hybmean_T = sqrt(nanmean((m_hybmean_T-signal).^2));

% Calculating the Average width of the 2 SD band

std_visu_H = sqrt(all_visu_H2/M - m_visu_H.^2);
std_visu_S = sqrt(all_visu_S2/M - m_visu_S.^2);
std_sure = sqrt(all_sure2/M - m_sure.^2);
std_hybsure = sqrt(all_hybsure2/M - m_hybsure.^2);
std_TI_H = sqrt(all_TI_H2/M - m_TI_H.^2);
std_TI_S = sqrt(all_TI_S2/M - m_TI_S.^2);
std_minimax_H = sqrt(all_minimax_H2/M-m_minimax_H.^2);
std_minimax_S = sqrt(all_minimax_S2/M - m_minimax_S.^2);
std_fdr_H = sqrt(all_fdr_H2/M- m_fdr_H.^2);
std_fdr_S = sqrt(all_fdr_S2/M - m_fdr_S.^2);
std_cv_H = sqrt(all_cv_H2/M-m_cv_H.^2);
std_cv_S = sqrt(all_cv_S2/M - m_cv_S.^2);
std_neighbl = sqrt(all_neighbl2/M - m_neighbl.^2);
std_blockJS_A = sqrt(all_blockJS2_A/M - m_blockJS_A.^2);
std_blockJS_T = sqrt(all_blockJS2_T/M - m_blockJS_T.^2);
std_penwav = sqrt(all_penwav2/M - m_penwav.^2);
std_scad = sqrt(all_scad2/M - m_scad.^2);
std_mixed = sqrt(all_mixed2/M - m_mixed.^2);
std_decompsh = sqrt(all_decompsh2/M - m_decompsh.^2);
std_thrda1 = sqrt(all_thrda12/M - m_thrda1.^2);
std_bams = sqrt(all_bams2/M - m_bams.^2);
std_blmed_A = sqrt(all_blmed2_A/M - m_blmed_A.^2); 
std_blmean_A = sqrt(all_blmean2_A/M - m_blmean_A.^2);
std_blmed_T = sqrt(all_blmed2_T/M - m_blmed_T.^2); 
std_blmean_T = sqrt(all_blmean2_T/M - m_blmean_T.^2); 
std_singlmed = sqrt(all_singlmed22/M - m_singlmed.^2);
std_singlmed2 = sqrt(all_singlmed222/M - m_singlmed2.^2); 
std_singlmean = sqrt(all_singlmean22/M -m_singlmean.^2);
std_singlmean2 = sqrt(all_singlmean222/M - m_singlmean2.^2);
std_singlhyp = sqrt(all_singlhyp2/M - m_singlhyp.^2); 
std_hybmed_A = sqrt(all_hybmed2_A/M - m_hybmed_A.^2);
std_hybmean_A = sqrt(all_hybmean2_A/M - m_hybmean_A.^2);
std_hybmed_T = sqrt(all_hybmed2_T/M - m_hybmed_T.^2);
std_hybmean_T = sqrt(all_hybmean2_T/M - m_hybmean_T.^2);
 
band_visu_H = (2/n)*nansum(std_visu_H);
band_visu_S = (2/n)*nansum(std_visu_S);
band_sure = (2/n)*nansum(std_sure);
band_hybsure = (2/n)*nansum(std_hybsure);
band_TI_H = (2/n)*nansum(std_TI_H);
band_TI_S = (2/n)*nansum(std_TI_S);
band_minimax_H = (2/n)*nansum(std_minimax_H);
band_minimax_S = (2/n)*nansum(std_minimax_S);
band_fdr_H = (2/n)*nansum(std_fdr_H);
band_fdr_S = (2/n)*nansum(std_fdr_S);
band_cv_H = (2/n)*nansum(std_cv_H);
band_cv_S = (2/n)*nansum(std_cv_S);
band_neighbl = (2/n)*nansum(std_neighbl);
band_blockJS_A = (2/n)*nansum(std_blockJS_A);
band_blockJS_T = (2/n)*nansum(std_blockJS_T);
band_penwav = (2/n)*nansum(std_penwav);
band_scad = (2/n)*nansum(std_scad);
band_mixed = (2/n)*nansum(std_mixed);
band_decompsh = (2/n)*nansum(std_decompsh);
band_thrda1 = (2/n)*nansum(std_thrda1);
band_bams = (2/n)*nansum(std_bams);
band_blmed_A = (2/n)*nansum(std_blmed_A);
band_blmean_A = (2/n)*nansum(std_blmean_A);
band_blmed_T = (2/n)*nansum(std_blmed_T);
band_blmean_T = (2/n)*nansum(std_blmean_T);
band_singlmed = (2/n)*nansum(std_singlmed);
band_singlmed2 = (2/n)*nansum(std_singlmed2);
band_singlmean = (2/n)*nansum(std_singlmean);
band_singlmean2 = (2/n)*nansum(std_singlmean2);
band_singlhyp = (2/n)*nansum(std_singlhyp);
band_hybmed_A = (2/n)*nansum(std_hybmed_A);
band_hybmean_A = (2/n)*nansum(std_hybmean_A);
band_hybmed_T = (2/n)*nansum(std_hybmed_T);
band_hybmean_T = (2/n)*nansum(std_hybmean_T);

%Stacking the Mean Squared Error Vectors
%Adding : mse_blhyp
mse=[mse_visu_H mse_visu_S mse_sure mse_hybsure mse_TI_H mse_TI_S mse_minimax_H mse_minimax_S mse_cv_H mse_cv_S mse_neighbl mse_blockJS_A mse_blockJS_T mse_thrda1 mse_fdr_H mse_fdr_S mse_penwav mse_scad mse_mixed mse_decompsh mse_blmed_A mse_blmed_T mse_hybmed_A mse_hybmed_T mse_blmean_A mse_blmean_T mse_hybmean_A mse_hybmean_T mse_singlmed mse_singlmed2 mse_singlmean  mse_singlmean2 mse_singlhyp mse_bams]; %

%Stacking the Maximum Deviation Vectors
%Adding : mxdv_blhyp
mxdv=[mxdv_visu_H mxdv_visu_S mxdv_sure mxdv_hybsure mxdv_TI_H mxdv_TI_S mxdv_minimax_H mxdv_minimax_S mxdv_cv_H mxdv_cv_S mxdv_neighbl mxdv_blockJS_A mxdv_blockJS_T mxdv_thrda1 mxdv_fdr_H mxdv_fdr_S mxdv_penwav mxdv_scad mxdv_mixed mxdv_decompsh mxdv_blmed_A mxdv_blmed_T mxdv_hybmed_A mxdv_hybmed_T mxdv_blmean_A mxdv_blmean_T mxdv_hybmean_A mxdv_hybmean_T mxdv_singlmed mxdv_singlmed2 mxdv_singlmean mxdv_singlmean2 mxdv_singlhyp mxdv_bams]; %
	
%Stacking the L1 Norm Vectors
%Adding : L1_blhyp
L1norm=[L1_visu_H L1_visu_S L1_sure L1_hybsure L1_TI_H L1_TI_S L1_minimax_H L1_minimax_S L1_cv_H L1_cv_S L1_neighbl L1_blockJS_A L1_blockJS_T L1_thrda1 L1_fdr_H L1_fdr_S L1_penwav L1_scad L1_mixed L1_decompsh L1_blmed_A L1_blmed_T L1_hybmed_A L1_hybmed_T L1_blmean_A L1_blmean_T L1_hybmean_A L1_hybmean_T L1_singlmed L1_singlmed2 L1_singlmean L1_singlmean2 L1_singlhyp L1_bams]; %

%Stacking the CPU time Vectors
%Adding : CPU_blhyp
CPU=[CPU_visu_H CPU_visu_S CPU_sure CPU_hybsure CPU_TI_H CPU_TI_S CPU_minimax_H CPU_minimax_S CPU_cv_H CPU_cv_S CPU_neighbl CPU_blockJS_A CPU_blockJS_T CPU_thrda1 CPU_fdr_H CPU_fdr_S CPU_penwav CPU_scad CPU_mixed CPU_decompsh CPU_blmed_A CPU_blmed_T CPU_hybmed_A CPU_hybmed_T CPU_blmean_A CPU_blmean_T CPU_hybmean_A CPU_hybmean_T CPU_singlmed CPU_singlmed2 CPU_singlmean CPU_singlmean2 CPU_singlhyp CPU_bams]; %

% amse_f[:,1] : Average Mean Squared Error
% amse_f[:,2] : Standard Deviation
amse_f = [nanmean(mse)' nanstd(mse)'];

% rmse_f : Root Mean Squarred Error
%Adding : rmse_blhyp
rmse_f = [rmse_visu_H rmse_visu_S rmse_sure rmse_hybsure rmse_TI_H rmse_TI_S rmse_minimax_H rmse_minimax_S rmse_cv_H rmse_cv_S rmse_neighbl rmse_blockJS_A rmse_blockJS_T rmse_thrda1 rmse_fdr_H rmse_fdr_S rmse_penwav rmse_scad rmse_mixed rmse_decompsh rmse_blmed_A rmse_blmed_T rmse_hybmed_A rmse_hybmed_T rmse_blmean_A rmse_blmean_T rmse_hybmean_A rmse_hybmean_T rmse_singlmed rmse_singlmed2 rmse_singlmean rmse_singlmean2 rmse_singlhyp rmse_bams]; %

% rmsb_f : Root Mean Squarred Bias
%Adding : rmsb_blhyp
rmsb_f = [rmsb_visu_H rmsb_visu_S rmsb_sure rmsb_hybsure rmsb_TI_H rmsb_TI_S rmsb_minimax_H rmsb_minimax_S rmsb_cv_H rmsb_cv_S rmsb_neighbl rmsb_blockJS_A rmsb_blockJS_T rmsb_thrda1 rmsb_fdr_H rmsb_fdr_S rmsb_penwav rmsb_scad rmsb_mixed rmsb_decompsh rmsb_blmed_A rmsb_blmed_T rmsb_hybmed_A rmsb_hybmed_T rmsb_blmean_A rmsb_blmean_T rmsb_hybmean_A rmsb_hybmean_T rmsb_singlmed rmsb_singlmed2 rmsb_singlmean rmsb_singlmean2 rmsb_singlhyp rmsb_bams]; %

% amxdv_f[:,1] : Average Maximum Deviation
% amxdv_f[:,2] : Standard Deviation
amxdv_f = [nanmean(mxdv)' nanstd(mxdv)'];

% band_f : Average width of the 2 SD band
%Adding : band_blhyp
band_f = [band_visu_H band_visu_S band_sure band_hybsure band_TI_H band_TI_S band_minimax_H band_minimax_S band_cv_H band_cv_S band_neighbl band_blockJS_A band_blockJS_T band_thrda1 band_fdr_H band_fdr_S band_penwav band_scad band_mixed band_decompsh band_blmed_A band_blmed_T band_hybmed_A band_hybmed_T band_blmean_A band_blmean_T band_hybmean_A band_hybmean_T band_singlmed band_singlmed2 band_singlmean band_singlmean2 band_singlhyp band_bams]; %

% aL1norm_f[:,1] : Average L1 Norm
% aL1norm_f[:,2] : Standard Deviation
aL1norm_f = [nanmean(L1norm)' nanstd(L1norm)'];

% aCPUtime_f[:,1] : Average CPU time
% aCPUtime_f[:,2] : Standard Deviation
aCPUtime_f = [nanmean(CPU)' nanstd(CPU)'];
fignum=0;

if show == 1,
% Producing Box-Plots
figure(1)
clf;
boxplot(mse);
ylabel('Mean Squared Error');
title([char(names(test)) ' : Various Smooth Estimates']);
gname=sprintf('t%in%ih%ir%iMSE',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(2)
clf;
boxplot(mxdv);
ylabel('Maximum Deviation');
title([char(names(test)) ' : Various Smooth Estimates'])
gname=sprintf('t%in%ih%ir%iMXDV',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(3)
clf;
boxplot(L1norm);
ylabel('L1 Norm');
title([char(names(test)) ' : Various Smooth Estimates']);
gname=sprintf('t%in%ih%ir%iL1',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(4)
clf;
bar(rmse_f);
ylabel('Root Mean Squared Error');
title([char(names(test)) ' : Various Smooth Estimates']);
gname=sprintf('t%in%ih%ir%iRMSE',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(5)
clf;
bar(rmsb_f);
ylabel('Root Mean Squared Bias');
title([char(names(test)) ' : Various Smooth Estimates']);
gname=sprintf('t%in%ih%ir%iRMSB',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(6)
clf;
bar(band_f);
ylabel('Average width of the 2 SD band');
title([char(names(test)) ' : Various Smooth Estimates']);
gname=sprintf('t%in%ih%ir%iBAND',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(7)
clf;
boxplot(CPU);
ylabel('CPU Time');
title([char(names(test)) ' : Various Smooth Estimates']);
gname=sprintf('t%in%ih%ir%iCPU',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(8)
clf;
xt=(1:n)/n;
subplot(2,2,1);
plot(xt,m_visu_H,xt,signal);
ylabel('Visu-H');
subplot(2,2,2);
plot(xt,m_visu_S,xt,signal);
ylabel('Visu-S');
subplot(2,2,3);
plot(xt,m_sure,xt,signal);
ylabel('sure');
subplot(2,2,4);
plot(xt,m_hybsure,xt,signal);
ylabel('hybsure');
gname=sprintf('t%in%ih%ir%iAVE1',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(9)
clf;
subplot(2,2,1);
plot(xt,m_TI_H,xt,signal);
ylabel('TI-H');
subplot(2,2,2);
plot(xt,m_TI_S,xt,signal);
ylabel('TI-S');
subplot(2,2,3);
plot(xt,m_minimax_H,xt,signal);
ylabel('minimax-H');
subplot(2,2,4);
plot(xt,m_minimax_S,xt,signal);
ylabel('minimax-S');
gname=sprintf('t%in%ih%ir%iAVE2',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(10)
clf;
subplot(2,2,1);
plot(xt,m_cv_H,xt,signal);
ylabel('cv-H');
subplot(2,2,2);
plot(xt,m_cv_S,xt,signal);
ylabel('cv-S');
subplot(2,2,3);
plot(xt,m_neighbl,xt,signal);
ylabel('neighbl');
subplot(2,2,4);
plot(xt,m_blockJS_A,xt,signal);
ylabel('blockJS-A');
gname=sprintf('t%in%ih%ir%iAVE3',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(11)
clf;
subplot(2,2,1);
plot(xt,m_blockJS_T,xt,signal);
ylabel('blockJS-T');
subplot(2,2,2);
plot(xt,m_thrda1,xt,signal);
ylabel('thrda1');
subplot(2,2,3);
plot(xt,m_penwav,xt,signal);
ylabel('penwav');
subplot(2,2,4);
plot(xt,m_scad,xt,signal);
ylabel('scad');
gname=sprintf('t%in%ih%ir%iAVE4',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(12)
clf;
subplot(2,2,1);
plot(xt,m_mixed,xt,signal);
ylabel('mixed');
subplot(2,2,2);
plot(xt,m_decompsh,xt,signal);
ylabel('decompsh');
subplot(2,2,3);
plot(xt,m_blmed_A,xt,signal);
ylabel('blmed-A');
subplot(2,2,4);
plot(xt,m_blmed_T,xt,signal);
ylabel('blmed-T');
gname=sprintf('t%in%ih%ir%iAVE5',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(13)
clf;
subplot(2,2,1);
plot(xt,m_hybmed_A,xt,signal);
ylabel('hybmed-A');
subplot(2,2,2);
plot(xt,m_hybmed_T,xt,signal);
ylabel('hybmed-T');
subplot(2,2,3);
plot(xt,m_blmean_A,xt,signal);
ylabel('blmean-A');
subplot(2,2,4);
plot(xt,m_blmean_T,xt,signal);
ylabel('blmean-T');
gname=sprintf('t%in%ih%ir%iAVE6',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)


figure(14)
clf;
subplot(2,2,1);
plot(xt,m_hybmean_A,xt,signal);
ylabel('hybmean-A');
subplot(2,2,2);
plot(xt,m_hybmean_T,xt,signal);
ylabel('hybmean-T');
subplot(2,2,3);
plot(xt,m_singlmed,xt,signal);
ylabel('singlmed');
subplot(2,2,4);
plot(xt,m_singlmed2,xt,signal);
ylabel('singlmed2');
gname=sprintf('t%in%ih%ir%iAVE7',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(15)
clf;
subplot(2,2,1);
plot(xt,m_singlmean,xt,signal);
ylabel('singlmean');
subplot(2,2,2);
plot(xt,m_singlmean2,xt,signal);
ylabel('singlmean2');
subplot(2,2,3);
plot(xt,m_singlhyp,xt,signal);
ylabel('singlhyp');
subplot(2,2,4);
plot(xt,m_bams,xt,signal);
ylabel('bams');
gname=sprintf('t%in%ih%ir%iAVE7',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

figure(16)
clf;
subplot(2,1,1);
plot(xt,m_fdr_H,xt,signal);
ylabel('fdr-H');
subplot(2,1,2);
plot(xt,m_fdr_S,xt,signal);
ylabel('fdr-S');
gname=sprintf('t%in%ih%ir%iAVE8',test,n,hint,rsnr);
fignum=fignum+1;
%laprint(fignum,gname, 'keepticklabels')
printFigure(2,gname)

end % end 'if show == 1'


fprintf('\n')
info.rnames = strvcat('Rows','visu_H','visu_S','sure', ...
 'hybsure','TI_H', 'TI_S','minimax_H', 'minimax_S', 'cv_H', 'cv_S', 'neighbl', ...
 'blockJS_A', 'blockJS_T', 'thrda1','fdr_H', 'fdr_S', 'penwav', 'scad', 'mixed', 'decompsh', 'blmed_A','blmed_T','hybmed_A','hybmed_T', 'blmean_A','blmean_T', 'hybmean_A','hybmean_T', ...
'singlmed','singlmed2','singlmean','singlmean2','singlhyp','bams'); % 

resultats=[amse_f amxdv_f aL1norm_f aCPUtime_f rmse_f' rmsb_f' band_f'];
[info.endr,info.endc]=size(resultats);
info.cnames = strvcat('Average MSE','St Dev','Average MXDV','St Dev','Average L1norm','St Dev', ...
'Average CPU','St Dev','RMSE','RMSB','BAND');
info.fid=fidname;
lprint(resultats,info);
fclose(fidname);
f = struct('amse_f', amse_f, ...
		   'rmse_f', rmse_f, ...
         'rmsb_f', rmsb_f, ...
		   'amxdv_f', amxdv_f, ...
		   'band_f', band_f, ...
		   'aL1norm_f', aL1norm_f, ...
         'aCPUtime_f', aCPUtime_f) ; 
      
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
