% This is an example of how to apply the SIMULATIONS procedure.

M=100; show=1;

for test=1:12,
		for n=[256 512 1024],
			for rsnr=[3 5 7],
				for hint=1:2,
f = simulations(test,M,n,rsnr,[0.5 1],[0.5 1 0.7],...
                [0.5 1 0.7],hint,show);
			end
		end
	end
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
