function [X output] = successive_normalize(X,options)
	%SUCCESSIVE_NORMALIZE implements normalization strategy inspired 
	% by Olshen and Rajaratnam (2010); Allen and Tibhsirani (2011); 
	%
	% USAGE: [X output] = successive_normalize(X)
	% INPUT
	% 	- X : a n x p transposable data matrix
	% OUTPUT: 
	% 	- X : normalized n x p data matrix
	% 	- output : Contains original row,column parameters
	% 
	%  Copyright 2017, Manjari Narayan
	%  BSD-2 Clause License
	% 

	stop_err = Inf;
	tol = 1e-12;
	n_iter = 1000; 
	k = 0;


	function [X mu sig]  = standardize(X)
		%STANDARDIZE is a helper function for SUCCESSIVE_NORMALIZE
		
		import standardize.standardize_cols
		[X mu sig] = standardize_cols(X);

	end

	while ((stop_err > tol) && (k<=n_iter))
	
		if(k==0)
			[X_colpolish mu_c sig_c] = standardize(X); 
			[X_rowpolish mu_r sig_r] = standardize(X');	
			
			oldXc = X;
			oldXr = X';
		else
			[X_colpolish] = standardize(X_rowpolish'); 
			[X_rowpolish] = standardize(X_colpolish');					
		end

	
		stop_err = (.5*norm(X_colpolish-oldXc,'fro') + ...
			 		.5*norm(X_rowpolish-oldXr,'fro'))/numel(X);
	
		if(k==n_iter)
			if(options.verbose)				
				disp(sprintf(...
				'Iter:%d, stop_err:%.8f. Exiting loop.',k,stop_err));
			end
		end
	
		k = k+1; 
		oldXr = X_rowpolish;
		oldXc = X_colpolish;

	end

	X = (X_colpolish + X_rowpolish')/2;
	output = struct('mu_r',mu_r,'sig_r', sig_r,'mu_c', mu_c, 'sig_c',sig_c,'iter',k,'err',stop_err);

end