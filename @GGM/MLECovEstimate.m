function [self results] = MLECovEstimate(self,varargin)
%MLECOVESTIMATE	
	
	if(ndims(self.Data)==3)
		[m p n] = size(self.Data); 
	elseif(ndims(self.Data)==2)
		[m p] = size(self.Data); 
		n = 1;
	end

	% Default options
	options.verbose = false;
	options.weight = eye(p);
	options.standardize = false;

	if(options.standardize)
		[X sn_results] = successive_normalize(self.Data,options); 
		results.succnorm = sn_results;
		results.Xclean = X;		
	else
		X = self.Data;
		results.succnorm = [];
		results.Xclean = [];		
	end
	
	Sighat = sampleWeightedCov(X,options.weight,options);
	[D R] = varcorr(Sighat);
	
	self.Sigma = Sighat;
	results.D = D; 
	results.R = R; 

	function [D varargout] = varcorr(Sigma)
	%varcorr - Variance Correlation Decomposition
	% Sigma = diag(D) * R * diag(D) 	
			
		D = diag(Sigma); 
		if(nargout>=2)
			Dinv = 1./sqrt(D); 
			R = diag(Dinv)*Sigma*diag(Dinv);
			varargout{1} = R;
		end
		
	end

	function Sighat = sampleWeightedCov(X,M,options)		
		
		if(ndims(X)==3)
			[m p n] = size(X); 
		elseif(ndims(X)==2)
			[m p] = size(X); 
			n = 1;
		end
					
		Sighat = zeros(p,p);
		
		if(M(1,1)~=1)
			M = M/M(1,1); 
			if(options.verbose)
				sprintf('Scaling Weight Matrix by first entry')
			end
		end

		if n>1
			for cc=1:n
				Sighat = Sighat + X(:,:,cc)'*M*X(:,:,cc)/m;
			end
			Sighat = Sighat/n;
		else
			Sighat = (X'*X)/m;
		end
	end
	
	
end


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
	tol = 1e-5;
	n_iter = 1000; 
	k = 0;

	oldXr = X';
	oldXc = X;


	function [X mu sig]  = standardize(X)
		%STANDARDIZE is a helper function for SUCCESSIVE_NORMALIZE
			 
		mu 	= mean(X); 
		X 	= bsxfun(@minus,X,mu); 
		sig = std(X); 
		assert(any(sig<1e-5)==0,'Check data matrix for near constant or 0 rows or colums'); 
		X 	= bsxfun(@rdivide,X,sig);  	

	end

	while ((stop_err > tol) && (k<=n_iter))
	
		if(k==0)
			[X_colpolish mu_c sig_c] = standardize(X); 
			%[X_rowpolish mu_r sig_r] = standardize(X_colpolish');
			
			[X_rowpolish mu_r sig_r] = standardize(X');	
			%[X_colpolish mu_c sig_c] = standardize(X_rowpolish'); 
			
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


% function SigmaHat = whitenCov(self, X, Omega)
% 	% Whitens the Column Covariance of X using Omega
% 	% Apply whitenCov to X' in order to whiten the Row Covariance
%
% 	[m p n] = size(X);
%
% 	if(self.verbose)
% 		size1 = size(X(:,:,1)')
% 		size2 = size(Omega)
% 	else
% 		size1 = size(X(:,:,1)') ;
% 		size2 = size(Omega);
% 	end
%
% 	try
% 		assertEqual(size1(2),size2(1),'Data Matrix incompatible with Row Precision');
% 	catch Exception
% 		disp('Using I since specified Whitening matrix compatible')
% 		Omega = eye(size1(2));
% 	end
%
%
% 	SigmaHat = zeros(p,p);
% 	% Check that Omega(1,1) = 1;
% 	Omega = Omega/Omega(1,1);
% 	for cc=1:n
% 		SigmaHat = SigmaHat + X(:,:,cc)'*Omega*X(:,:,cc)/m;
% 	end
% 	SigmaHat = SigmaHat/n;
% 	if(self.useCorr)
% 		% SigmaHat = p*SigmaHat/trace(SigmaHat);
% 		What2 = sqrt(diag(diag(SigmaHat)));
% 		SigmaHat = What2^(-1)*SigmaHat*(What2^(-1));
% 	end
%
%
% end