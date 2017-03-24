function [Sigma results] = tylermle_sample_covariance(X, varargin)
	%TYLERMLE_SAMPLE_COVARIANCE implements Tyler's M-estimator for a Covariance Matrix. This is a robust covariance estimate to replace the standard MLE for the sample covariance. Uses Steinian type shrinkage for high dimensional case.
	% 
	% USAGE: 
	% function [Sigma results] = tylerMLE(X,varargin)
	% 
	% INPUTS:
	% 	X 		- Data matrix of size m x p. Rows are observations
	% 	options -
	%   options.standardize = true (default), false
	% 	
	% 
	% S_i = X_i/||X_i||_2
	% 
	% Fixed of the following estimators 
	% Sighat_(k+1) = (1-rho) * (p/n * \sum_{i=1}^{n} (S_i' * S_i)/(S_i' Thethat_{k} S_i))
	% 						+ rho * I
	% Thethat_(k+1) = inv(Sighat_(k+1))
	% 
	% plugin estimator for rho
	% R - sample covariance using X 
	% 
	% rho = (Tr^2(R) + (1-2/p)Tr(R^2))/( (1-n/p - 2n/p^2)Tr^2(R) + (n+1+2(n-1)/p)Tr(R^2))
	
	max_iter = 100; 
	converged=0;
	eps_c = 1e-10;
	useShrinkage = 1;
	shrinkage_method = 'Chen-Hero'; % Chen-Hero, Wiesel Shrinkage, KL shrinkage
	
	if(ndims(X)==3)
		[m p n] = size(X); 
	elseif(ndims(X)==2)
		[m p] = size(X)
		n = 1;
	end
	
	if(n>1)
		warning('Multiple Subjects not supported here'); 
		X = X(:,:,1);		
	end
	
	narginchk(1,2);	
	if(nargin==2)
		options = varargin{1}; 
	else
		options = [];
	end
	
	if(isfield(options,'standardize'))
		if(options.standardize)
			[Y succnorm] = standardize.successive_normalize(X');
			X = Y';		
		else
			warning('This method will fail without standardizing data matrix')
			[X succnorm] = standardize.standardize_cols(X);					
			assert(sum(abs(var(X',1))-1>5*eps)~=0,...
										'Samples not normalized');				
		end
	else
		assert(sum(abs(var(X',1))-1>5*eps)~=0,'Samples not normalized');	
		fprintf('Standardizing rows and columns of X');		
		[Y succnorm] = standardize.successive_normalize(X');
		X = Y';				
	end
	results.options = options;
	results.succnorm = succnorm; 
	
	switch shrinkage_method
		
	case 'Chen-Hero'
		SampleCov = cov(X); 
		TargetCov = diag(diag(SampleCov));
		
	case 'Wiesel'
		SampleCov = cov(X); 
		TargetCov = diag(diag(SampleCov)); 		
	case 'KL'	
		SampleCov = cov(X); 
		TargetCov = diag(diag(SampleCov)); 		
	end
	
	if(m<10*p)
		useShrinkage = 1;
	end
	
	if(useShrinkage)
		trCov2 = trace(SampleCov)^2;
		tr2Cov = trace(SampleCov^2);
		rho = (tr2Cov + (1-2/p)*trCov2) /((1-m/p - 2*m/p^2)*tr2Cov + (m+1+2*(m-1)/p)*trCov2);
	else
		% if no shrinkage
		rho = 0;
	end
	
	initSigma = eye(p); 
	initTheta = eye(p);
	Sigma_k = eye(p); 
	Theta_k = eye(p);
	iter_score = zeros(2,max_iter); 	
	
	verbose=1;
	kk=1;
	Rk = zeros(1,m);  
	
	% for each iteration
	while((kk<max_iter)&&(converged==0))
		if(verbose)
			kk
		end
		% Store old iterates
		initSigma = Sigma_k;
		initTheta = Theta_k; 
		
		Ck = zeros(p,p);
		for ii=1:m
			Rk(ii) = (X(ii,:)*Theta_k*X(ii,:)');
			Ck = Ck + (X(ii,:)'*X(ii,:))/Rk(ii);
						
		end
		
		switch shrinkage_method
		case 'Chen-Hero'
			Sigma_k = (1-rho) * p/m * Ck + rho*TargetCov; % LW shrinkage
			Sigma_k = Sigma_k/trace(Sigma_k); 
		case 'Wiesel'
			% This can be reparametrized. Maybe faster to do this.
			% Sigma2 = sqrtm(inv(TargetCovf))*Sigma*sqrtm(inv(TargetCovf))
			% Would ensure that trace(Theta_k*TargetCov) = p; 
			Sigma_k = (1-rho) * p/m * Ck + ...
				 	rho*p*TargetCov/trace(Theta_k*TargetCov); 
			%Sigma_k = Sigma_k/(trace(Sigma_k)); 
		case 'KL'
			Sigma_k = (1-rho) * p/m * Ck + rho*TargetCov; 
			Sigma_k = Sigma_k/(trace(Sigma_k)); 			
		end
		
		if(rho~=0)
			Theta_k = inv(Sigma_k); % no need for pinv as automatically regularized
		else
			Theta_k = pinv(Sigma_k); % as precau
		end
		% if(strcmp(shrinkage_method,'Wiesel'))
		% 	Theta_k = Theta_k/(trace(Theta_k)/p);
		% end	
		
		% check convergence condition
		iter_score(1,kk) = sum(sum(abs(Sigma_k-initSigma).^2)); 
		iter_score(2,kk) = sum(sum(abs(Theta_k-initTheta).^2)); 
	
		if(verbose)
			disp(['Fixed point error']);
			iter_score(1,kk)
		end
	
		if(iter_score(1,kk)<eps_c)
			converged=1;
		end	
		
		kk=kk+1;
		
	end
	
	
	Sigma = Sigma_k/(trace(Sigma_k)/p);	
	tmpD = diag(1./sqrt(diag(Sigma)));
	Sigma = tmpD*Sigma*tmpD;
	
	results.Rk = Rk;
	results.W = TargetCov; 
	results.Sigma = Sigma;
	results.Theta = Theta_k;
	results.kk = kk;
	results.rho = rho;
	results.eps_c = eps_c;
	results.iter_score = iter_score;
	results.converged = converged;
	
end
	