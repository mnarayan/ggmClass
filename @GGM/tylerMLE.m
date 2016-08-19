function [self results] = tylerEstimate(self,varargin)
	% function [self results] = tylerEstimate(varargin)
	% Use robust covariance estimate instead of standard sample covariance. Based on Tyler's M-estimation for covariance matrix. 
	% 
	% Inputs
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
	eps_c = 1e-15;
	useShrinkage = 0;
	
	[m p n] = size(self.Data); 
	if(n>1)
		warning('Multiple Subjects not supported here'); 
	end
	
	SampleCov = cov(self.Data(:,:,1)); 
	
	
	if(m<p)
		useShrinkage = 0;
	end
	
	if(useShrinkage)
		trCov2 = trace(SampleCov)^2;
		tr2Cov = trace(SampleCov^2);
		rho = (tr2Cov + (1-2/p)*trCov2)/((1-m/p - 2*m/p^2)tr2Cov + (m+1+2*(m-1)/p)trCov2);
	else
		% if no shrinkage
		rho = 0;
	end
	
	initSigma = eye(p); 
	initTheta = eye(p);
	Sigma_k = eye(p); 
	Theta_k = eye(p);
	iter_score = zeros(2,max_iter); 
	
	X = self.Data(:,:,1);
	X = bsxfun(@rdivide, X', std(X',1))'; 
	assert(sum(abs(var(X',1))-1>5*eps)==0,'Samples not normalized');
	
	verbose=1;
	kk=1;
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
			% size(Ck)
			% size(X(ii,:)'*X(ii,:))
			Ck = Ck + (X(ii,:)'*X(ii,:))/(X(ii,:)*Theta_k*X(ii,:)');			
		end
		
		Sigma_k = (1-rho) * p/n * Ck + rho*eye(p); 
		Sigma_k = Sigma_k/(trace(Sigma_k)/p); 
		
		if(rho~=0)
			Theta_k = inv(Sigma_k); % no need for pinv as automatically regularized
		else
			Theta_k = pinv(Sigma_k); % as precau
		end
		
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
	
	results.kk = kk;
	results.rho = rho;
	results.eps_c = eps_c;
	results.iter_score = iter_score;
	results.converged = converged;
	
end
	