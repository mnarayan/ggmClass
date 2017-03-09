function [Theta results] = constrainMLE(SigmaHat,InitialTheta,varargin);
	% Relaxed or Constrained or Re-fitted MLE
	% Maybe make this a twostageMLE of any kind?
	% Inputs:
	% Sigma - Sample Covariance or Correlation
	% InitialTheta - Initial Theta estimate by any method;
	
	
	p = size(InitialTheta,1);
	nlambdas = size(InitialTheta,3); 
	
	if(nargin==3)
		Lrange = varargin{1};
	else
		Lrange = [];
	end
	
	getmaxlambda = @(x)(max(max(triu(abs(x)))));
	lambda_max = getmaxlambda(SigmaHat); 
	
	
	function Lw = getLambdaWeight(Theta,method)
		
		Lw  = [];
		
		switch method 
			
		case 'binary'
			Lw = ones(p,p); 
		case 'inv'	
			Lw = 1./abs(Theta);		Lw = Lw/trace(Lw);
		case 'invsq'
			Lw = 1./abs(Theta).^2; 	Lw = Lw/trace(Lw);
			
		end
		
		nonzero_supp =  1.*(abs(Theta)>1e-5); 
		Lw(abs(nonzero_supp)==0) = 100*lambda_max; Lw(find(eye(p))) = 0;
	end
	
	
	if(nlambdas==1)
		
		% Adaptive
		Lambda = Lrange(1);
		Lw = getLambdaWeight(InitialTheta,'inv'); 
		[Theta SigmaHat2] = QUIC('default', SigmaHat, Lambda*Lw, 1e-6, 0, 1000);
	
		results.Lw = Lw;
		results.Sigma = SigmaHat2;
	
	elseif(nlambdas>1)

		% Constrained
		Lw = getLambdaWeight(InitialTheta,'inv'); 
		[Theta SigmaHat2] = QUIC('path', SigmaHat, Lw, Lrange, 1e-6, 0, 1000);
	
		results.Lw = Lw;
		results.Sigma = SigmaHat2;
	
	end
	
	
end