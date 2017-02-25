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
	
	if(nlambdas==1)
		Lambda = Lrange(1);
		Lw = zeros(p,p);
		Lw = 1./abs(InitialTheta); Lw = Lw/trace(Lw);
		Lw(InitialTheta==0) = 1000; Lw(find(eye(p))) = 0;

		[Theta SigmaHat2] = QUIC('default', SigmaHat, Lambda*Lw, 1e-6, 0, 1000);
	
		results.Lw = Lw;
		results.Sigma = SigmaHat2;
	
	elseif(nlambdas>1)
		
		% Constrained/Adaptive MLE 
		Lw = ones(p,p);
		Lw(self.InitialTheta==0) = 1000; Lw(find(eye(p))) = 0;

		[Theta SigmaHat3] = QUIC('path', SigmaHat, Lw, Lrange, 1e-6, 0, 1000);;
	
		results.Lw = Lw;
		results.Sigma = SigmaHat3;
	
	end
	
	
end