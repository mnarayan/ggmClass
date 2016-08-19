function [Theta results] = constrainMLE(Sigma,InitialTheta,varargin);
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
		L_w = zeros(p,p);
		L_w = 1./abs(InitialTheta); L_w = L_w/trace(L_w);
		L_w(InitialTheta==0) = 1000; L_w(find(eye(p))) = 0;

		[Theta SigmaHat2] = QUIC('default', SigmaHat, Lambda*L_w, 1e-6, 0, 1000);
	
		results.Lw = Lw;
		results.Sigma = SigmaHat2;
	
	elseif(nlambdas>1)
		
		% Constrained/Adaptive MLE 
		L_w = ones(p,p);
		L_w(self.InitialTheta==0) = 1000; L_w(find(eye(p))) = 0;

		[Theta SigmaHat3] = QUIC('path', SigmaHat, L_w, Lrange, 1e-6, 0, 1000);;
	
		results.Lw = Lw;
		results.Sigma = SigmaHat3;
	
	end
	
	
end