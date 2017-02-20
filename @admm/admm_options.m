function options = admm_options(varargin)
	
	options = {};
	
	% Global Constants
	options.QUIET 		= 1; 
	options.MAX_ITER 	= 1000; 
	options.ABSTOL		= 1e-4; 
	options.RELTOL		= 1e-3;
	
	% Default values
	options.lambda		= 0.1;
	options.rho			= 1; 
	options.alpha		= .5; 	
	
	% Function handles
	options.objective = @(Sigma,Theta,Z,lambda)(admm.objective(Sigma,Theta,Z,lambda)); 
	options.shrinkage = @(a,kappa)(shrinkage(a,kappa)); 
	
	% Mode 
	options.pathwise = false;
	
end


function y = shrinkage(a, kappa)
	if(ismatrix(a))
		a_diag = diag(a); 
	end
    y = max(0, a-kappa) - max(0, -a-kappa);
	if( ismatrix(a) && ismatrix(y) )
		y = setdiagLS(y,a_diag); 
	end
end