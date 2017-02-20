function [Z_new X_new] = update_Z(X,Z,U,options)
	
	p 		= size(X,1);
	lambda	= options.lambda;
	rho 	= options.rho;
	alpha 	= options.alpha;
	if(~ismatrix(alpha))
		alpha = alpha.*ones(n,n);
		alpha = setdiagLS(alpha,1);
	end
	
    X_new = alpha.*X + (1 - alpha).*Z;
    Z_new = options.shrinkage(X_new + U, lambda./rho);	
	
end