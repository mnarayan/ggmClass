function options = initializeXZU(Sigma,options)
	
	narginchk(2,2) 
	
	p = length(Sigma); 
	D = diag(Sigma);
	D(D==0) = 1;
	
	% Initialization of X,Z,U
	if(options.pathwise)
		lambda = sort(options.lambda,'ascend');
		lambda = lambda(1);
		% Initialize to dense Ridge estimate at smallest lambda
		if(exist('inv_posdef'))
			options.X = inv_posdef(Sigma.*(1-lambda)+lambda*eye(p));
		else
			options.X = pinv(Sighat.*(1-lambda)+lambda*eye(p));
		end
		% options.X = eye(p,p)
	else
		options.X = diag(1./D);
		% options.X = eye(p,p); 
	end		
	options.Z = options.X;
	options.U = zeros(p,p);

end