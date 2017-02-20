function [output] = admm_solver(Sigma,varargin)
	
	narginchk(1,2)
	if nargin==1
		options = admm.admm_options();
	else
		options = varargin{1}; 
	end
	n = size(Sigma,1);	
	
	%% Global constants and defaults	
	QUIET 		= options.QUIET; 		
	MAX_ITER 	= options.MAX_ITER; 	
	ABSTOL		= options.ABSTOL;	
	RELTOL		= options.RELTOL;		
		
	%% Default options
	lambda 		= options.lambda;
	lambda		= sort(lambda,'ascend'); 
	rho			= options.rho; 
	alpha 		= options.alpha; 
	if(~ismatrix(alpha))
		alpha = alpha.*ones(n,n);
		alpha = setdiagLS(alpha,1);
	end
	
	%% Initialization
	options = admm.initializeXZU(Sigma,options); 
	X 			= options.X; 
	Z 			= options.Z; 
	U 			= options.U; 
	
	%% ADMM solver
	if ~QUIET
	    fprintf('%3s\t%10s\t%10s\t%10s\t%10s\t%10s\n', 'iter', ...
	      'r norm', 'eps pri', 's norm', 'eps dual', 'objective');
	end

	if(~options.pathwise)
		
		% check options.lambda
		if(length(lambda)>1)
			lambda = lambda(1); 
			options.lambda = lambda;
		end 
			
		t_start = cputime; 
		for k = 1:MAX_ITER
    
		    % z-update with relaxation
		    Zold = Z;
			[Z X_hat] = admm.update_Z(X,Zold,U,options);		
	
		    % x-update		
			X = admm.update_X(Sigma,X,Z,U,options);    

			% u-update
		    U = admm.update_U(X_hat,Z,U); 
		

		    % diagnostics, reporting, termination checks
		    history.objval(k)  = options.objective(Sigma, X, Z, lambda);
    
		    history.r_norm(k)  = norm(X - Z, 'fro');
		    history.s_norm(k)  = norm(-rho*(Z - Zold),'fro');
    
		    history.eps_pri(k) = sqrt(n*n)*ABSTOL + RELTOL*max(norm(X,'fro'), norm(Z,'fro'));
		    history.eps_dual(k)= sqrt(n*n)*ABSTOL + RELTOL*norm(rho*U,'fro');


		    if ~QUIET
		        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
		            history.r_norm(k), history.eps_pri(k), ...
		            history.s_norm(k), history.eps_dual(k), history.objval(k));
		    end

		    if (history.r_norm(k) < history.eps_pri(k) && ...
		       history.s_norm(k) < history.eps_dual(k))
		         break;
		    end
		end

		if ~QUIET
		    history.timing = cputime-t_start
		else
		    history.timing = cputime-t_start;
			if QUIET==1
				disp(sprintf('Convergence in %.4f secs',history.timing));
			end
		end
		
		output.Theta = X;
		output.history = history;
		output.options = options; 
		
	else
		
		output = {};
		
		for ll=1:length(options.lambda)
		
			% options warm starts
			options_ws = options;
			lambda_ll = options.lambda(ll); 			
			options_ws.lambda = lambda_ll;
			
			% Set X to previous output
			if(ll>1)
				X = output{ll-1}.Theta;
				Z = zeros(n,n); 
				U = zeros(n,n); 
			end
			
			t_start = cputime;		
			for k = 1:MAX_ITER    
			    % z-update with relaxation
			    Zold = Z;
				[Z X_hat] = admm.update_Z(X,Zold,U,options_ws);		
	
			    % x-update		
				X = admm.update_X(Sigma,X,Z,U,options_ws);    

				% u-update
			    U = admm.update_U(X_hat,Z,U); 		

			    % diagnostics, reporting, termination checks
			    history.objval(k)  = options.objective(Sigma, X, Z, lambda_ll);
    
			    history.r_norm(k)  = norm(X - Z, 'fro');
			    history.s_norm(k)  = norm(-rho*(Z - Zold),'fro');
    
			    history.eps_pri(k) = sqrt(n*n)*ABSTOL + RELTOL*max(norm(X,'fro'), norm(Z,'fro'));
			    history.eps_dual(k)= sqrt(n*n)*ABSTOL + RELTOL*norm(rho*U,'fro');

			    if ~QUIET
			        fprintf('%3d\t%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.2f\n', k, ...
			            history.r_norm(k), history.eps_pri(k), ...
			            history.s_norm(k), history.eps_dual(k), history.objval(k));
			    end

			    if (history.r_norm(k) < history.eps_pri(k) && ...
			       history.s_norm(k) < history.eps_dual(k))
			         break;
			    end
			end

			if ~QUIET
			    history.timing = cputime-t_start
			else
			    history.timing = cputime-t_start;
				if QUIET==1
					disp(sprintf('Convergence for lambda=%.2f in %.4f secs',lambda_ll,history.timing));
				end
			end				
		
			output{ll}.Theta = X;
			output{ll}.history = history;
			output{ll}.options = options_ws; 	
			
		end

	end
	
end



