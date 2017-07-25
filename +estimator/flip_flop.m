function [results] = flip_flop(X,row_options,col_options,varargin)		
%FLIP_FLOP_ESTIMATORS - Fit kronecker-separable estimators
% USAGE: [results] = flip_flop(X,row_options,col_options
% 
% INPUTS
% - X is a m x p x n data matrix 
% - row_options, col_options are estimator.create_options for row and column covariance of X 

    ff_defaults = create_options();
    
    if(nargin>=4)
        iterative = varargin{1};
        ff_defaults.iterative = iterative;
    else
        iterative = ff_defaults.iterative;
    end
    
    if(~iterative)
        ff_defaults.max_iter = 1;
    end
    
    if(ndims==2)
        [m,p] = size(X);
    elseif(ndims==3)
        [m,p,n] = size(X);
    end
    
    max_iter = ff_defaults.max_iter;
    stop_err = ff_defaults.convergence.tol;
    k = 0;

	while ((stop_err > tol) && (k<=max_iter))
	
		if(k==0)
			
            row_M = eye(p);
            row_options.covariancefun = ...
                            @(rX)(covariance.mle_sample_covariance( ...
                                    rX,...
                                    struct('standardize','cols','M', row_M)...
                                    ));
            [rSig, rThet rS row_fit_opts] = estimate.fit(X',row_options);
			
            
            if(ndims(rThet)==2 & row_options.refit)
                col_M = rThet;
            else
                error('Row estimator needs to return single Theta');
            end
            
            col_options.covariancefun = ...
                            @(cX)(covariance.mle_sample_covariance( ...
                                    cX,...
                                    struct('standardize','cols','M', col_M)...
                                    ));
            [cSig, cThet cS col_fit_opts] = estimate.fit(X,col_options);
                
				
			oldSigc = cSig;
			oldSigr = rSig;
            oldThetc = cThet;
            oldThetr = rThet;
            
		else

            if(ndims(oldThetc)==2 & col_options.refit)
                row_M = oldThetc;
            else
                error('Col estimator needs to return single Theta');
            end
            row_options.covariancefun = ...
                            @(rX)(covariance.mle_sample_covariance( ...
                                    rX,...
                                    struct('standardize','cols','M', row_M)...
                                    ));
            [rSig, rThet rS row_fit_opts] = estimate.fit(X',row_options);
			
            
            if(ndims(oldThetr)==2 & row_options.refit)
                col_M = oldThetr;
            else
                error('Row estimator needs to return single Theta');
            end            
            col_options.covariancefun = ...
                            @(cX)(covariance.mle_sample_covariance( ...
                                    cX,...
                                    struct('standardize','cols','M', col_M)...
                                    ));
            [cSig, cThet cS col_fit_opts] = estimate.fit(X,col_options);

			
		end
	 
        row_err = (.5*norm(rSig-oldSigr,'fro'))/numel(rSig);
        col_err = (.5*norm(cSig-oldSigc,'fro'))/numel(cSig);
    
		stop_err = row_err + col_err;
	
		if(k==n_iter)
			if(options.verbose)				
				disp(sprintf(...
				'Iter:%d, stop_err:%.8f. Exiting loop.',k,stop_err));
			end
		end			
	
	
		k = k+1; 
		err(k) = stop_err;
		
		oldSigc = cSig;
		oldSigr = rSig;
        oldThetc = cThet;
        oldThetr = rThet;

	end


    results.row.options = row_options;
    results.row.covariance_estimate = rSig;
    results.row.inverse_covariance_estimate = rThet;
    results.row.sample_covariance = rS;
    results.row.fit_options = row_fit_opts;

    results.col.options = col_options;
    results.col.covariance_estimate = cSig;
    results.col.inverse_covariance_estimate = cThet;
    results.col.sample_covariance = cS;
    results.col.fit_options = col_fit_opts;

    results.ff_options = ff_defaults;
    results.ff_options.err = err;
    results.ff_options.row_err = row_err;
    results.ff_options.col_err = col_err;
    results.ff_options.k = k;
    
end


function ff_options = create_options()
   
    ff_options.iterative = false;
    ff_options.convergence.tol = 1e-5;
    ff_options.max_iter = 1;
    
end
