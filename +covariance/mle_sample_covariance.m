function [SigHat results] = mle_sample_covariance(X,varargin)		
%MLE_SAMPLE_COVARIANCE
% 
	
	
	if(ndims(X)==3)
		[m p n] = size(X); 
	elseif(ndims(X)==2)
		[m p] = size(X); 
		n = 1;
	end
	
	narginchk(1,2); 
	if(nargin==2)
		options = varargin{1};
		if(~isfield(options,'M'))
			options.M = eye(m);
		end
		if(~isfield(options,'verbose'))
			options.verbose = true;
		end
		if(~isfield(options,'standardize'))
			options.standardize = 'rowcolumn';
		end
	else
		% default
		options.M = eye(m);
		options.verbose = true;
		options.standardize = 'rowcolumn';
	end	

	if(islogical(options.standardize))
        if(options.standardize)
            options.standardize = 'rowcolumn';
        else
            options.standardize = 'cols';
        end
    end
    
	% Initialize
	SigHat = zeros(p,p);
	M = options.M;
	
	if(M(1,1)~=1)
		M = M/M(1,1); 
		if(options.verbose)
			sprintf('Scaling Weight Matrix by first entry')
		end
	end

	if n>1
		for cc=1:n
			SigHat = SigHat + X(:,:,cc)'*M*X(:,:,cc)/m;
		end
        results.succnorm = [];
		SigHat = SigHat/n;
	else
        
    	switch options.standardize
        case 'rowcolumn'        
    		[Y succnorm] = standardize.successive_normalize(X');
        	results.succnorm = succnorm; 
            
    		X = Y';		
        case 'cols'
    		[X results.mu_c results.sig_c] = ...
                             standardize.standardize_cols(X);		
        otherwise
            warning('Using Sample Covariance');
            [X results.mu_c results.sig_c] = ...
                            standardize.standarize_cols(X,'center');
        end
        
		SigHat = (X'*M*X)/m;
	end
	
	results.X = X; 
	results.whitenMatrix = M; 	
	results.options = options;
	
end