function [SigHat ThetaHat S options] = penalized_mle_inverse_covariance(X,varargin)		
%PENALIZED_MLE_INVERSE_COVARIANCE

    if(nargin==1)
        options = estimator.create_options();
    else
        options = varargin{1};
    end

    if(options.nlambdas>1 & options.Lambda & isscalar(options.Lambda))
        warning('Inconsistent options. Scalar Lambda takes precedence')
        disp('Setting options.nlambdas = 1'); 
        options.nlambdas = 1;
    end
    
    S = options.covariancefun(X);
    
    if(options.nlambdas>1)

        lrange = options.get_lambda(S);
        options.lambda_min = lrange{1};
        options.lambda_max = lrange{2};
        options.path = options.lambdafun(lrange{1},lrange{2});

    elseif(options.nlambdas==1)

        if(options.Lambda==0 | ...
            sum(sum(options.Lambda))==0 | ...
            isempty(options.Lambda) ...
            )
            warning('options.Lambda cannot be 0 if nlambdas = 1')
            lrange = options.get_lambda(S);
            options.lambda_min = lrange{1};
            options.lambda_max = lrange{2};
            options.Lambda = .1*lrange{2};            

        else

            lrange = options.get_lambda(S);
            options.lambda_min = lrange{1};
            options.lambda_max = lrange{2};
        
            if(options.Lambda>options.lambda_max | ...
                 options.Lambda < options.lambda_min)        
                warning('Lambda outside of allowed range')
                disp(['Lambda min' num2str(lrange{1})]);    
                disp(['Lambda max' num2str(lrange{2})]);    
                options.Lambda = .1*lrange{2};
                options.Lambda
            end
        end
    else        
       error('options.nlambdas cannot be zero') 
    end
 
    
    if(options.refit)
        
        switch options.init_estimator
        
        case 'lasso_mle'
            [SigHat0 ThetaHat0] = lasso_mle(X,options);
        otherwise
            warning('Only lasso_mle supported') 
        end
        
        [SigHat ThetaHat options] = constrain_mle(X,ThetaHat0,options);

    else

        switch options.estimator

        case 'lasso_mle'
           [SigHat ThetaHat] = lasso_mle(X,options);
        otherwise
           warning('Only lasso_mle supported') 
        end
    end

end


function [SigHat ThetaHat] = lasso_mle(X,options)
    
    S = options.covariancefun(X);
    
    switch options.solver
    case 'QUIC'
        [SigHat ThetaHat] = quic_helper(S,options);
    otherwise
        warning('Only QUIC solver currently supported'); 
    end
    
    % If in path mode, perform model selection
    % Get grphs
    
end

function [Sighat ThetaHat] = lasso_clime(X,options)
    
    S = options.covariancefun(X);
    
    
    
end


function [SigHat ThetaHat options] = constrain_mle(X,ThetaHat0,options)
    
    
    if(~isempty(options.W))
        warning('options.refit is turned on. Overwriting options.W');
    end 
    if(ndims(ThetaHat0)==2)       
        options.W = options.lambda_max*100.*(ThetaHat0==0) + ...
                            1.0.*(ThetaHat0~=0);
    else
       warning(['Need to perform model selection first. ...' ...
               ' Cannot refit to a path of estimates']);
    end

    switch options.estimator

    case 'lasso_mle'
       [SigHat ThetaHat] = lasso_mle(X,options);
    otherwise
       warning('Only lasso_mle supported') 
    end

    
end


function [SigHat ThetaHat] = chol_mle(X,options)
    
    

    
end


function [SigHat ThetaHat] = banded_mle(X,options)
    
    

    
end

function W = create_reg_matrix(Lambda,p)
   
   if(isscalar(Lambda))
       W =  ones(p,p).*Lambda;
       % make zero diagonals
       W(find(eye(p))) = 0;
   else   
       % symmetricize and set diagonals to zero
       if(ndims(Lambda)==2)
           W = triu(Lambda,1); 
           W = W + W';
       else
           warning('Lambda is not a matrix')
           Lambda
       end
   end
   
end


function [SigHat ThetaHat] = quic_helper(S,options)
%QUIC_HELPER

    p = size(S,1);
    if(isscalar(options.Lambda))
        if(options.Lambda==0 | isempty(options.Lambda))
            pathmode = true;
            options.Lambda = 1.0;
        else
            pathmode = false;
        end
    else
        if(isempty(options.Lambda))
            pathmode = true;
            options.Lambda = 1.0;
        end
    end
    
    if(isempty(options.W))
        if(pathmode)
            W = create_reg_matrix(1.0,p);
        else
            W = create_reg_matrix(options.Lambda,p);
        end
   else
       if(pathmode)
           W = create_reg_matrix(options.W,p); 
       else
           W = create_reg_matrix(options.Lambda,p); 
           W = W.*options.W; % Enforce constraints
       end
   end
    
    % Make sure diag(S) + (S.*L)_off is posdef.  
    function [chkposdef Lw] = isLwPSD(S,Lw)
        
        diagS = diag(diag(S)); 
        
        SLw = S-diagS;
        SLw = SLw.*Lw;        
        mineig = min(eig(diagS+SLw));       
        if(mineig<eps)
           warning('Regularization fails PSD check');
           chkposdef = false;
       else
           chkposdef = true;
       end
        
        if(abs(imag(mineig))>1e-03)
            warning('Imaginary eigenvalues for S.*W. Symmetrizing weights'); 
            Lw = (Lw + Lw')/2;
        end
    end
    
    if(pathmode)
        
        [chkpsd W] = isLwPSD(S,W);        
    	[ThetaHat SigHat] = QUIC('path', ...
    							S, ...
    							W, ...
    							options.path, ...
    							1e-6, ...
    							0, ...
    							1000 ...
    							);
    else
                
    	[ThetaHat SigHat] = QUIC('default', ...
    							S, ...
    							W, ...
    							1e-6, ...
    							0, ...
    							1000 ...
    							);
    end

end