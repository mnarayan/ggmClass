function [shrinkage] = stein_shrinkage_sample_covariance(Shat,varargin)		
%STEIN_SHRINKAGE_SAMPLE_COVARIANCE
% 
% References: 
% Ledoit-Wolf Linear Shrinkage
%   
% Lin-Perlman Stein Shrinkage
% 
% Condition Number Shrinkage
% 
% GPL-2 License
% 2017, Manjari Narayan
	
    % if(ndims(X)==3)
    %     [m p n] = size(X);
    % elseif(ndims(X)==2)
    %     [m p] = size(X);
    %     n = 1;
    % end
    
   
    if(nargin>=2)
        opts = varargin{1};
    else
        opts = estimator.create_options();
    end
    
    [V D] = eig(Shat); 
    eigcoefs = diag(D);
    mineig = min(eigcoefs);
    tol = 1e-6;
    if(mineig<0)
        mineig
        eigcoefs = eigcoefs + abs(mineig) + tol;
    end
    
    [sorteigs sorteig_idx] = sort(eigcoefs,'descend'); 
    Vsort = V(:,sorteig_idx);

    % Condition number shrinkage
    % Set regularization path kappa_max = 1 to condnum(Shat)
    lambda_min = 1.01;
    lambda_max = max(sorteigs)/min(sorteigs); 
    lambda_path = logspace(log10(lambda_min*min(sorteigs)),...
                            log10(1),opts.nlambdas) .* lambda_max;
    
    disp('Regularization Parameters')
    disp(lambda_path)
    
    shrinkage = estimator.cnshrinkage(sorteigs,lambda_path);
    shrinkage.path = lambda_path;
    
    for pathno=1:length(lambda_path)
        
        shrinkage.covariance_estimate{pathno} = ...
            Vsort*diag(shrinkage.eig_shrinkage{pathno})*Vsort';
        
        R = covariance.var_corr(...
              shrinkage.covariance_estimate{pathno});              
        shrinkage.correlation{pathno} = R;
        clear R;
            
        shrinkage.inverse_covariance_estimate{pathno} = ...
            Vsort*diag(1./shrinkage.eig_shrinkage{pathno})*Vsort';
        
        R = covariance.var_corr(...
                    shrinkage.inverse_covariance_estimate{pathno});
        
        shrinkage.partial_correlation{pathno} = 2*eye(length(R))-R;
    end
    
end


function fobj = CNobjective(eigcoefs,tau)
%CNobjective is the condition number regularized MLE objective
% 
% 
% Implements equations (14)-(16) from Won et. al. (2014)
    
% objfun = @(eigvalue,mu)(eigvalue*mu-log(mu)); 
% use mu = tau if l_i <= tau; => fobj = l_i*tau - log(tau)
% use mu = tau if l_i >= kappa*tau => fobj = l_i*tau - log(tau)
% use mu = 1/l_i if l_alpha < l_i < l_beta => fobj = 1 + log(l_i)

    
end