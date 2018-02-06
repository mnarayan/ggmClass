function fitobj = concord(X,varargin)
    
    fitobj = {};
    
    [Sigma Theta opts] = coordinatewise_concord(X, varargin{:}); 
    
    fitobj.covariance = Sigma;
    fitobj.inverse_covariance = Theta;
    Theta_Dinv = diag(1./sqrt(diag(Theta))); 
    fitobj.partial_correlation = eye(size(Theta,1)) - ...
                                     Theta_Dinv*Theta*Theta_Dinv;
    fitobj.opts = opts;

end


function [Sigma Theta opts] = coordinatewise_concord(X,varargin)
    
    if(nargin>=2)
        opts = varargin{1}; 
    else
        opts = create_options();
    end
    
    addpath('../');
    [n p] = size(X); % n samples x p features;
    Shat = opts.covfun(X); 
    % X = standardize.standardize_cols(X,'center');
    % Shat = cov(X);
    % Shat = covariance.mle_sample_covariance(X, ...
    %                             struct('standardize','cols'));
    if(opts.kroot>1)
        [V D] = eig(Shat);
        Shat = V * diag(diag(D).^(1/(opts.kroot))) * V';  
    end                  
    Shat_D = reshape(diag(Shat), [p 1]); % vector p x 1
    Shat_X = Shat - diag(Shat_D); % matrix p x p
    if(isvector(opts.lambda))
        lambda = opts.lambda .* max(max(abs(triu(Shat))));
    else
        lambda = opts.lambda;
    end
    % Initialize
    iter = 0; 
    converged = false;
    %Theta = randn(p,p); Theta(find(eye(p))) = 1.0; 
    %Theta = diag(1./(diag(Shat))); 
    %Theta = eye(p); 
    Theta = pinv(Shat); 
    Theta_D = diag(diag(Theta)); Theta_D_update = Theta_D; 
    Theta_X = Theta - Theta_D; Theta_X_update = Theta_X; 
    
    while (iter < opts.max_iter && ~converged )
        % Off-Diagonal Update
        Theta_D = Theta_D_update;
        Theta_X = Theta_X_update;

        for pp=1:p
            for qq=pp+1:p
                symTrST(pp,qq) = Shat(pp,:)*Theta(:,qq) + ...
                                 Theta(pp,:)*Shat(:,qq); 
                symTrST(pp,qq) = symTrST(pp,qq) - ...
                            Shat(pp,pp)*Theta(pp,qq) - ...
                            Theta(pp,qq)*Shat(qq,qq);
                symTrST(pp,qq) = symTrST(pp,qq) ./ ...
                            (Shat_D(pp) + Shat_D(qq));
                if(isscalar(lambda))
                    normLambda = lambda ./ ((Shat_D(pp) + Shat_D(qq)));
                elseif(all(size(lambda)==size(Theta)))
                    normLambda = lambda(pp,qq) ./ ((Shat_D(pp) + Shat_D(qq)));
                end
                Theta_X_update(pp,qq) =  ...
                        opts.shrinkage( - symTrST(pp,qq), normLambda) ./ ...
                            (Shat_D(pp) + Shat_D(qq));
                Theta_X_update(qq,pp) = Theta_X_update(pp,qq);
                Theta(pp,qq) = Theta_X_update(pp,qq); 
                Theta(qq,pp) = Theta_X_update(qq,pp);
            end
        end
        
        if(iter<=2 && opts.verbose)
            imagesc(-symTrST); colorbar;
            pause(1);
        end
        
        % % Diagonal Update
        %Theta_D_update = zeros(p,1);
        Theta_D_update = diag(Theta_D_update);
        for pp=1:p
            TrST = Shat_X(pp,:)*Theta(:,pp);
            TrST = TrST - Shat_X(pp,pp)*Theta(pp,pp);
            Theta_D_update(pp) =  ( - TrST + ...
                                sqrt(TrST.^2 + 4*Shat_D(pp))) ./ (2*Shat_D(pp));
            Theta(pp,pp) = Theta_D_update(pp);                    
        end
        Theta_D_update = diag(Theta_D_update);
        
        Theta = Theta_X_update + Theta_D_update;
        
        err(iter+1,1) = norm(diag(Theta_D-Theta_D_update),'fro')/p; 
        err(iter+1,2) = norm(Theta_X-Theta_X_update,'fro')/(p*(p-1));
        
        objective(iter+1) = 1/n * loss.concord(X,Shat,Theta) + ...
                                 sum(sum(lambda .*abs(triu(Theta))));
         
        iter = iter + 1;                            
        if(sum(err(iter,:)) < opts.tol)
            converged = true;
        end
        if(opts.verbose && mod(iter,5)==0 )
            disp(sprintf('Iteration: %d, Conv. Err: %.6f, Obj: %.4f', ...
                               iter, sum(err(iter,:)), objective(iter)));
        end
        
    end
    
    Theta = Theta_D_update + Theta_X_update; 
    
    if(opts.kroot>1)
        Theta = Theta^(opts.kroot);
        Shat = V * diag(diag(D)) * V'; 
    end
        
    Sigma = Shat;
    opts.lambda = lambda;
    opts.iter = iter;
    opts.err = err; 
    opts.objective = objective;

end



function options = create_options()
    
    options = {};
    options.solver_type = 'coordinate';
    options.covfun = @(X)(covariance.mle_sample_covariance(X, ...
                                 struct('standardize','cols')));
    options.shrinkage = @(a,kappa)(shrinkage(a,kappa));
    options.kroot = 1;
    options.lambda = .2;
    options.max_iter = 1000; 
    options.tol = 1e-6;
    options.verbose = false;
    
end


function y = shrinkage(a, kappa)
	if(ismatrix(a) && ~isscalar(a))
		a_diag = diag(a);
        savediag = true;
    else
        savediag = false;
    end
    y = max(0, a-kappa) - max(0, -a-kappa);
    % if(a>0)
    %     y = max(0,a - kappa);
    % elseif(a<=0)
    %     y = min(0,a + kappa);
    % end
	if(~isscalar(y) && savediag)
        p = size(a,1); 
		y(find(eye(p))) = a_diag; 
	end
    %disp(sprintf('Val: %.2f, lambda: %.2f, STVal: %.2f',a,kappa,y));
end