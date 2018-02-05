function loss = concord(X, Sigma, Theta, varargin)
    % LOSS.CONCORD returns the loss function for the CONCORD pseudolikelihood estimator
    % 
    %  - sum_ii n * log(theta_ii) + 1/2 * trace(Theta*X', X*Theta)
    
    [n p] = size(X); 
    X = standardize.standardize_cols(X,'center'); 
    
    Theta_D = diag(Theta); 
    Theta_R = Theta - diag(Theta_D); 
    Sigma_D = diag(Sigma); 
    Sigma_R = Sigma - diag(Sigma_D); 
    
    logdetST = - n * sum(log(Sigma_D)) ; 
    
    % trST = 0;
    % for pp=1:p
    %     trST = trST +  .5 * (Theta(:,pp)' * Sigma * Theta(:,pp)) ;
    % end
    trST = .5 * n * trace(Theta' * Sigma * Theta); 
    
    loss = logdetST + trST; 
    
end