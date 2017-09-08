function lk = mvn_loglikelihood(Theta,Sigma,n)
% Loglikelihood for Multivariate Gaussian/Normal 
% 
%   LL = n * (1/2) * {log(det(Theta)) - trace(Sigma*Theta)}
    
    U = chol(Theta);
    logdet = 2*sum(log(diag(U))); 
    lk = n/2 * (logdet - trace(Sigma*Theta));

end