function [results] = ebic(ThetaList,X,options)	
% EBIC - returns extended BIC scores and best regularization parameter from  a family of graphs along the regularization path.
%    EBIC(E(\lambda)) = -2 * loglikelihood + |E| log(n) + 4 |E(\lambda)| gamma log(p), 
%       where E is the non-zero edgelist at regularization lambda
% Input
%   - ThetaList : inverse_covariance_estimate p x p x nlambdas
%   - X : n x p data matrix
%   - options : options structure that follows estimator.create_options()

    results = options.ebic;
    
    [n p] = size(X); 
    [~,~,nlambdas] = size(ThetaList); 
    S = options.covariancefun(X);
    scores = nan(1,nlambdas);
    sparsity = nan(1,nlambdas);
    lambda_opt = NaN;

    for lambdano=1:nlambdas
        sparsity(lambdano) = (sum(sum(abs(ThetaList(:,:,lambdano))~=0))-p);
        scores(lambdano) = ...
                    -2 * ...
                    options.ebic.loglikefun(ThetaList(:,:,lambdano),S,n) + ...
                    sparsity(lambdano) * log(n) + ...
                    4 * sparsity(lambdano) * options.ebic.gamma * log(p);
    end
    
    [min_val,min_idx] = min(scores);
     
    results.sparsity = sparsity;
    results.scores = scores;
    results.lambda_opt_idx = min_idx;
    results.lambda_opt = options.path(min_idx);

end