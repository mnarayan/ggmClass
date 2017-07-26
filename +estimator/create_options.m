function options = create_options(varargin)
    
    % Default what is the model 
    options.model = 'gaussian';
    
    % Default options for estimator
    options.covariancefun = @(X)(covariance.mle_sample_covariance(X,...
                                        struct('standardize','cols') ...
                                        ));
    options.init_estimator = 'lasso_mle';                                  
    options.estimator = 'lasso_mle';
    options.solver = 'QUIC';
    options.refit = false;
    
    % Default options for the regularization path
    options.get_lambda = @(Sigma)({...
                            .1*min(abs(Sigma(abs(triu(Sigma,1))>0))),...
                            max(max(abs(triu(Sigma,1)))) ...
                            });
    options.lambda_min = .01;
    options.lambda_max = 1.0;
    options.nlambdas = 30;
    options.lambdafun = @(lmin,lmax)(fliplr(linspace(...
                                            lmin,...
                                            lmax,...
                                            options.nlambdas...
                                            )));
    options.Lambda = 0;
    options.W = [];
    options.path = options.lambdafun(options.lambda_min,options.lambda_max);
    
    % Default options for creating resamples or subsamples
    options.resampler.options = @generator.create_options;
    options.resampler.run = @(n)(generate_resamples(...
                                            options.resampler.options(n) ...
                                            ));
    
    % Default options for model selection, 
    % - default StARS
    options.selection.score = @(Sigma,Theta)(stars(Theta));
    options.graphlet = 'dyad'; % edgelet
    options.instability_beta = .15;
    
    
end