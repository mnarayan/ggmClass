function options = create_options(varargin)
    
    % Default options for estimator
    options.estimator = @sparseMLE;
    options.solver = 'QUIC';
    
    % Default options for the regularization path
    options.get_lambda = @(Sigma)({...
                            .01*min(min(abs(triu(Sigma(abs(Sigma)>0),1)))),...
                            max(max(abs(triu(Sigma,1)))); ...
                            });
    options.lambda_min = .01;
    options.lambda_max = 1.0;
    options.nlambdas = 30;
    options.lambdafun = @(lmin,lmax)(linspace(lmin,lmax,options.nlambda));
    
    % Default options for creating resamples or subsamples
    options.resampler.options = @generator.create_options
    options.resampler.run = @(n)(generate_resamples(...
                                            options.resampler.options(n) ...
                                            ));
    
    % Default options for model selection, 
    % - default StARS
    options.selection.score = @(Sigma,Theta)(stars(Theta));
    options.graphlet = 'dyad'; % edgelet
    options.instability_beta = .15;
    
end