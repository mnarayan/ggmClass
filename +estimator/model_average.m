function output = model_average(X,varargin)
% ESTIMATOR.MODEL_AVERAGE
% Model average estimation given a single data matrix
% 
    
    if(nargin==1)
        estimate_options = estimator.create_options();
    else
        estimate_options = varargin{1};
    end
    
    % Get resampling iterator: 
    % options include bootstrapping, subsampling, CP intersections
    n_samples = size(X,1);
    resampling_options = estimate_options.resampler.options(n_samples);
    resampler = estimate_options.resampler.run(resampling_options);
    n_resamples = resampler.options.B;
    
    % Step 1: Resampling
    ensemble_results = {};
    graphs = {};
    for bb=1:n_resamples
        % - a) for each resample X*
        % - b) Modify weights for re-weighted estimators
        % - c) call estimator.fit() on X*
        ensemble_results = estimator.fit(...
                                        X(resampler.samples(bb,:),:),...
                                        estimate_options...
                                        );
        graphs{bb,1} = ensemble_results.inverse_covariance_estimate;
        cp_samples = setdiff([1:n_samples],resampler.samples(bb,:));
        ensemble_results = estimator.fit(...
                                        X(cp_samples,:),...
                                        estimate_options...
                                        );
        graphs{bb,2} = ensemble_results.inverse_covariance_estimate;
    end
    % Step 2: Aggregation
    %
    % - a) Ensemble Average Model Selection 
    [stars_scores stars_lambda stability_graphs] = ...
                             estimator.stars(graphs,estimate_options);
    output.stars.scores = stars_scores;
    output.stars.lambda_opt_idx = stars_lambda;
    output.stars.lambda_opt = estimate_options.path(stars_lambda);
    % - b) Get final re-fitted/relaxed estimator
    output.stability_graphs = stability_graphs;
    refit_options = estimate_options;
    refit_options.refit = true;
    refit_options.nlambdas = 1;
    refit_options.Lambda = output.stars.lambda_opt;
    output.refit_estimator = estimator.fit(X,refit_options);
    
    output.fit_options = estimate_options;
    output.resampler = resampler;
    
end