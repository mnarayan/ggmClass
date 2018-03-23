function output = model_average_populations(S,varargin)
% ESTIMATOR.MODEL_AVERAGE_POPULATIONS
% Model average estimation given an array of sample covariance or correlation matrices
% 
% INPUT
%  - S is a cell array (length n) or matrix of p x p  covariance matrices S1,S2,...,Sn
% 

    if(nargin==1)
        estimate_options = estimator.create_options();
    else
        estimate_options = varargin{1};
    end

    if(ndims(S)~=3)
        error('S should consist of more than 1 correlation matrix');
    else
        if(size(S,3)<10)
            error('Insufficient sample size n (no. of covariances)')
        end
    end

    % Get resampling iterator: 
    % options include bootstrapping, subsampling, CP intersections
    n_samples = size(S,3);
    resampling_options = estimate_options.resampler.options(n_samples);
    resampler = estimate_options.resampler.run(resampling_options);
    n_resamples = resampler.options.B;

    % Step 1: Resampling
    ensemble_results = {};
    graphs = {};

    function pooledShat = group_pooled_covariance(Sn)
        pooledShat = mean(Sn,3);
    end
    esimate_options.covariancefun = @group_pooled_covariance;

    for bb=1:n_resamples
        % - a) for each resample X*
        % - b) Modify weights for re-weighted estimators
        % - c) call estimator.fit() on X*
        ensemble_results = estimator.fit(...
                                        S(:,:,resampler.samples(bb,:)),...
                                        estimate_options...
                                        );
        graphs{bb,1} = ensemble_results.inverse_covariance_estimate;
        cp_samples = setdiff([1:n_samples],resampler.samples(bb,:));
        ensemble_results = estimator.fit(...
                                        S(:,:,cp_samples),...
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
    output.refit_estimator = estimator.fit(S,refit_options);

    output.fit_options = estimate_options;
    output.resampler = resampler;
end