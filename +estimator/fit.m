function [results] = fit(X,options)
%ESTIMATOR.FIT() Fits graphical model estimator specified by options
% 

    results = {};
	
    if(isempty('options'))
        options = create_options();
    end
    results.options = options;
    
    switch options.model
        
	case 'gaussian'
        [covariance_estimate, inverse_covariance_estimate, sample_covariance, fit_options] = estimator.penalized_mle_inverse_covariance(X,options);
	otherwise
		disp('Not yet supported');
	end
    
    if((options.Lambda == 0) && (options.nlambdas > 1))
        results.ebic = estimator.ebic(inverse_covariance_estimate,...
                                        X,...
                                        fit_options...
                                        );
    end
    results.covariance_estimate = covariance_estimate;
    results.inverse_covariance_estimate = inverse_covariance_estimate;
    results.sample_covariance = sample_covariance;
    results.fit_options = fit_options;
    
end