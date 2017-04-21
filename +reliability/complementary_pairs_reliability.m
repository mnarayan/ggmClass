function output = complementary_pairs_reliability(X,options)
%COMPLEMENTARY_PAIRS_RELIABILITY repeatedly applies some custom metric to
% disjoint data splits such as 2-fold CV splits using a custom function and
% returns reliability of the metric using a custom function
% 
% Requires  GENERATOR.CREATE_REPEATEDKFOLD
% 
% USAGE:   output = complementary_pairs_reliability(X, ...
% 					struct(...
% 					'metricfun', @(x)(sum(sum(corr(x)))),...
% 					'reliabilityfun', @kendallsW));
% INPUTS
% X - measurements x variables data matrix. 
% options.metricfun -  Custom function that takes n x p data matrix 
%	as input and computes a metric.
% options.reliabilityfun - Custom function takes samples x raters
%	and returns reliability across raters for each feature
% 	Example using ICC: options.reliabilityfun = @(x)(IPN_icc(x,2,'k'));
%
% OUTPUTS
% output.cp_scores - contains metric obtained over each replicate x fold
% output.cp_reliability - contains reliability scores of the same length as the output from metricfun

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Copyright 2017 Manjari Narayan %%
%% BSD-2-Clause License           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if(ndims(X)==2)
		[n p] = size(X);
	elseif(ndims(X)>2)
		error('Only 2-Dimensional Data Matrix supported')
	end

	if(~isfield(options,'metricfun')|~isfield(options,'reliabilityfun'))
		error('please specify valid metricfun and reliabilityfun in options parameter');
	end
	
	sampling_options = generator.create_options(n); 
	sampling_options.folds = 2; 
	sampling_options.use_replacement = false;
	output = generator.create_repeatedkfold(sampling_options); 
	partition = output.partition;
	
	try
		tmp_metric = options.metricfun(X);
	catch me
		fprintf('metricfun cannot process input data matrix')
		disp(me)
	end
	n_features = length(tmp_metric);
	cp_scores = zeros(length(partition),options.folds,n_features);
	for trialno = 1:length(partition)
		for foldno = 1:options.folds
			samples = partition{trialno}.samples(:,foldno);
			cp_scores(trialno,foldno,:) = options.metricfun(X(samples,:));
		end
	end
	
	cp_reliablity = zeros(1,n_features); 
	for featureno=1:n_features
		cp_reliability(featureno) = options.reliabilityfun(cp_scores(:,:,featureno));
	end
	
	output.metricfun = options.metricfun;
	output.reliabilityfun = options.reliabilityfun;
	output.cp_scores = cp_scores;
	output.cp_reliability = cp_reliability;
	
end
		
		