function output =  create_repeatedkfold(options)
%CREATE_REPEATEDKFOLD creates repeated k-fold partitions
% 
% INPUTS
% 	options.n 			- sample size
%
% 	(Optional)
% 	options.folds		- Number of folds. Default 2. 
%  	options.B			- Number of trials/repetitions
% 	options.stratify 	- true or false
% 	options.group 		- group indicators for stratification

	output.state = rng;
	partition = {};
	if(exist('cvpartition'))	
		for ii=1:options.B
			partition{ii}.state = rng;			
			tmp_partition = cvpartition(options.n,'KFold',options.folds);
			for foldno=1:options.folds
				partition{ii}.samples(:,foldno) = tmp_partition.training(foldno);
				% partition{ii,foldno}.test(:,foldno) = tmp_partition.Test(foldno);
			end
			clear tmp_partition;
		end
	else
		disp('Call custom partitioner')
	end
	
	output.partition = partition;
	output.options = options;
	
end