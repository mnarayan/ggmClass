function [samples varargout] =  create_complementarypairs(options)
%CREATE_COMPLEMENTARYPAIRS converts resamples or subsamples into corresponding disjoint sets 
	
    disp('Uses Half Sample Splitting')
    options.use_replacement = false;
    options.folds = 2;

    output = generator.create_repeatedkfold(options);
    samples = output.partition;

    if(nargout==2)
        varargout{1} = output;
    end
	
end