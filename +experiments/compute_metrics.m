function metrics = compute_metrics(G,G0,varargin)  
% Computes Relative Frobenious Norm, True Positive Rate, False Positive Rate, False Negative Rate and F-1 Score
% Inputs
% 	G 	- Estimated Graph
% 	G0 	- Population (True) Graph
% 	verbose (optional) - display metrics
% 	
% Outputs
%	metrics - table of metrics
	
	if(nargin==2)
		verbose = 1;
	else
		switch nargin
		case 3
			verbose = varargin{1};
		end
	end
	
	nunits = size(G,3);
    
    offDiag = eye(size(G,1))==0;
    
    InfNorm = squeeze(max(max(offDiag.*abs(bsxfun(@minus,G,G0)),[],1),[],2));
    
	TP = squeeze(sum(sum(bsxfun(@and,G~=0,G0~=0),1),2)/sum(sum(G0~=0,1),2));

	FP = squeeze(sum(sum(bsxfun(@and,G~=0,G0==0),1),2)/sum(sum(G0==0,1),2));

	FN = squeeze(sum(sum(bsxfun(@and,G==0,G0~=0),1),2)/sum(sum(G0~=0,1),2));

	F1 = squeeze(2*TP./(2*TP + FP + FN));

	S = squeeze(sum(sum(G~=0,1),2));
	
    if(nunits==1)
        froberr = norm(G-G0,'fro')/trace(G0);
        spectral = norm(G-G0,2); 
        condnum = cond(G); 
        graphNorm = trace(G*G)/trace(G);
        metrics.spectralErr = spectral;
        metrics.condnum = condnum; 
        metrics.graphNorm = graphNorm;
        if(verbose)
            disp(sprintf('Spectral Err. %.2f',spectral));
            disp(sprintf('Off diag. Norm %.2f',graphNorm));
            disp(sprintf('Cond. Num %.2f',condnum));
        end            
    else
    	froberr = squeeze(...
    	      sum(sum(bsxfun(@minus,G,G0).^2,1),2)/sum(sum(G0.^2,1),2));        
        metrics.spectralErr = nan(nunits,1);
        metrics.condnum = nan(nunits,1);
        metrics.graphNorm = nan(nunits,1);
        
    end
	
    metrics = table();
    
    metrics.frobeniusErr = reshape(froberr,[nunits 1]);
    metrics.matrixinfErr = reshape(InfNorm,[nunits 1]);
    metrics.truepos = reshape(TP,[nunits 1]);
    metrics.falsepos = reshape(FP,[nunits 1]); 
    metrics.falseneg = reshape(FN,[nunits 1]);
    metrics.f1 = reshape(F1,[nunits 1]);
    metrics.sparsity = reshape(S,[nunits 1]);
    
    if(verbose)
        disp(sprintf(['Relative Frobenius Error %.2f ' ...
                     '\nTrue Positive Rates %.2f'],  ...
                    mean(froberr),mean(TP)));
        disp(sprintf('Max. Entrywise Err. %.2f',mean(InfNorm)));
        disp(sprintf('False Positive Rates %.2f',mean(FP)));
        disp(sprintf('False Negative Rates %.2f',mean(FN)));
        disp(sprintf('F-1 Score %.2f',mean(F1)));
    end

end