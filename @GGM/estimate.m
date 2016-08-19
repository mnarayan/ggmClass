function self = Estimate(self,varargin)
	% Choose between different GGM estimators 
	% 
	% Implement Single Stage
	% sparse MLE (via QUIC)
	% sparse PLE (Pseudo-Likelihood Estimator/Neighborhood Selection)
	% CLIME
	% DtraceLasso
	% Implement Two-stage
	% Apply any single stage, then re-fit constrained MLE
	% Apply any single stage, then adaptive MLE
	% Apply any single stage, use in SCAD/MCP LLE approximation
	% 
	
	
	
	switch func2str(self.methodfun)
		
	case 'sparseMLE'
		[self results] = self.sparseMLE();
	otherwise
		disp('Not yet supported');
	end
	
end