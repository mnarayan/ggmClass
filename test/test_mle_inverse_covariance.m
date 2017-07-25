function test_suite = test_GGM
	buildFunctionHandleTestSuite(localfunctions);

end


function test_diagonal_matrix
    
    load('test/test_lattice')
    
	% load simulation parameters
	m = simulation.simopts.m;
	p = simulation.simopts.p; 
    n = 1;  
    Data = simulate.matrixvariate(m,p,eye(m),eye(p)); 
    
    
    
end

	
% function test_GGM_chain
%
%     load('test/test_lattice')
%
%     % load simulation parameters
%     m = simulation.simopts.m;
%     p = simulation.simopts.p;
%     n = 1;
%     popTheta = simulate.make_psd(simulation.Theta);
%     popTheta = simulate.varcorr(popTheta);
%     popChol = simulate.mcd(popTheta);
%     popVar = eye(p);
%     popTheta = popChol*popVar*popChol';
%     Data = simulate.matrixvariate(m,p,eye(m),popTheta);
%
%     % Initialize
%     GGMobj = GGM(Data);
%     results = {};
%
%     % Check sample covariance self.Sigma
%     [GGMobj results.mlecov] = GGMobj.MLECovEstimate();
%
%     % Check variance self.W
%
%     % Check Theta at given lambda value
%     [GGMobj results.sparsemle] = GGMobj.sparseMLE();
%     [ThetaRefit results.constrainmle] = GGM.constrainMLE( ...
%                                         GGMobj.Sigma, ...
%                                         GGMobj.Theta, ...
%                                         1.0);
%
%     % repeat
%     GGMobj.estimate();
%
%     % Check Theta over regularization path
%     clear GGMobj;
%     GGMobj = GGM(Data,1,0);
%     [GGMobj results.sparsemle] = GGMobj.sparseMLE();
%
%
% end