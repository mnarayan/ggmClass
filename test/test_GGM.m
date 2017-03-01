function test_suite = test_GGM
	buildFunctionHandleTestSuite(localfunctions);

end
	
function test_GGM_chain

	load('test/test_lattice'); 
	
	% load simulation parameters
	m = simulation.simopts.m;
	p = simulation.simopts.p; 
	n = simulation.simopts.n;
	
	% Check test paramters
	assert(n==1,'Multiple subjects not supported');
	
	Data = mvnSim.generateNormal(1000,p,simulation.Theta);
	
	% Initialize
	GGMobj = GGM(Data);
	
	% Check sample covariance self.Sigma
	
	% Check variance self.W
	
	% Check Theta at given lambda value
	GGMobj.sparseMLE();
	% repeat
	GGMobj.estimate();
		
	% Check Theta over regularization path
	clear GGMobj;
	GGMobj = GGM(Data,1,0);
	GGMobj.sparseMLE();

	
	% UNTESTED
	% if(GGMobj.verbose)
	% 	gr_sparsity = squeeze(sum(sum(GGMobj.ThetaPath~=0,1),2))';
	% 	n_lambdas = min(25,floor(.25*numel(gr_sparsity)));
	% 	gr_sparsity(1:n_lambdas)
	% 	gr_sparsity = squeeze(sum(sum(GGMobj.AdaptTheta~=0,1),2))';
	% 	gr_sparsity(1:n_lambdas)
	% end
	% % Check plotting of regularization path
	% GGMobj.plot();
	
	% Check Model Selection
	currstate = rng;
	GGMobj.randState = currstate;
	GGMobj.verbose = 0;
	[GGMobj grphs] = GGMobj.stableEstimate(currstate,@stars,10);
	
	% Check plotting of stability path
	
		
	% Check tylerMLE
	GGMobj = GGM(Data);
	[GGMobj results] = GGMobj.tylerMLE();
	
	% Check debias
	GGMobj = GGM(Data); 
	GGMobj.debias();
	
	clear n_lambda;	
end

% 	% Initialize Variables
% 	% m = 100;
% 	% p = 50; % 100 and 150
% 	% n = 1;
% 	load testcase1
%
% 	%Initialize non-zero patterns;
% 	M1 = eye(m);
% 	M2 = eye(p);
% 	%M2 = smallw(p,2,.25,1);
%
% 	SimulateMatrixVariate;
% 	X = zeros(m,p,n);
% 	% Matrix-Variate
% 	for ii=1:n
% 		Z = randn(m,p);
% 		X(:,:,ii) = A*Z*B';
% 	end
% 	Lambda = .1;
%
% 	disp('Testing Diag-Glasso Estimator ...');
% 	try
% 		[Grph1a Pi1a Omega1a lambdas_a Ka] = QuasiFF_Gemini_Model1(X,.5,n,0,0,Lambda,1);
% 		metrics = ComputeMetrics(Grph1a,M2,1);
%
% 		% Test Row Supp
% 		assertElementsAlmostEqual(Omega1a,M1, 'absolute',.015,sqrt(eps), 'Incorrect Row Supp for testcase 1, M1 = I');
%
% 		% Test F1-Score
% 		assertAlmostEqual(metrics.f1, 1,.2, 'Large Difference in F1 for testcase 1, M1 = I');
% 	catch
% 		%error(sprintf('QuasiFF_Gemini under Model 1, Row Diag option failed with Lambda %.2f',Lambda));
% 		error(['GEMINI for Model1, Row Diag: ' ...
% 			' QuasiFF_Gemini_Model1(X,.5,s1 = %d,s2 = %d,%d,Lambda = %.2f,isRowDiag=%d)'], ...
% 			n,0,0,Lambda,1);
% 	end
%
%
% 	disp('Testing Oracle-Gemini Estimator ...');
% 	try
% 		[Grph1b Pi1b Omega1b lambdas_b Kb] = QuasiFF_Gemini_Model1(X,.5,n,0,0,Lambda,0,M1);
% 		metrics = ComputeMetrics(Grph1a,M2,1);
%
% 		% Test Row Supp
% 		assertElementsAlmostEqual(Omega1b,M1, 'absolute',.015,sqrt(eps), 'Incorrect Row Supp for testcase 1, M1 given');
%
% 		% Test F1-Score
% 		assertAlmostEqual(metrics.f1,1,.2, 'Large Difference in F1 for testcase 1, M1 = I');
%
% 	catch
% 		error(['GEMINI for Model1, AR1 Diag: ' ...
% 			' QuasiFF_Gemini_Model1(X,.5,s1 = %d,s2 = %d,%d,Lambda = %.2f,isRowDiag=%d)'], ...
% 			n,0,0,Lambda,0);
% 	end
%
% 	disp('Testing KGGM class constructor with no banding ...')
% 	try
% 		testcase1 = KGGM(X);
% 		testcase1.verbose=0;
% 		testcase1.fit_banding(0,1);
%
% 	catch
% 		error('KGGM(X) constructor with data matrix X failed');
% 	end
%
% 	% Test Row Supp
% 	assertElementsAlmostEqual(1*((testcase1.Theta1~=0)), abs(M1~=0), 'absolute',.015,sqrt(eps), ...
% 							'Incorrect Row Supp for testcase 1, M1 = AR-1');
%
% 	% Test F1 Score on Row Supp
% 	metrics = ComputeMetrics(testcase1.Theta1,M1,0);
% 	assertAlmostEqual(metrics.f1,1,.01, 'Large Difference in F1 for testcase 1, M1 = I');
%
% end