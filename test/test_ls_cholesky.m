function test_suite = test_ls_cholesky
	buildFunctionHandleTestSuite(localfunctions);

end


% function test_diagonal_matrix
%
%     load('test/test_lattice')
%
%     % load simulation parameters
%     m = simulation.simopts.m;
%     p = simulation.simopts.p;
%     n = 1;
%     Data = simulate.matrixvariate(m,p,eye(m),eye(p));
%
%     options = estimator.create_options();
%     options.refit = false;
%
%
% end


function test_stationary_arsim

    addpath('external/arfit'); 

    load('test/test_lattice')    
	% load simulation parameters
	m = simulation.simopts.m;
	p = 10*simulation.simopts.p; 
    n = 1;  
       
    E = eye(p);
    mu = zeros(p,1);
    lags = 1;
    rho = .3;
    A = eye(p)*rho;
    Y = arsim(mu,A,E,m);
    
    [T D] = estimator.ls_cholesky(Y); 
    assertElementsAlmostEqual(-T(2,1),rho,'absolute',1e-2);
    assertElementsAlmostEqual(D(1,1),1,'absolute',1e-2);
    
    
end