function test_suite = test_flip_flop
	buildFunctionHandleTestSuite(localfunctions);

end


function test_diagonal_matrix
    
    load('test/test_lattice');
    
	% load simulation parameters
	m = simulation.simopts.m;
	p = simulation.simopts.p; 
    n = 1;  
    Data = simulate.matrixvariate(m,p,eye(m),eye(p)); 
    
    options.nlambdas = 1;
    options.Lambda = .05;
    options.refit = true;
    
    row_options = options;
    col_options = options;
    iterative = false;
    results = estimator.flip_flop(Data,row_options,col_options,iterative);
    
    
end


function test_eye_chain_matrix
    
    load('test/test_lattice');
    
	% load simulation parameters
    m = 100;
	p = simulation.simopts.p; 
    n = 1;  
    popTheta = simulate.make_psd(simulation.Theta);
    popTheta = simulate.varcorr(popTheta);
    popChol = simulate.mcd(popTheta);
    popVar = eye(p);
    popTheta = popChol*popVar*popChol';
    Data = simulate.matrixvariate(m,p,eye(m),popTheta);
    
    options.nlambdas = 1;
    options.Lambda = .05;
    options.refit = true;
    
    row_options = options;
    col_options = options;
    iterative = false;
    results = estimator.flip_flop(Data,row_options,col_options,iterative);
    
end