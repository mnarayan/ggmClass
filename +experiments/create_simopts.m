function options = create_simopts(varargin)
    
    options = {};
    
    % Dimensionality and sample size
    options.m = 100; 
    options.p = 25;
    options.n = 1;
    
    % Monte Carlo Trials
    options.ntrials = 2;
    options.seed = 0; 
    options.getstate = @(seed)(RandStream('mt19937ar','Seed',seed));
    options.resetstate = @()(RandStream.setGlobalStream(...
                                            options.getstate(0) ...
                                            ));
    
	if(options.seed==0)
	    s = RandStream('mt19937ar','Seed',0);
	    RandStream.setGlobalStream(s);
		options.streamstate = get(s); 
		options.rngstate = rng('default'); 
	else
	    s = RandStream('mt19937ar','Seed',options.seed); 	
		options.streamstate = getGlobalStream(s); 
		options.rngstate = rng(options.rngstate);	
	end 
    
    % Adjacencies
    options.type = 'rowcolumn';  % or 'col' only.     
    options.simulator.row_options = simulate.create_options;
    options.simulator.col_options = simulate.create_options;
    options.simulator.row_options.ndims = options.m;
    options.simulator.row_options.nsamples = options.p;   
    options.simulator.col_options.ndims = options.p;
    options.simulator.col_options.nsamples = options.m;
    % Adjacency simulation for pggm-sims
    options.row_adjacency = options.simulator.row_options.adj_fun(options.m); 
    options.col_adjacency = options.simulator.col_options.adj_fun(options.p);
    
    % Estimator
    options.estimator.row_options= estimator.create_options();
    options.estimator.col_options = estimator.create_options();    
    
    
end