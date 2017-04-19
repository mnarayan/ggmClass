function [samples options] = generate_samples(options)
%GENERATE_SAMPLES is internal helper function
	
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

	rng(options.rngstate); 
	samples = options.generatorfun(); 
	
end