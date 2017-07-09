function options = create_options(n,varargin)
%GENERATOR.CREATE_OPTIONS creates options for resampling or subsampling
% 
% SEE ALSO GENERATOR.GENERATE_RESAMPLES, GENERATOR.GENERATE_SUBSAMPLES

	options = {};
	options.seed = 0;
	options.streamstate = {};
	options.rngstate = rng('default'); 	% state before calling resamples
	options.B = 50; 					% Number of subsamples or resamples
	options.sizefun = ... 				% To set Bsize for general m of n
		@(m,n)( ...
		(m<n)*ceil(m) + (m>=n)*n ...
		);
	options.blocklength = 2*ceil(n^(1/3)); 	% Default block length when use_blocks is true;
	options.folds = 2;
	options.groups	= ones(1,n); 		% Single group for all observations 
    
    %%%% Flag Options %%%%%
	options.use_replacement = true; 	% Bootstrapping or Subsampling
	options.use_mofn = true;			% Use m out of n
	options.use_blocks	= false;		% Call block bootstrap/subsamples
	options.stratify = false;
	options.use_fixedblock = true; 		% Fixed vs. variable block length
	options.use_circular  = false; 	    % Make blocks circular, wrap around 
    options.verbose = false;
    setflags = {'use_replacement',...
                'use_mofn',...
                'use_blocks',...
                'stratify',...
                'use_fixedblock',...
                'use_circular',...
                'verbose'}; 
    
    default_opts = options;    
	if(nargin==1)
		disp('Using default options')	
	else
		disp('Processing custom options not tested. Use with caution');
        if(exist('parseArgs'))
            options = parseArgs(varargin,default_opts,setflags);
        elseif(exist('ezparse'))
            default_flag_opts = default_opts;
            for flagNo=1:length(setflags)
                default_flag_opts = setfield(default_flag_opts, ...
                                setflags{flagNo},'%FLAG%');
            end
            [options,~] = ezparse(default_flag_opts,varargin);
            for flagNo=1:length(setflags)
                options.(setflags{flagNo}) = options.(setflags{flagNo}) | ...
                                             default_opts.(setflags{flagNo});
            end
        end
        if(options.verbose)
            disp('Custom options used')
            options	
        end
	end	
	
	if(options.use_replacement)
		if(options.use_mofn)
			options.Bsize = options.sizefun(.8*n,n); 
		else
			options.Bsize = options.sizefun(n,n); 
		end
	else
		% Subsampling is always m out of n
		options.use_mofn = true;
		if(n>=100)
			options.Bsize = options.sizefun(n/options.folds,n); 
		else
			options.Bsize = options.sizefun(.8*n,n); 
			warning('Might not have adequate subsample size'); 
		end
	end
	
	if(options.use_fixedblock)
		options.blockfun = @generator.create_fixedblocks;
	else
		options.blockfun = @generator.create_fixedblocks;
	end
	options.generatorfun = '';		
	
	options.n = n;		
	options.currstate = {}; % state after calling resamples
			
end