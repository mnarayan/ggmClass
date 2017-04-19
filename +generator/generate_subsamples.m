function resampler = generate_subsamples(varargin)

	if(nargin==0)
		opts = generator.create_options(100);
	else
		opts = varargin{1}; 
	end
	
	if(opts.use_replacement)
		warning('Replacement is on. Setting to false')
		opts.use_replacement = false;
	end
	
	% Subsampling is always m out of n. Default to half-sampling
	opts.use_mofn = true;
	if(opts.n>=100)
		disp('Half-sampling')
		opts.Bsize = opts.sizefun(opts.n/2,opts.n); 
	else
		opts.Bsize = opts.sizefun(.8*opts.n,opts.n); 
		warning('Might not have adequate subsample size'); 
	end
	
	function [submat] = gensubsamples(N,n,b)

		submat = zeros(N,b);
		for ii = 1:N
			submat(ii,:) = sort(randsample(n,b,0),'ascend');
		end

	end
	
	function [submat] = fixedblocksampler(N,Bsize,blockmat);

		[n_blocks blocklength] = size(blockmat); 
		
		% Ensuring submat Bsize requirement is met
		nblks = ceil(Bsize/(blocklength+1));	% # blocks to resample
		submat = zeros([N nblks*blocklength]);	% Ensure Bsize is met
		
		submat = [];
		blockresamples = gensubsamples(N,n_blocks,nblks);
		for ii = 1:N
			tmprs = blockmat(blockresamples(ii,:),:);
			submat(ii,:) = sort(reshape(tmprs,[1 numel(tmprs)]));
		end
		
	end
	
	function [submat] = variableblocksampler(N,Bsize,blockmat,var_window)
		
	    [tot_blocks blocklength] = size(blockmat); 		
		n_blocks = ceil(Bsize/(blocklength+1));	% # blocks to resample
		%submat = zeros([N nblks*blocklength]);	% Ensure Bsize is met
		submat = cell(1,N); 
		
		if(var_window>blocklength)
			warning('Length of variable block component too large'); 
			disp('Setting to max variable length to .25 of blocklength'); 
			var_window = max(3,ceil(blocklength/3);
		end
		
		blockresamples = gensubsamples(N,tot_blocks,n_blocks);
		
		for ii = 1:N
			tmp_resamples = [];
			random_offset = randsample(ceil(var_window)/2,n_blocks,1); 
			random_center = randsample(ceil(blocklength/2)-random_offset:ceil(blocklength/2)+random_offset,n_blocks,1);
			for nn=1:n_blocks
				tmp_win = floor((blocklength-var_window)/2);
				start_idx = max(1, ...
								floor(random_center(nn)-tmp_win)); 
				stop_idx = min(blocklength, ...
								ceil(random_center(nn)+tmp_win));
				tmp_blkresample = blockmat(blockresamples(ii,nn),:);
				tmp_blkresample = sort(reshape(tmp_blkresample,[1 numel(tmp_blkresample)]));
				tmp_resamples = [tmp_resamples tmp_blkresample(start_idx:stop_idx)];
			end
			
			submat{ii} = tmp_resamples;
		end
		
	end
	
	% Choose relevant generator function 
	if(opts.use_blocks)
		disp('Variable Block subsampling not yet implemented')
		blockmat = opts.blockfun(opts.n,opts.blocklength);
		opts.blockmat = blockmat;								
		if(opts.use_fixedblock)
			opts.generatorfun =  @()( ...
									fixedblocksampler(opts.B, ...
												opts.Bsize, ...
												blockmat) ...
									); 
		else
			opts.variable_window = ceil((opts.Bsize/size(blockmat,2))/3);
			opts.generatorfun = @()( ...
									variableblocksampler(opts.B, ...
												opts.Bsize,  ...
												blockmat, ...
												opts.variable_window) ...
									); 
		end							
	else
		opts.generatorfun = @()( ...
								gensubsamples(opts.B, ...
											opts.n, ...
											opts.Bsize) ...
								); 
	end
	
	[samples opts] = generator.generate_samples(opts); 
	
	resampler.options = opts;
	resampler.samples = samples;
	resampler.iterator = [];
	
end