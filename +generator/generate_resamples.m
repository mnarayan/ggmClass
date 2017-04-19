function resampler = generate_resamples(varargin)
	
	if(nargin==0)
		opts = generator.create_options(100);
	else
		opts = varargin{1}; 
	end
	
	if(~opts.use_replacement)
		warning('Replacement is off. Setting to true')
		opts.use_replacement = true;
	end
	
	function [submat] = genbootstraps(N,n,b)

		submat = zeros(N,b);
		for ii = 1:N
			submat(ii,:) = sort(randsample(n,b,1),'ascend');
		end

	end
	
	function [submat] = fixedblocksampler(N,Bsize,blockmat);

		[n_blocks blocklength] = size(blockmat); 
		
		% Ensuring submat Bsize requirement is met
		nblks = ceil(Bsize/(blocklength+1));	% # blocks to resample
		submat = zeros([N nblks*blocklength]);	% Ensure Bsize is met
		
		submat = [];
		blockresamples = genbootstraps(N,n_blocks,nblks);
		for ii = 1:N
			tmprs = blockmat(blockresamples(ii,:),:);
			submat(ii,:) = sort(reshape(tmprs,[1 numel(tmprs)]));
		end
		
	end
	
	function [submat] = variableblocksampler(N,Bsize,blockmat,var_window)
		
		% for b = 1:B
		% 	tmp_resamples = [];
		% 	blks_idx = round(linspace(TR_idx(80),TR_idx(end-300),n_blks)) + randsample(200,n_blks,0)';
		% 	for nn = 1:n_blks
		% 		tmp_resamples = [tmp_resamples TR_idx(blks_idx(nn)-blks_win:blks_idx(nn)+blks_win)];
		% 	end
		% 	resamples(b,:) = tmp_resamples;
		% end
		
	    [tot_blocks blocklength] = size(blockmat); 		
		n_blocks = ceil(Bsize/(blocklength+1));	% # blocks to resample
		%submat = zeros([N nblks*blocklength]);	% Ensure Bsize is met
		submat = cell(1,N); 
		
		if(var_window>blocklength)
			warning('Length of variable block component too large'); 
			disp('Setting to max variable length to .25 of blocklength'); 
			var_window = max(3,ceil(blocklength/3));
		end
		
		blockresamples = genbootstraps(N,tot_blocks,n_blocks);
		
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
		opts.generatorfun = 	@()( ...
								genbootstraps(opts.B, ...
											opts.n, ...
											opts.Bsize) ...
								); 
	end

	[samples opts] = generator.generate_samples(opts); 

	resampler.options = opts;
	resampler.samples = samples;
	resampler.iterator = [];

end