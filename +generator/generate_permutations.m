function permuter = generate_permutations(n,varargin)
    %GENERATE_PERMUTATIONS
    % 
    % USAGE: 
    %   permuter = generate_permutations(n)
    %   permuter = generate_permutations(n,options)
    % 
    % INPUTS: 
	%   n - vector of sample sizes. Typically length 1 or length 2
    %   options - output of create_options. (optional) 
      
    
	if(nargin==1)
		opts = create_options();
    elseif(nargin==0)
        warning('At least 1 argument must be provided. Setting n = 10')
        opts = create_options();
        n = 10;
	else
        warning('2nd argument must be previously used options structure')
		opts = varargin{1}; 
	end
    permuter = {};
    permuter.n = n;
    permuter.options = opts;
    permuter.generatorfun = @()(opts.perm_design(n)); 
    disp(sprintf('Generating %d `%s` mode permutations', opts.nperms,opts.design));
    for permno=1:opts.nperms
        permuter.samples{permno} = permuter.generatorfun();
    end
    permuter.options.endstate = rng;
    
	function options = create_options(varargin)
        options = {}; 
        % default options 
        options.seed = 0; 
        options.nperms = 2000; 
        options.currstate = rng('default'); 
        options.perm_dimension = 1;     % permute first dimension of X. 
        options.design = 'signrank';   % onesample, signrank
        [perm_design] = permdesign(options);  % 
        options.perm_design = perm_design;
        options.statistic = ''; %default    
	end
	
    
    function resamples = permute_groups(n)
    %PERMUTE_GROUPS - returns permutation of integers of length n
    % INPUTS
    %   - n is a vector of length 2 n=[n1 n2]; 
    % OUTPUT
    %   - resamples is a cell array. resamples{1} contain permuted samples of length n1 and resamples{2} contains permuted indices of length n2. 
    % 
        ngroups = length(n); 
        if(ngroups>2)
            warning('More than 2 groups not supported. Using first 2 groups');
        elseif(ngroups<=1)
            warning('Less than 2 groups not supported. Use randperm for 1 group'); 
        end
        
        resamples = cell(1,ngroups);
        curr_samples = [];  
        for ii=1:max(1,ngroups-1)
            tmp_n = setdiff(1:sum(n),curr_samples);
            resamples{ii} = randsample(tmp_n,n(ii));
            curr_samples = cat(1,curr_samples,resamples{ii});  
        end
        resamples{ngroups} = setdiff(1:sum(n),curr_samples); 
        
        resamples
    end
    
    
    function perm_design = permute_twosample(n)
        
        ngroups = length(n); 
        if(ngroups>2)
            warning('More than 2 groups not supported. Using first 2 groups');
        end                

        resamples = permute_groups(n(1:2)); 
        
        perm_design = zeros(sum(n(1:2)),2);
        for ii=1:length(resamples)
            perm_design(resamples{ii},ii) = 1; 
        end
        
    end
    
    function perm_design = permute_signrank(n)
        
        if(length(n)>1)
            warning('Too many groups. Only matched pairs for n(1)');
        end                

        %disp('running permute_signrank(n)')
        resamples{1} = randperm(n(1)); 
        
        perm_design = zeros(n(1),2);
        first_half = 1:floor(n(1)/2); 
        sec_half  = floor(n(1)/2)+1:n(1); 
        
        perm_design(resamples{1}(first_half),1) = 1; 
        perm_design(resamples{1}(first_half),2) = -1; 
        perm_design(resamples{1}(sec_half),1) = -1; 
        perm_design(resamples{1}(sec_half),2) = 1; 
                        
    end
    
    
    function [perm_design options]   =   permdesign(options)
        
        switch options.design
            
        case 'twosample'
            perm_design = @(n)(permute_twosample(n)); 
        case 'signrank'
            perm_design = @(n)(permute_signrank(n)); 
        end
        
       options.perm_design = perm_design;
       
    end
    
    
end