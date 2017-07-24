function [self results] = stableEstimate(self,varargin)
	% Creates an ensemble of graph estimaties using a variety of methods
	% Resulting ensemble/model Averaged estimate then stored self.Pi
	% 
	% 
	
	results = {};
	
	switch nargin
		
	case 2
		currstate = varargin{1};
		modselfun = @stars;
		N = 25;
	case 3
		currstate = varargin{1};
		modselfun = varargin{2}; 
		N = 25;
		
	case 4
		currstate = varargin{1};
		modselfun = varargin{2}; 
		N = varargin{3};
	end
	
	% Resampling type
	% 	: boot
	% 	: mofnboot
	% 	: subsample
	[m p n] = size(self.Data); 
	grphs = cell(N,2);
	
	% random generator state;
	rng(currstate);
	if(m<5*p) % subsample
		if(m<120)
			b = min(ceil(10*sqrt(m)), ceil(.6*m));
		else
			b = ceil(m/2); 
		end
		%submat = gensubsamples(N,m,b);
        [N,m,b]
		submat = genblockresamples(N,m,max(b,60),0); 
		
	elseif(m>5*p) % bootstrap
		b = ceil(.7*m);		
		submat = genbootstraps(N,m,b);				
	end
	
	
	if(isa(modselfun,'function_handle'))
		
		switch func2str(modselfun)

		case 'stars'					
			assert(size(submat,1)>=2,'Atleast 2 resamples must be used')
			assert(size(submat,2)>=ceil(m/2),'Subsample size should be at least 10\sqrt(n) or m/2')
			grphs = resample_helper(self,submat,1); 
            [scores best_lambda] = estimator.stars(grphs)
            if(self.verbose)
                %self.plot(scores)
            end
			tmpGGM = GGM(self.Data(:,:,1),1,0); 
            [tmpGGM mleresults] = tmpGGM.estimate();
            mleresults.sparsemle.sparsity'
			self.ThetaPath = tmpGGM.ThetaPath;
			self.Theta = tmpGGM.ThetaPath(:,:,best_lambda);			
			warning('Change to take in alternative grid of lambdas')		
							
		case 'stsel'
			assert(size(submat,1)>=50,'Atleast 50 resamples must be used')
			assert(size(submat,2)>=floor(m/2),'Subsample size should be at least 10\sqrt(n) or m/2')
	
		case 'unionint'
			disp('Union of Intersection Double Resampling not Implemented')
	
		case 'resample'					
			assert(size(submat,1)>=2,'Atleast 2 resamples must be used')
			assert(size(submat,2)>=floor(m/2),'Subsample size should be at least 10\sqrt(n) or m/2')
			cp_pairs = true;
			grphs = resample_helper(self,submat,cp_pairs); 
			tmpPi = zeros(p,p,size(grphs{1,1},3)); 
			if(cp_pairs)
				for ii=1:length(grphs)
					tmpPi = tmpPi + (abs(grphs{ii,1}(:,:,:,1))~=0).*(abs(grphs{ii,2}(:,:,:,1))~=0); 
				end
			else
				for ii=1:length(grphs)
					tmpPi = tmpPi + (abs(grphs{ii,1}(:,:,:,1))~=0); 
				end
			end
			self.Pi = tmpPi; 
			
		otherwise 					
			[GrphPath] = self.Glasso_Helper(X,0,n,0,0); 
			self.ThetaPath = GrphPath;
			self.Theta = GrphPath(:,:,end-10);
			warning(sprintf('Model Selection modselfun=%s not implemented',func2str(modselfun))); 
			disp('Returning a path of graphs without model selection ... ')							
		end
	
	else
							
		[GrphPath] = self.Glasso_Helper(X,0,n,0,0); 
		self.ThetaPath = GrphPath;				
		self.Theta = GrphPath(:,:,end-10);
		warning('Model Selection modselfun is none'); 
		disp('Returning a path of graphs without model selection ... ')	
						
	end
	
	results.grphs = grphs; 
	results.modselfun = modselfun;
    if(exist('scores'))
        results.scores = scores;
    end
    if(exist('best_lambda'))
        results.best_lambda = best_lambda;
    end
	results.currstate = currstate; 
	results.submat = submat;
			
end


function [grphs] = resample_helper(self,submat,varargin)
	
	[nresamples] = size(submat,1); 
	
	if(nargin>=3)
		usePairs = varargin{1}; 
	else
		usePairs = false;
	end
	if(length(unique(submat(1,:)))~=size(submat,2))
		usePairs = false;
	end

    if(usePairs)
    	grphs = cell(2*nresamples,1); 		
    else
        grphs = cell(nresamples,1);
    end

	m = size(self.Data,1);


	for ii=1:nresamples
		tmpgrphs = cell(1,1);		
		if(self.Lambda==0)
			tmpGGM = GGM(self.Data(submat(ii,:),:,:), 1, 0);
			tmpGGM.verbose = self.verbose; 
			tmpGGM.estimate();
			if(ii==1)
				tmpGGM
			end			
			tmpgrphs(1) = {tmpGGM.ThetaPath};
		else
			tmpGGM = GGM(self.Data(submat(ii,:),:,:),1,self.Lambda);
			tmpGGM.verbose = self.verbose; 			
			tmpGGM.estimate();
			if(ii==1)
				tmpGGM
			end
			tmpgrphs(1) = {tmpGGM.AdaptTheta};			
		end
		grphs(ii) = tmpgrphs
	end	
	
	if(usePairs)
		cp_submat = []; 
		for ii=1:nresamples
			cp_submat(ii,:) = setdiff([1:m],submat(ii,:));
		end
	end
	
	if(usePairs)
		for ii=1:nresamples
			tmpgrphs = cell(1,1);
			if(self.Lambda==0)
				tmpGGM = GGM(self.Data(cp_submat(ii,:),:,:), 1, 0);
				tmpGGM.verbose = self.verbose; 				
				tmpGGM.estimate();
				if(ii==1)
					tmpGGM
				end			
				tmpgrphs(1) = {tmpGGM.ThetaPath};
			else
				tmpGGM = GGM(self.Data(cp_submat(ii,:),:,:),1,self.Lambda);
				tmpGGM.verbose = self.verbose; 				
				tmpGGM.estimate();
				if(ii==1)
					tmpGGM
				end
				tmpgrphs(1) = {tmpGGM.AdaptTheta};			
			end
			grphs(ii+nresamples) = tmpgrphs;
		end		
    	grphs = reshape(grphs, [nresamples 2]);     
	end
		
end

function [submat] = genbootstraps(N,n,b)
	
	submat = zeros(N,b);
	for ii = 1:N
		submat(ii,:) = randsample(n,b,1);
	end
	
end

function [submat] = gensubsamples(N,n,b)
	
	submat = zeros(N,b);
	for ii = 1:N
		submat(ii,:) = randsample(n,b);
	end
end

function [submat] = genblockresamples(N,n,b,replacement)
% Generate BlockBootstrap Samples
% Inputs
% 	N - Number of Resamples to Generate
%	n - size of the number of available observations
% 	b - size of the bootstrap
	%y = randsample(population,k)
	b2 = ceil(2*n^(1/3));
	%% moving block 
	%N2 = n-b2+1;
	% submat2 = zeros([N2 b2]);
	% for i = 1:N2
	% 	submat2(i,:) = [i:min(i+b2-1,n)];
	% end	
	%% fixed block length
	N2 = floor(n/(b2));
	submat2 = zeros([N2 b2]);
	for ii = 1:N2
		submat2(ii,:) = [max(1,(ii-1)*b2+1):min(ii*b2,n)];
	end	
		
	blklength = ceil(b/(b2-1));
	submat = zeros([N blklength*b2]);
	
	if(replacement)
		blockresamples = genbootstraps(N,N2,blklength);
	else
		blockresamples = gensubsamples(N,N2,blklength);
	end
	for ii = 1:N
		tmprs = submat2(blockresamples(ii,:),:); 
		submat(ii,:) = sort(reshape(tmprs,[1 numel(tmprs)]));
	end
end