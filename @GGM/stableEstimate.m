function [self grphs] = stableEstimate(self,varargin)
	% Creates an ensemble of graph estimaties using a variety of methods
	% Resulting ensemble/model Averaged estimate then stored self.Pi
	% 
	% 
	
	switch nargin
		
	case 2
		currstate = varargin{1};
		modselfun = @stars;
		N = 10;
	case 3
		currstate = varargin{1};
		modselfun = varargin{2}; 
		N = 10;
		
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
	if(m<p) % subsample
		b = max(ceil(10*sqrt(m)), ceil(m/2));
		submat = gensubsamples(N,m,b);
		
	elseif(m>p) % bootstrap
		b = ceil(.7*m);		
		submat = genbootstraps(N,m,b);				
	end
	
	
	if(isa(modselfun,'function_handle'))
		
		switch func2str(modselfun)

		case 'stars'					
			assert(size(submat,1)>=2,'Atleast 2 resamples must be used')
			assert(size(submat,2)>=floor(m/2),'Subsample size should be at least 10\sqrt(n) or m/2')
			grphs = resample_helper(self,submat,1); 
			%[score ] = stars(self,...
			%				struct('n',m,'p',p,'s',1,'estMethod','QUIC','dm',X,'nlambda',100,'beta',.1),...
			%				grphs); 
			% if(self.verbose)
			% 	self.plot(score)
			% end			
			warning('Change to take in alternative grid of lambdas')		
							
		case 'stsel'
			assert(size(submat,1)>=50,'Atleast 50 resamples must be used')
			assert(size(submat,2)>=floor(m/2),'Subsample size should be at least 10\sqrt(n) or m/2')
	
		case 'unionint'
			disp('Union of Intersection Double Resampling not Implemented')
	
		case 'resample'					
			assert(size(submat,1)>=10,'Atleast 10 resamples must be used')
			assert(size(submat,2)>=floor(m/2),'Subsample size should be at least 10\sqrt(n) or m/2')
			grphs = resample_helper(self,submat); 
			tmpPi = zeros(p,p,size(grphs{1},3)); 
			for ii=1:length(grphs)
				tmpPi = tmpPi + (abs(grphs{ii,1}(:,:,:,1))~=0); 
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
			
end


function [grphs] = resample_helper(self,submat,varargin)
	
	[nresamples] = size(submat,1); 
	grphs = cell(nresamples,2); 		
	
	if(nargin>=3)
		usePairs = varargin{1}; 
	else
		usePairs = false;
	end

	m = size(self.Data,1);
	
	if(usePairs)
		for ii=1:nresamples
			warning('Will be slow. No parfor being used')
			tmpgrphs = cell(1,2);
			if(self.Lambda==0)
				tmpGGM = GGM(self.Data(submat(ii,:),:,:), 1, 0);
				tmpGGM.estimate();
				if(ii==1)
					tmpGGM
				end			
				tmpgrphs{1,1} = tmpGGM.ThetaPath;
				clear tmpGGM
				cp_submat = setdiff([1:m],submat(ii,:));
				tmpGGM = GGM(self.Data(cp_submat,:,:), 1, 0);
				tmpGGM.estimate();
				tmpgrphs{1,2} = tmpGGM.ThetaPath;
			else
				tmpGGM = GGM(self.Data(submat(ii,:),:,:));
				tmpGGM.estimate();
				if(ii==1)
					tmpGGM
				end
				tmpgrphs{1,1} = tmpGGM.Theta;			
				clear tmpGGM
				cp_submat = setdiff([1:m],submat(ii,:));
				tmpGGM = GGM(self.Data(cp_submat,:,:));
				tmpGGM.estimate();
				tmpgrphs{1,2} = tmpGGM.Theta;
			end
			grphs(ii,:) = tmpgrphs;
		end		
	else
		parfor ii=1:nresamples
			tmpgrphs = cell(1,2);
			if(self.Lambda==0)
				tmpGGM = GGM(self.Data(submat(ii,:),:,:), 1, 0);
				tmpGGM.estimate();
				if(ii==1)
					tmpGGM
				end			
				tmpgrphs{1,1} = tmpGGM.ThetaPath;
			else
				tmpGGM = GGM(self.Data(submat(ii,:),:,:));
				tmpGGM.estimate();
				if(ii==1)
					tmpGGM
				end
				tmpgrphs{1,1} = tmpGGM.Theta;			
			end
			grphs(ii,:) = tmpgrphs;
		end	
	end
end

function [submat] = genbootstraps(N,n,b)
	
	submat = zeros(N,b);
	for i = 1:N
		submat(i,:) = randsample(n,b,1);
	end
	
end

function [submat] = gensubsamples(N,n,b)
	
	submat = zeros(N,b);
	for i = 1:N
		submat(i,:) = randsample(n,b);
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
	N2 = floor(n/(b2+1));
	submat2 = zeros([N2 b2]);
	for i = 1:N2
		submat2(i,:) = [max(1,(i-1)*b2+1):min(i*b2,n)];
	end	
	submat = zeros([N N2*b2]);
	
	blklength = floor(b/(b2+1));
	if(replacement)
		blockresamples = genbootstraps(N,N2,blklength);
	else
		blockresamples = gensubsamples(N,N2,blklength);
	end
	for ii = 1:N
		submat(ii,:) = reshape(submat2(blockresamples(ii,:),:)',[1 blklength*b2])';
	end
end