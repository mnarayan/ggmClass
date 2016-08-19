function [self grphs] = stableEstimate(self,varargin)
	% Creates an ensemble of graph estimaties using a variety of methods
	% Resulting ensemble/model Averaged estimate then stored self.Pi
	% 
	% 
	
	switch nargin
		
	case 2
		currstate = varargin{1};
		modselfun = @stars;
		N = 100;
	case 3
		currstate = varargin{1};
		modselfun = varargin{2}; 
		N = 100;
		
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
	grphs = cell(1,N); 
	
	% random generator state;
	rng(currstate);
	if(m<p) % subsample
		b = max(ceil(10*sqrt(m)), ceil(m/2));
		submat = gensubsamples(N,n,b)
		
	elseif(m>p) % bootstrap
		b = ceil(.7*m);		
		submat = genbootstraps(N,n,b);				
	end
	
	
	if(isa(modselfun,'function_handle'))
		
		switch func2str(modselfun)

		case 'stars'					
			assert(size(submat,1)>=5,'Atleast 50 resamples must be used')
			assert(size(submat,2)>=floor(m/2),'Subsample size should be at least 10\sqrt(n) or m/2')
			grphs = resample_helper(self,submat); 
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
	
		case 'rsample'					
			assert(size(submat,1)>=50,'Atleast 50 resamples must be used')
			assert(size(submat,2)>=floor(m/2),'Subsample size should be at least 10\sqrt(n) or m/2')
			grphs = resample_helper(self,X,submat); 
			tmpPi = zeros(p,p,size(grphs{1},3)); 
			for ii=1:length(grphs)
				tmpPi = tmpPi + (abs(grphs{ii}(:,:,:,1))~=0); 
			end
			self.Pi2 = tmpPi; 
			
		otherwise 					
			[GrphPath] = self.Glasso_Helper(X,0,n,0,0); 
			self.ThetaPath2 = GrphPath;
			self.Theta2 = GrphPath(:,:,end-10);
			warning(sprintf('Model Selection modselfun=%s not implemented',func2str(modselfun))); 
			disp('Returning a path of graphs without model selection ... ')							
		end
	
	else
							
		[GrphPath] = self.Glasso_Helper(X,0,n,0,0); 
		self.ThetaPath2 = GrphPath;				
		self.Theta2 = GrphPath(:,:,end-10);
		warning('Model Selection modselfun is none'); 
		disp('Returning a path of graphs without model selection ... ')	
						
	end
			
end


function [grphs] = resample_helper(self,submat)
	
	[nresamples] = size(submat,1); 
	grphs = cell(1,nresamples); 
	
	parfor ii=1:nresamples
		
		if(self.Lambda==0)
			tmpGGM = GGM(self.Data(submat(ii,:),:,:), 1, 0);
			tmpGGM.estimate();
			tmpGGM
			grphs{ii} = tmpGGM.ThetaPath;
		else
			tmpGGM = GGM(self.Data(submat(ii,:),:,:));
			tmpGGM.estimate();
			tmpGGM
			grphs{ii} = tmpGGM.Theta;
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