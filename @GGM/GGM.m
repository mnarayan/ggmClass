classdef GGM < handle
	% Data structure 
	
	properties 	
		
		% Data
		Data; 
		
		% Options
		verbose = 1;
		useCorr = true;
		randState;
		% Input Parameters
		LagOp;  % Only used for ordered variables
		Lambda; % Regularization parameter for Theta				
		Lrange;
		h; 
		
		% Parameters (To be Estimated)	
		methodfun = @sparseMLE;
		Sigma; % Covariance
		Theta; % Concentration
		W;	   % Variance
		ThetaPath; % Theta estimates along Regularization Path
		AdaptTheta; % Adaptive/Multi-Step Estimates of Theta
		Pi;	% Stability Scores for Theta
		Wt; % Initial Weights for Theta structure 
		
	end
	
	methods


		function self = GGM(varargin)
		
			p = inputParser;
			p.addRequired('self',@(x) true);
			p.addOptional('Data',randn(100,10),@(x) sum(size(x)~=0)>=2);
			p.addOptional('useCorr',true);
			p.addOptional('Lambda',.25,@isVector);
			
			switch nargin
			
			case 0
				p.parse(self);				
				self.Data = p.Results.Data;
				self.useCorr = p.Results.useCorr;
				self.Lambda = p.Results.Lambda;
				
				
			case 1
				
				Data = varargin{1};
				p.parse(self,Data);
				self.Data = Data;
				self.useCorr = p.Results.useCorr;
				self.Lambda = p.Results.Lambda;
				
			case 2
								
				Data = varargin{1};
				p.parse(self,Data);
				self.Data = Data;
				useCorr = varargin{2};
				p.parse(self,useCorr);
				self.useCorr = useCorr;				
				self.Lambda = p.Results.Lambda;
				
			case 3
								
				Data = varargin{1};
				p.parse(self,Data);
				self.Data = Data;
				useCorr = varargin{2};
				p.parse(self,useCorr);
				self.useCorr = useCorr;	
				Lambda = varargin{3};
				p.parse(self,Lambda);			
				self.Lambda = Lambda;
				
			end
			
			% Empirical/Sample Covariance
			[m p n] = size(self.Data);
			Sigma = zeros(p,p);
			for cc=1:n
				Sigma = Sigma + self.Data(:,:,cc)'*self.Data(:,:,cc)/m;
			end
			Sigma = Sigma/n;
			
			if(self.useCorr)
				self.W = diag(diag(Sigma)); 
				tmpD = inv(sqrt(self.W));
				Sigma = tmpD*Sigma*tmpD';
			end
			self.Sigma = Sigma;	
			
			if(self.Lambda==0)
				maxL = max(max(abs(triu(self.Sigma,1))));
				minL = .5*min(min(abs(triu(self.Sigma,1))));
				if(minL==0)
					%minL = .1*sqrt(log(p)/(m*n));
					minL = max([1e-3*maxL,1/(m*n)]);
				end				
				% Lrange = logspace(log10(minL),log10(maxL),100);
				% self.Lrange = fliplr(Lrange);
				% Lrange = logspace(log10(maxL),log10(minL),100);
				Lrange = linspace((maxL),(minL),50);				
				self.Lrange = (Lrange);
			else
				self.Lrange = [];
			end
	
		end
		
		% loadobj
		function A = loadobj(S)
			
			if isstruct(S)
				
				A = GGM();
				A.Data = S.Data;
				A.verbose = S.verbose;
				A.useCorr = S.useCorr;
				A.LagOp = S.LagOp;
				A.Sigma = S.Sigma;
				A.Theta = S.Theta;
				A.W = S.W; 
				A.ThetaPath = S.ThetaPath;
				A.Pi = S.Pi;
				A.Lambda = S.Lambda;
				A.Wt = S.Wt;

			elseif isa(s,'GGM')
			   	A = S;
			else 
			    A = [];
			end
			
		end
		
		% saveobj
		function S = saveobj(A)
			
			S = struct(A);
			
		end
			
	    function newGGM = copy(obj) 
		% function newGGM = copyobj(ggmObj)			
	    % Create a shallow copy of the calling object.
		% Following YAGTOM Manual		
		
	        newGGM = eval(class(obj));
	        meta = eval(['?',class(obj)]);
	        for p = 1: size(meta.Properties,1)
	            pname = meta.Properties{p}.Name;
	            try
	                eval(['newGGM.',pname,' = obj.',pname,';']);
	            catch
	                fprintf(['\nCould not copy ',pname,'.\n']);
	            end
	        end
	    end		
		
		function self = plot(self,varargin)
			% function self = plot(self,varargin)
			% Inputs
			% score (optional): Model selection score across path
			
			switch nargin
				
			case 2
				model_score = varargin{1}
			otherwise
				model_score = [];
			end
			
			minL = min(min(triu(abs(self.Sigma),1))); 
			maxL = max(max(triu(abs(self.Sigma),1))); 
			
			if(length(size(self.Pi))==3)
				% Call plot_path
				h = helper_plot_path(struct('G',self.Pi,'model',self.Theta,'score',model_score,...
					'lambda',self.Lambda,'lambdaList',self.Lrange)); 	
				self.h = h;
			elseif(length(size(self.ThetaPath))==3)
				% Call plot_path
				h = helper_plot_path(struct('G',self.AdaptTheta,'model',self.Theta,'score',model_score,...
					'lambda',self.Lambda,'lambdaList',self.Lrange)); 	
				self.h = h;
			end
	   end
	
	end
	
	methods(Static)
	
		[Theta results] = constrainMLE(SigmaHat,InitialTheta,varargin);
					[W] = adaptiveWeights(Theta,varargin); 
					Rho = rankCovEstimate(X,varargin); 
	
	end
	
end