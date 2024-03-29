function [self results] = sparseMLE(self,varargin)
	% 
	
	
	[m p n] = size(self.Data);
	SigmaHat = self.Sigma;
	Lambda = self.Lambda;
	

	if(self.verbose)
		SigmaHat(1:5,1:5)
		disp('Maximum off-diagonal in S1')
		max(max(triu(abs(SigmaHat),1)))
		%disp('340')
		pause(.01)	
	end

	if(Lambda~=0)
		
		if(self.verbose)
			disp(sprintf('Lambda = %.2f is provided ... ',Lambda));
		end
		
		if(self.verbose)
			disp('Estimating Theta ...')
			pause(.01)
		end
		
		[Thetahat0 SigmaHat2] 		=  QUIC(...
										'default', 	...
										SigmaHat,  	...
										Lambda,		...
										1e-6,  		...
										0,  		...
										1000);
		self.Theta = Thetahat0;
		results.SigmaHat = SigmaHat2;

		% % Relaxed/Refitted MLE
		% L_w = zeros(p,p);
		% L_w = 1./abs(Thetahat0); L_w = L_w/trace(L_w);
		% L_w(Thetahat0==0) = 1000; L_w(find(eye(p))) = 0;
		%
		% [Thetahat1 SigmaHat2] = QUIC('default', SigmaHat, Lambda*L_w, 1e-6, 0, 1000);
		%
		% self.AdaptTheta = Thetahat1;
		% self.Wt = L_w;
		[self.AdaptTheta results] = GGM.constrainMLE(...
									SigmaHat, ...
									Thetahat0, ...
									1.0); 
		self.Wt = results.Lw;
	
		results.constrainMLE = results;
				
		clear Thetahat0 Thetahat1 SigmaHat2;
		
		if(self.verbose)
			disp('Done ...')
			pause(.01)
		end

	elseif(Lambda==0)
		
		assert(Lambda==0,'Lambda should be 0 in path mode');
		 
		if(self.verbose) 
			disp('Lambda not provided ... ');
			disp('Estimating Full Regularization Path ....')
		end
		
		% minL = min(min(triu(abs(SigmaHat),1)));
		% maxL = max(max(triu(abs(SigmaHat),1)));
		% %Lrange = fliplr(logspace(-3,log10(maxL),100));
		% Lrange = fliplr(linspace(max(minL,1e-3),maxL,100));
		
		assert(isfield(struct(self),'Lrange')==1,'Lrange not defined');
		assert(length(self.Lrange)>5,'Grid of lambda values not defined');
		Lrange = self.Lrange;
		
		if(self.verbose)
			disp('Checking regularization parameter values ...')
			Lrange(1:10)
	%		pause(.01)('392')
		end
	
		if(isreal(Lrange(1)))
			[Thetahat0 SigmaHat2] = QUIC('path',  	...
										SigmaHat,  	...
										1.0,		...
										Lrange,  	...
										1e-6,  		...
										0,  		...
										1000);			
										
			% gr_sparsity = squeeze(sum(sum(Thetahat0~=0,1),2));
			% allzero_idx = min(find((gr_sparsity==0)));
			% self.Lrange = self.Lrange(1:allzero_idx);
			% self.ThetaPath = Thetahat0(:,:,1:allzero_idx);
			self.ThetaPath = Thetahat0;
			self.Theta = [];
			results.sparsity = squeeze(sum(sum(Thetahat0~=0,1),2));
			results.Lrange = Lrange;
			results.SigmaHat = SigmaHat2; 
		else
			error('Regularization values are complex')
			self.ThetaPath = []
			self.Theta = [];
		end
	
	
	end
	
end