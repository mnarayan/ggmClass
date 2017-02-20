function W = adaptiveWeights(Theta,varargin)
	
	narginchk(1,2); 
	p = length(Theta);
	
	if(nargin==3)
		Lambda = varargin{1}; 
	else
		Lambda = .01*ones(p,p);
	end
	nonzeroW = 1.*(abs(Theta)>1e-5); 
	
	A = diag(diag(Theta)) - Theta; 	
	t_start = cputime;
	if(exist('mexeig'))
		[V,D] = mexeig(A);		
	else 
		[V,D] = eig(A);
		D = diag(D);  
	end	
	expA = V*diag(exp(D))*inv(V);			 
	timing = cputime-t_start;
	disp(sprintf('expm took %.4f sec',timing)); 
	
	% % Convert expA into a distance metric
	% diagA = diag(expA);
	% % A_i,i - A_i,j
	% Gpos = bsxfun(@minus,diagA,expA);
	% % A_j,j - A_j,i
	% Gneg = bsxfun(@minus,diagA',expA');
	% distA = (Gpos + Gneg)/2;
	% sqrtD = 1./(diagA);
	% distA = scale_cols(distA,sqrtD);
	% distA = scale_rows(distA,sqrtD);
	
	% distA = 1./abs(expA);		
	% W = ((nonzeroW==1) + (nonzeroW==0).*(distA)).*Lambda;  
	
	W = (nonzeroW==0).*Lambda;
	W = setdiagLS(W,0); 
	
end