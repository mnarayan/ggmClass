function [Sigma results] = conditional_sample_covariance_separate(X,options)
%CONDITIONAL_SAMPLE_COVARIANCE_SEPARATE
% 
% Schur complement using correlation decomposition
% Schur complement using covariance decomposition
% Residuals after regressing out confounds and modified cholesky decomposition 
% 	
% USAGE: 
% 
% INPUT
% 	- X : A single n x p Data matrix
% 	- options.scaled: Columns of X scaled to unit diagonals
% 	- options.verbose: Plot and print outputs
% 	- options.outputdir: Defaults to 'tmp/'
% 	- options.filename: Defaults to 'conditionalCov'
	
	
	if(~exist('options','var'))
		options = []
	end
	
	if(isempty(options))
		options.verbose = false;
		options.outputdir = fullfile('tmp',datestr(now,'dd-mmm-yyyy-HHMM')); 
		options.filename = mfilename;
	end
	% TBD: Add process_options or processArgs
	% ~isfield(options,'verbose')
	
	if(~isfield(options,'corrfun'))
		options.corrfun = @(X)(covariance.mle_sample_covariance( ...
									X, ...
									struct('standardize',false) ...
									));
	end
	
	function Sig = my_xcorr(X,y,options)
	% Returns p x q cross correlation between n x p and n x q	
	
		[n1 p] = size(X); 
		[n2 q] = size(y); 			
		assert(n1==n2,'X and y need to have same row length')
		assert(n2>q,'Too many nuisance variables'); 
	
		if(~exist('options','var'))
			options = [];
		end
		if(~isfield(options,'standardize'))	
			options.standardize = false;
		end
		
		if(options.standardize)
			Xt = standardize.successive_normalize(X');
			X  = Xt'; 
			if(q>1)
				y = standardize.successive_normalize(y); 	
			else
				y = standardize.standardize_cols(y); 	
			end
			Sig = X'*y/n1;
					
		else
			[~, mu_x, sig_x] = standardize.standardize_cols(X); 
			[~, mu_y, sig_y] = standardize.standardize_cols(y); 
			%y = zscore(y);
			
			if(q>1)
				Sig = bsxfun(@minus,X,mu_x)'*bsxfun(@minus,y,mu_y)/n1; 
			elseif(q==1)
				Sig = bsxfun(@minus,X,mu_x)'*(y-mu_y)/n1; 
			end
			
			Sig = diag(1./sig_x)*Sig*diag(1./sig_y);			
			
		end
		
		
	end
	
	if(~isfield(options,'xcorrfun'))
		
		options.xcorrfun = @(X,y)(my_xcorr(X,y,struct('standardize',false)));
		
	end
	
	if(~isfield(options,'nuisance'))
		options.nuisance = mean(X,2); 
	end
	
	y = options.nuisance;
	Sig_xx = options.corrfun(X); 		
	Sig_xy = options.xcorrfun(X,y); 
	Sig_yy = options.corrfun(y);
	SigmaL = Sig_xy*inv(Sig_yy)*Sig_xy';
	Sigma = Sig_xx - SigmaL; 
	
	results.SigXX = Sig_xx;
	results.SigXY = Sig_xy;
	results.SigYY = Sig_yy;
	results.nCov  = SigmaL;
	results.Sigma = Sigma;
	results.NSR   = trace(SigmaL)/trace(Sig_xx);
	results.Y = y;
	
	function Xy =  Xy_orthogonalize(X,Y)
	% Nuisance signal regression
		
		% p = size(X,2);
		% Xy = zeros(size(X));
		%
		% for ii=1:p
		% 	gsProj = gsr'*X(:,ii)/(gsr'*gsr)*gsr;
		% 	Xy = X(:,ii) - gsProj;
		% end
		proj = @(X)(X*pinv(X'*X)*X');
		Xy = X - proj(Y)*X;
	end

	results.X_perpY = Xy_orthogonalize(X,options.nuisance);
	
	
	if(options.verbose)
		% Only perform for comparison purposes;
		gsr = options.nuisance; 				
		Xgsr = Xy_orthogonalize(X,gsr); 
		
		figure; 
		subplot(2,2,1); 
		imagesc(Sig_xx); colorbar; axis equal image; 
		title('Correlation'); 
		subplot(2,2,2); 
		imagesc(Sigma); colorbar; axis equal image; 
		title('Conditional Corr.'); 
		subplot(2,2,3);  
		imagesc(SigmaL); colorbar; axis equal image; 
		title('Nuisance Factor'); 
		subplot(2,2,4); 
		imagesc(options.corrfun(Xgsr)); colorbar; axis equal image; 
		title('Correlation (NSR)');
		print('-dpng','-r300',fullfile(options.outputdir,options.filename)); 			
		close all;
	end	
	
	
	
end