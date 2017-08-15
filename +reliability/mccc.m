function [mccc Mccc] =  MCCC(X,Y,varargin)
%MCCC - Multivariate Concordance Correlation Coefficient	
% 
% USAGE: MCCC
% 
% INPUTS:
% 	- X : A data matrix of n i.i.d measurements and p features
% 	- Y : A data matrix of n i.i.d measurements and p features
% 	- options (optional)	
% 
% REFERENCES: 
% Hiriote and Chinchilli, "Matrix-based concordance correlation coefficient for repeated measures", Biometrics. (2011)	
% 
% Copyright, 2017, Manjari Narayan
% GPL-2 License
% 

	assert(size(X,1)==size(Y,1),'X and Y must have same number of rows'); 
	assert(size(X,2)==size(Y,2),'X and Y must have same number of cols'); 

    crosscov = @(X,Y)(X'*Y/size(X,1));

    [mccc Mccc] = matrix_ccc(X,Y,@cov,crosscov); 
    

	function [mccc M_ccc] = matrix_ccc(X,Y,covfun,crosscov)
        
        n = size(X,1); p = size(X,2);
        
		mu_x = mean(X,1)*n; 
		mu_y = mean(Y,1)*n;

		X = bsxfun(@minus,X,mu_x/n); 
		Y = bsxfun(@minus,Y,mu_y/n);
        		
		V_ind = covfun(X) + covfun(Y) + (mu_x - mu_y)*(mu_x - mu_y)'/(n*(n-1));
        V_dep = V_ind - crosscov(X,Y) - crosscov(Y,X);
        
        [Veig Deig] = eig(V_ind);
        nonneg_eig = find(diag(Deig)>0);
        Vposdef = Veig(:,nonneg_eig);
        Dposdef = Deig(nonneg_eig,nonneg_eig);
        Vsqinv = Vposdef*diag(sqrt(1./diag(Dposdef)))*Vposdef';
        
        M_rho = Vsqinv * V_dep * Vsqinv;
		M_ccc = eye(p) - M_rho;
        
        mnorm = @(A)(norm(A,'fro')); 
        mccc = 1 - mnorm(M_rho)/mnorm(eye(p));
        
	end

	
end