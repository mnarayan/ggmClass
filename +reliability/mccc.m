function [mccc Mccc varargout] =  MCCC(X,Y,varargin)
%MCCC - Multivariate Concordance Correlation Coefficient	
% 
% The matrix_ccc statistic is 
%       I - {V_indep}^-1/2 * V_dep * {V_indep}^-1/2
% and the corresponding coefficient is 
%       1 - matrixnorm(matrix_ccc)/matrixnorm(I_p)
% where,  
%   V_dep = sum_i (xi-yi)'*(xi-yi)
%   V_indep = sum_i sum_j (xi-yj)'(xi-yj)
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

    n = size(X,1); p = size(X,2);

    covfun = @cov;
    pair_crosscov = @(X,Y)(pairwise_crosscov(X,Y));

    [mccc Mccc] = matrix_ccc(X,Y,covfun,pair_crosscov); 
	
    if nargout>=3    
        mccc_ci = matrix_ccc_ci(X,Y,covfun,pair_crosscov);
        varargout{1} = mccc_ci;
    end
    
end

function SigmaPair = pairwise_crosscov(X,Y)
    
    [n p] = size(X);
    
    SigmaPair = zeros(p,p);
    for ii=1:n
        for jj=setdiff(1:n,ii)
            SigmaPair = SigmaPair + X(ii,:)'*Y(jj,:);
        end
    end
    
    SigmaPair = .5*SigmaPair/(size(X,1)*(size(Y,1)));
end

function [mccc M_ccc] = matrix_ccc(X,Y,covfun,pair_crosscov)
    
    p = size(X,2);
    
    crosscov = @(X,Y)(X'*Y/(size(X,1)));
            
    mu_x = mean(X,1);
    mu_y = mean(Y,1);

    X = bsxfun(@minus,X,mu_x);
    Y = bsxfun(@minus,Y,mu_y);
    		
	V_ind = covfun(X) + covfun(Y) - ...
                        .5*pair_crosscov(X,Y) - ...
                        .5*pair_crosscov(Y,X); 
                        % + (mu_x - mu_y)*(mu_x - mu_y)'/(n*(n-1));
    V_dep = covfun(X) + covfun(Y) - crosscov(X,Y) - crosscov(Y,X);
    
    [Veig Deig] = eig(V_ind);
    nonneg_eig = find(diag(Deig)>0);
    Vposdef = Veig(:,nonneg_eig);
    Dposdef = Deig(nonneg_eig,nonneg_eig);
    Vsqinv = Vposdef*diag(sqrt(1./diag(Dposdef)))*Vposdef';
    
    M_rho = Vsqinv * V_dep * Vsqinv;
	M_ccc = eye(p) - M_rho;
    
    mnorm = @(A)(trace(A)); % Alternative, @(A)(norm(A,'fro')); 
    mccc = 1 - mnorm(M_rho)/mnorm(eye(p));
    
end

% Generate confidence intervals    
function [mccc_ci mccc_boot] = matrix_ccc_ci(X,Y,covfun,pair_crosscov)
    
    p = size(X,2);
    
    function [mccc] = bootfun(AllData)
        XX = AllData(:,1:p); 
        YY = AllData(:,p+1:2*p); 
        mccc = matrix_ccc(XX,YY,covfun,pair_crosscov);
    end
    
    [mccc_ci mccc_boot] = bootci(100,{@bootfun,cat(2,X,Y)},'type','per');
    
end