function [Sigma Sigma_i] = common_lowrank_decomposition(A,options)
%Common_LOWRANK_DECOMPOSITION decomposes a covariance matrix into two covariances such that the group covariance belongs to common subspace across all covariance and the separate covariance is orthogonal to the common subspace. 
% 
% Model: 
%   Sigma_i = \tilde{Sigma}_i + Sigma , i = 1,...,s
%                   = (A_i)^T Lambda_i (A_i) + (A_0)^T Lambda_0 (A_0)
%              where (A_i, A_0) form a p x p orthogonal basis
% 
% Algorithm: Implements a fixed point algorithm that alternates between pooling and projection operations until convergence. 
%   Sigma_i = Sigma + proj(alpha(Sigma))^T (Sigma_i - Sigma) proj(alpha(Sigma))
% 
	
% USAGE: 
% 
% INPUT
% 	- A : consists of multiple data matrices of n x p x s or multiple observed covariance or correlation matrices
%   - 
% 	- options.scaled: Columns of X (dimension 2) scaled to unit diagonals
%   - options.mode: 'X' if input A consists of data matrices (default) or 'S' if sample covariances. 
%   - options.sampleCovariance : Specify function for to compute sample covariance if input is a data matrix
%   - options.poolCovariance: Function for pooling covariance. Default is ('pooledMLE') but could also be geometric mean. 
%   - options.projectCovariance: Function for projecting given covariance onto subspace orthogonal to the shared subspace. 
% 	- options.verbose: Plot and print outputs
% 	- options.outputdir: Defaults to 'tmp/'
% 	- options.filename: Defaults to 'clr_decomposition'



% Initialize default options   
options.scaled = true;
options.mode = 'X'; 
options.sampleCovariance = @(X)(covariance.mle_sample_covariance(X)); 
options.poolCovariance = @(Si,ni)(sum(Si.*repmat(ni,[size(Si,1) size(Si,2) 1]),3)/sum(ni)); 
options.projectCovariance = @covariance.projection_tangentspace;
options.max_iter = 100;
options.err_tol = 1e-5;


%  Fixed point algorithm: 
lr_dim = 1;
options.lr_dim = lr_dim;
% Initialize
Si = []; ni=[];
if(strcmp(options.mode,'X'))
    [n p s] = size(A);
    for cc=1:s
        Si(:,:,cc) = options.sampleCovariance(A(:,:,cc)); 
        ni(1,1,cc) = size(A,1); 
    end
else
    [~, p, s] = size(A);
    Si = A;
    % What is ni in this case? 
    ni = 1/size(Si,3).*ones(1,1,cc); 
end

%clear A; 

kk=0;
conv_err = Inf;
oldSig = options.poolCovariance(Si,ni)- eye(size(Si,2)); 
oldSig_i = Si; 
newSi = zeros(size(Si)); 

while((kk<options.max_iter) && (conv_err>=options.err_tol))
    
    hatSigma = oldSig + options.poolCovariance(newSi,ni); 
    [V D] = eig(hatSigma); 
    [~,lr_idx] = sort(D,'descend'); 
    B = V(:,lr_idx(1:lr_dim)) * diag(sqrt(D(lr_idx(1:lr_dim)))); % Shared low rank

    options.projectfun = @(S)(B*inv(B'*S*B)*B'*S); 
    for cc=1:s
        [hatSigma_i(:,:,cc) ResProjSigma] = options.projectCovariance(...
                                        options.projectfun,  ...
                                        oldSig_i(:,:,cc),  ...
                                        hatSigma);
    end
    
    % Check convergence    
    err_Sig = norm(hatSigma-oldSig,'fro')/numel(oldSig);
    err_Sigi = sum(sum(sum((hatSigma_i-oldSig_i).^2,3)))/numel(oldSig_i); 
     
    if(kk<=5 | mod(kk,5)==0)
        disp(sprintf('Iter:%d, err:%.5f',kk+1,err_Sig)); 
    end    
    conv_err = err_Sig;
    
    % Update Individual Covariances;
    %projS = options.projectfun(hatSigma); 
    %newSi = Si-repmat(projS'*hatSigma*projS, [1 1 s]) + hatSigma_i;
    newSi = Si-hatSigma_i; %+ repmat(hatSigma,[1 1 s]);
    oldSig_i = hatSigma_i;
    oldSig = hatSigma;
    kk=kk+1;

end

Sigma = hatSigma;
Sigma_i = hatSigma_i;