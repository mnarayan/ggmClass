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
options.poolCovariance = @(Si,ni)(sum(Si.*repmat(ni,[size(Si,1) size(Si,2)]),3)/sum(ni)); 
options.projectCovariance = @projection_tangentspace;

%  Fixed point algorithm: 
lr_dim = 1;
options.lr_dim = lr_dim;
% Initialize
Si = []; ni=[];
if(strcmp(options.mode,'X'))
    for cc=1:s
        Si(:,:,cc) = options.sampleCovariance(A(:,:,cc)); 
        ni(cc) = size(A,1); 
    end
else
    Si = A;
    % What is ni in this case? 
    ni = 1/size(Si,3); 
end
clear A; 

hatSigma = options.poolCovariance(Si,ni); 
[V D] = eig(hatSigma); 
[~ lr_idx] = sort(D,'descend'); 
B = V(:,lr_idx(1:lr_dim)) * diag(sqrt(D(lr_idx(1:lr_dim)))); % Shared low rank


options.projectfun = @(S)(B*inv(B'*S*B)*B'*S); 

for cc=1:s
    hatSigma_i(:,:,cc) = options.projectCovariance(...
                                    projectfun,  ...
                                    S(:,:,cc),  ...
                                    hatSigma);
end


