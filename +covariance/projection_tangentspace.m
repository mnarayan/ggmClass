function [OrthProjSigma] = projection_tangentspace(projectfun, Si, S); 
%PROJECTION_TANGENTSPACE: Given matrices and their frechet mean, returns projection of each matrix onto tangent direction to frechet mean
% 
% INPUT
%   projectfun  - specific projection matrix
%   Si          - covariance matrix to be projected
%   S           - frechet mean covariance or covariance in some shared subspace
% 
    projS = projectfun(S);
    OrthProjSigma = projS'*(Si-S)*projS; 

end
