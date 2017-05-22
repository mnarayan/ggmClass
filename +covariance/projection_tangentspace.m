function [OrthProjSigma ResProjSigma] = projection_tangentspace(projectfun, Si, S); 
%PROJECTION_TANGENTSPACE: Given matrices and their frechet mean, returns projection of each matrix onto tangent direction to frechet mean
% 
% INPUT
%   projectfun  - specific projection matrix
%   Si          - covariance matrix to be projected
%   S           - frechet mean covariance or covariance in some shared subspace
% 
    projS = projectfun(S);
    %OrthProjSigma = projS'*(Si-S)*projS; 
    OrthProjSigma = zeros(size(Si)); 
    if(ndims(Si)==3)
        for cc=1:size(Si,3)
            OrthProjSigma(:,:,cc) = Si(:,:,cc)-projS'*(Si(:,:,cc))*projS;
        end
    else
        %disp('Single projection')
        %OrthProjSigma = Si-projS'*(Si)*projS;
        ResProjSigma = projS'*(Si)*projS;
        OrthProjSigma = Si - ResProjSigma;
    end
    
end
