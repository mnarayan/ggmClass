function [T D] = ls_cholesky(X,varargin)
%ESTIMATOR.LS_CHOLESKY - returns the modified cholesky factor by solving autoregression and innovation errors. 
%     
% where Cov(X)^{-1} = TD^{-2}T' i.e. TX = E such that TCov(X)T' = D^2
%
% Requires the arfit toolbox 
%
% constant lag vs. variable lag 

    X = permute(X,[1 3 2]);
    [m p n] = size(X); 

    fitobj = {};
    pmin = 0;
    pmax = 1;
    T = eye(m); 
    D = zeros(m,m); 

    % stationary AR fitting
    [mu,phi,Sig,~,~,~] = arfit(X,pmin,pmax,'zero');
    lags = length(phi);
    for tt=1:m
        T(tt+1:min(tt+lags,m),tt) = -phi;
    end
    for nn=1:n
        % Using cov(epislon) = cov(Ty);
        D = D + cov(T*X(:,1,nn));
    end
    D = diag(diag(D))/nn;
    L = inv(T); 
    
    % % non-stationary AR fitting
    % for shiftNo=1:m
    %     shift_idx = [ fliplr(shiftNo:m) ]
    %     [mu,phi,Sighat,~,~,~] = arfit(X,pmin,pmax,'zero');
    %     lags = length(phi);
    %     T(shiftNo,1:lags) = -phi;
    %     D(shiftNo,shiftNo) = Sighat(1,1);
    % end
    
    fitobj.params.mu = mu;
    fitobj.params.phi = phi;
    fitobj.params.Sig = Sig;
    fitobj.params.lags = lags;
    fitobj.D = D;   % Prediction variance, innovation variance 
    fitobj.T = T;   % Cholesky factor of inverse covariance
    
end

function [beta] = banded_ols(Y,X)
    
    
    
end


function [k] = banded_BIC(fitobj,varargin)
    
    
    
end


