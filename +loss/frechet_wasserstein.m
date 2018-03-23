function L = frechet_wasserstein(S1,S2)
% LOSS.FRECHET_WASSERSTEIN
% 
% Given X ~ N(0,S1) and Y ~ N(0,S2)
% d(X,Y) = [trace(S1) + trace(S2) - 2 trace((S1 * S2)^1/2)]^(1/2)
% 
% Reference
% The Fréchet Distance between Multivariate Normal Distributions 
% JOURNAL OF MULTIVARIATE ANALYSIS 12, 450-455 (1982)
% 
% To estimate the Frechet-Wasserstein barycenter, use the following fixed point
% algorithm 
% https://arxiv.org/pdf/1511.05355.pdf
% 
% Needs to be applied to regularized sample covariance estimators. 

    T1 = trace(S1); 
    T2 = trace(S2); 
    T3 = trace(sqrtm(S1*S2));
    
    L = sqrt(T1 + T2 - 2 * T3);
    
end