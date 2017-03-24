function [Rho] = rank_sample_covariance(X,varargin)
%RANK_SAMPLE_COVARIANCE estimates modified spearman rank correlation as a nonparametric replacement for the sample covariance from standard MLE. 
% 
% REFERENCES: 
% "Regularized Rank Estimation of High Dimensional Nonparanormal Graphical Models", Xue an Zhou (2013)
% 
% 
% SEE ALSO GGM.MLECOVESTIMATE GGM.RANKCOVESTIMATE
% 

	% % Use spearman's
	% Xranks = tiedrank(X,1);
	% InitRho = corr(Xranks);
	% Rho  = 2*sin(pi/6.*InitRho);
	
	% Use Kendall's
	if(exist('kendalltau') & numel(X)<1e6)
		Rho = kendalltau(X); 
	else
		Rho = corr(X,'type','kendall'); 
	end
	
	Rho = sin((pi/2).*Rho); 	
	
end