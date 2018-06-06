function [ba varargout] =  Bland_Altman_LOA(M,D,varargin)
%Bland_Altman_LOA - Limits of Agreement for Bland Altman Difference/Mean Statistics

% Currently this implementation assumed paired data where each row of X, Y corresponds to the same unit of measurement i.e. subject. 
% 
%   d_ij = x_ij - y_ij
%         = mu_d + a_i + e_ij
%        where, 
%        i = 1, 2, ...,  n  subjects (rows, independent units of measurement)
%        j = 1, 2, ...., m  replicates or sessions (any repeated measurements)
%        k = 1, 2, ...., p  features
% 
% USAGE:
%       % Given X and Y which are two methods or two raters
%       [ba] = bland_altman_MD(X,Y); 
%       % Given D and M for features x subjects per session
%       D = ba.D; M = ba.M;
%       % Find limits of agreement per feature
%       bastats = bland_altman_loa(M(feature,:)',D(feature,:)'); 
% 
% INPUTS:
% 	- M : Mean data matrix of n i.i.d measurements and m sessions (Method 1)
% 	- D : Difference data matrix of n i.i.d measurements and m sessions (Method 2)
% 	- options (optiona)
% 
% OUTPUTS:
%   - ba.D : Difference matrix X - Y
%   - ba.M : Mean matrix (X + Y)/2
%   - ba.mode : 'raw' or 'percent' 
%   - ba.stats.mu : Estimates of mu_d
%   - ba.stats.loa : Agreement interval for mu_d. (:,1) is lower and (:,2) is upper
%   - ba.stats.alpha : Tail cutoff for agreement interval. alpha = 5% by default. 

% REFERENCES: 
% Bland and Altman 1983, The Statistician
% Bland and Altman 1986, The Lancet
% Bland and Altman 1999, Statistical Methods in Medical Research
% Zou, 2013
% 
% Copyright, 2018, Manjari Narayan
% GPL-2 License
% 

    ba = {}; 

    ba.D = D; 
    ba.M = M;
    ba.stats.mu = bland_altman_stats(M,D); 
    
    rng('default');
    nresamples = 100;

    [ci_boot bstats] = bootci(nresamples,{@bland_altman_stats,M,D}, ...
                                         'type', 'bca');
    ba.stats.ci = ci_boot;
    ba.stats.bstats = bstats;

end

function [mu sd] = bland_altman_stats(M,D)
    
    mu = mean(D); %./std(M);
    sd =  std(D);
    
end



