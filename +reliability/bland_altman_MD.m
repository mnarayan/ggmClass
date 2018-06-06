function [ba varargout] =  Bland_Altman_MD(X,Y,varargin)
%Bland_Altman_MD - Bland Altman Difference/Mean Statistics
% 
% Currently this implementation assumed paired data where each row of X, Y corresponds to the same unit of measurement i.e. subject. 
% 
%   d_ij = x_ij - y_ij
%        = mu_d + a_i + e_ij
%        where, 
%        i = 1, 2, ..., n  subjects (rows, independent units of measurement)
%        j = 1, 2, ...., m replicates or sessions (any repeated measurements)
% 
% USAGE: Bland_Altman_LOA(X,Y) 
% 
% INPUTS:
% 	- X : A data matrix of n i.i.d measurements and m sessions (Method 1)
% 	- Y : A data matrix of n i.i.d measurements and m sessions (Method 2)
% 	- use_logtransform : Use log(X), log(Y)
%   - use_relative : Use Relative/Percent Difference as in D./M
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

    switch nargin
    case 3
        use_logtransform = varargin{1}; 
        use_relative = false; 
    case 4
        use_logtransform = varargin{1}; 
        use_relative = varargin{2}; 
    otherwise
        use_logtransform = false; 
        use_relative = false; 
    end
    
    ba = {}; 

    if(use_logtransform)
        X = log(X);
        Y = log(Y); 
    end
    
    ba.D = X - Y;
    ba.M = (X + Y)/2;
    
    if(use_relative)
        ba.D = ba.D./ba.M;
    end
    
end

