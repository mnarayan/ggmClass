function L = bhattacharya(S1,S2,varargin)
% LOSS.BHATTACHARYA
% a) Root-Stein Loss = ( log | .5 S1 + .5 S2| - .5 log | S1*S2 | )^{1/2}
% 
%                         log | .5 S1 + .5 S2| 
% b) Bhattacharya Loss =  ____________________
%                           .5 log | S1*S2 |
%
% INPUTS
%   - S1 : p x p PSD matrix (population)
%   - S2 : p x p PSD matrix (estimator)
%

    if nargin<3
        methodtype = 'root-stein';
    else
        methodtype = varargin{1};
    end
    
    % Arithmetic Mean
    AM = (S1 + S2)/2;
    GM = S1 * S2;

    switch methodtype
    case {'root','root-stein','rootstein','sdiv','Stein'}
        L = sqrt(fastlogdet(AM) - .5*fastlogdet(GM));
    case {'bhattacharya','Bhattacharya'}
        L = sqrt(fastlogdet(AM)/(.5*fastlogdet(GM)));
    end
    
end



function y =  fastlogdet(A)
    
    try
        U = chol(A);
        y = 2*sum(log(diag(U)));
    catch
        [V D] = eig(A);
        y = 2*sum(log(diag(D)));
    end
    
end