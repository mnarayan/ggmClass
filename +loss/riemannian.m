function L = riemannian(S1,S2,varargin)
% LOSS.RIEMANNIAN
% Equivalent to the Fisher Information Metric
%
%  Loss = || log ((S1)^{-1/2} S2 (S1)^{-1/2}) ||_{frobenius}
%   
% SEE ALSO LOSS.STEIN, LOSS.BHATTACHARYA 
%
% INPUTS
%   - S1 : p x p PSD matrix (population)
%   - S2 : p x p PSD matrix (estimator)

    
    S1inv = inv(S1);
    S1sqinv = sqrtm(S1inv);
    
    Fisher = S1sqinv * S2 * S1sqinv;
    L = norm(logm(Fisher));


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