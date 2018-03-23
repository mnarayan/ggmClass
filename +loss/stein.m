function L =  stein(S1,S2)
    % LOSS.STEIN
    %   trace(inv(S1)*S2) - logdet(inv(S1)*S2) - p
    % INPUTS
    %   - S1 : p x p PSD matrix (population)
    %   - S2 : p x p PSD matrix (estimator)
   
    S1inv = inv(S1);
    S12 = S1inv*S2;
    p = size(S1,1);
    
    L = trace(S12) - fastlogdet(S12) - p;
    
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