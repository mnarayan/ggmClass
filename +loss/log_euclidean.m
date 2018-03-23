function L = log_euclidean(S1,S2)
    % LOSS.LOG_EUCLIDEAN
    % Loss = || log(S1) - log(S2) ||_{frobenius}
    %
    % INPUTS
    %   - S1 : p x p PSD matrix (population)
    %   - S2 : p x p PSD matrix (estimator)
    %
    L = norm(logm(S1)-logm(S2),'fro')
end