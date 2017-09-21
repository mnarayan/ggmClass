function shrinkage =  cnshrinkage(eigcoefs, kappa_max)
    %CNSHRINKAGE performs eigenvalue shrinkage for condition number regularized MLE where largest and smallest eigenvalues of the sample covariance undergo double-sided shrinkage
    % 
    %   eig_shrinkage_i =   min(max(tau*,l_i),kappa_max*tau*)
    %   
    %                       tau*  if l_i ≤ tau*
    %                   =   l_i   if tau* ≤ l_i ≤ kappa_max * tau*
    %                       kappa_max * tau* if l_i ≥ kappa_max * tau* 
    %       
    %   where l_1 > l_2 > .... > l_p are the ordered eigenvalues of Shat
    % 
    % tau* = 1/(kappa_max u*) = 
    %           sum_{i=1 to alpha*}l_i/kappa_max + sum_{i=beta* to p}l_i
    %           --------------------------------------------------------
    %                           alpha* + p - beta* + 1
    % 
    %   where, 
    %    alpha* is largest index in (1,...,p) s.t. l_{alpha*} < tau*
    %    beta* is smallest index in (1,...,p) s.t. l_{beta*} > kappa_max * tau*
    % 
    % kappa_max is the regularization parameter that varies from 1 to condition_number(Shat)
    % 
    % USAGE: 
    %
    % SEE ALSO STEIN_SHRINKAGE_SAMPLE_COVARIANCE 
    % 
    
    shrinkage = {};
    shrinkage.Shat = [];
    shrinkage.eig = eigcoefs;
    shrinkage.eig_shrinkage = {};
    shrinkage.tau = [];
    shrinkage.objective = [];
    
    n_path = length(kappa_max);
    n_eigs = length(eigcoefs);
    p = n_eigs;
    
    % For each eigenvalue
    % Univariate shrinkage operator
    
    objfun = @(eigvalue,mu)(eigvalue.*mu-log(mu)); 
    
    
    for kappano=1:n_path
        tmp_shrinkage = zeros(n_eigs,n_eigs-1);
        tmp_obj =zeros(1,n_eigs-1);
        for eigno=2:n_eigs
            
            tau_hi = eigcoefs((eigno-1))
            tau_lo = tau_hi/kappa_max(kappano)
            
            lo_idx = max(find(eigcoefs>=tau_hi));
            hi_idx = min(find(eigcoefs<=tau_lo));
            lo_idx = max([1 lo_idx]);
            hi_idx = min([n_eigs hi_idx]);
            % disp('Largest eigenvalues');
            % eigcoefs(lo_idx)'
            % disp('Smallest eigenvalues');
            % eigcoefs(hi_idx)'
            
            disp(['alpha*: ' num2str(lo_idx) ...
                    ' beta*: ' num2str(hi_idx)]);
            disp(['tau_hi: ' num2str(eigcoefs(lo_idx)) ...
                    ' tau_lo: ' num2str(eigcoefs(hi_idx))]); 
            
            tau = sum(eigcoefs(1:lo_idx)./kappa_max(kappano)) ...
                         + sum(eigcoefs(hi_idx:end));
            
            tau = tau/(lo_idx + (p + 1 - hi_idx));
            
            tmp_shrinkage(:,eigno-1) = ...
                 eig_shrinkage(eigcoefs,tau,kappa_max(kappano));
                 
            tmp_obj(eigno-1) = sum(objfun(eigcoefs,1./tmp_shrinkage(:,eigno-1)));
            
        end
        
        [fminval fmin] = min(tmp_obj);
        disp(['Min. objective function: ' num2str(fminval)])
        
        % lo_idx = eigno;
        % hi_idx = min(find(eigcoefs>kappa_max(kappano)*eigcoefs(lo_idx)));
        
        % tau_lo = eigcoefs(n_eigs-(fmin-1))
        % tau_hi = tau_lo*kappa_max(kappano)
        
        tau_hi = eigcoefs((fmin));
        tau_lo = tau_hi/kappa_max(kappano);      
        lo_idx = max(find(eigcoefs>=tau_hi));
        hi_idx = min(find(eigcoefs<=tau_lo));
        tau = sum(eigcoefs(1:lo_idx)./kappa_max(kappano)) ...
                     + sum(eigcoefs(hi_idx:end));
        tau = tau/(lo_idx + (p + 1 - hi_idx));
        
        lo_idx = max(find(eigcoefs>=tau*kappa_max(kappano)));
        hi_idx = min(find(eigcoefs<=tau));
        tau_opt = sum(eigcoefs(1:lo_idx)./kappa_max(kappano)) ...
                     + sum(eigcoefs(hi_idx:end));
        tau_opt = tau_opt/(lo_idx + (p + 1 - hi_idx));
        
        shrinkage.eig_shrinkage{kappano} = ...
         eig_shrinkage(eigcoefs,tau_opt,kappa_max(kappano));
        
        shrinkage.tau(kappano) = tau_opt; 
        shrinkage.objective{kappano} = tmp_obj;
    end
    
end


function coefs = eig_shrinkage(eigcoefs,tau,kappa)
    % 
    %
    % 
    disp(sprintf('Range: %0.2f to %0.2f ',tau,kappa*tau))
    
    coefs = eigcoefs;
    
    coefs_ceiling =  eigcoefs>=kappa*tau; 
    coefs_floor = eigcoefs<=tau;
    coefs(find(coefs_ceiling)) = kappa*tau;
    coefs(find(coefs_floor)) = tau;
    
end