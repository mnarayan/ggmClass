function [mccc Mccc varargout] =  MCCC(X,Y,varargin)
%MCCC - Multivariate Concordance Correlation Coefficient	
% 
% The matrix_ccc statistic is 
%       I - {V_indep}^-1/2 * V_dep * {V_indep}^-1/2
% and the corresponding coefficient is 
%       1 - matrixnorm(matrix_ccc)/matrixnorm(I_p)
% where,  
%   V_dep = sum_i (xi-yi)'*(xi-yi)
%   V_indep = sum_i sum_j (xi-yj)'(xi-yj)
% 
% USAGE: [mccc] = MCCC(X,Y);
%        [mccc, Mccc,  mccc_ci] = MCCC(X,Y);
%        [mccc, Mccc, mccc_ci, mccc_boot] = MCCC(X,Y);
% 
% INPUTS:
% 	- X : A data matrix of n i.i.d measurements and p features
% 	- Y : A data matrix of n i.i.d measurements and p features
% 	- options (optional)	
% 
% REFERENCES: 
% Hiriote and Chinchilli, "Matrix-based concordance correlation coefficient for repeated measures", Biometrics. (2011)	
% 
% Copyright, 2017, Manjari Narayan
% GPL-2 License
% 

	assert(size(X,1)==size(Y,1),'X and Y must have same number of rows'); 
	assert(size(X,2)==size(Y,2),'X and Y must have same number of cols'); 

    n = size(X,1); p = size(X,2);

    covfun = @cov;
    pair_crosscov = @(X,Y)(pairwise_crosscov(X,Y));

    [mccc M_ccc precision accuracy normVind normVdep] = ...
                     matrix_ccc(X,Y,covfun,pair_crosscov); 
	
    Mccc = {};
    Mccc.Mccc = M_ccc;
    Mccc.precision = precision;
    Mccc.accuracy = accuracy;
    Mccc.normVind = normVind;
    Mccc.normVdep = normVdep;
    
    s = RandStream('mt19937ar','Seed',1243);
    RandStream.setGlobalStream(s);
    currstate = rng;  
    
    if(nargout==3)    
        [mccc_ci] = matrix_ccc_ci(X,Y,covfun,pair_crosscov);
        varargout{1} = mccc_ci;
    elseif(nargout>=4)
        [mccc_ci mccc_boot] = matrix_ccc_ci(X,Y,covfun,pair_crosscov);
        varargout{1} = mccc_ci;
        varargout{2} = mccc_boot;
    end
    
end

function SigmaPair = pairwise_crosscov(X,Y)
    
    [n p] = size(X);
    
    SigmaPair = zeros(p,p);
    for ii=1:n
        for jj=setdiff(1:n,ii)
            SigmaPair = SigmaPair + X(ii,:)'*Y(jj,:);
        end
    end
    
    SigmaPair = .5*SigmaPair/(size(X,1)*(size(Y,1)-1));
end


function [mccc varargout] = matrix_ccc(X,Y,covfun,pair_crosscov)
    
    [n p] = size(X);
    verbose = false;
    crosscov = @(X,Y)(X'*Y/(size(X,1)-1));

    mu_x = mean(X,1);
    mu_y = mean(Y,1);

    Xc = bsxfun(@minus,X,mu_x);
    Yc = bsxfun(@minus,Y,mu_y);
    
    mean_agreement = (mu_x - mu_y)'*(mu_x - mu_y);
    
	V_ind = covfun(Xc) + covfun(Yc) - ...
                        .5*pair_crosscov(Xc,Yc) - ...
                        .5*pair_crosscov(Yc,Xc) ...
                        + mean_agreement;
    V_dep = covfun(Xc) + covfun(Yc) ...
                        - crosscov(Xc,Yc) - crosscov(Yc,Xc) ...
                        + mean_agreement;
    
    %mnorm = @(A)(norm(A,'fro'));
    mnorm = @(A)(trace(A));
    %mnorm = @(A)(norm(A,2));
    if(nargout>=2)
        tol = 1e-5;
        [Veig Deig] = eig(V_ind);
        nonneg_eig = find(diag(Deig)>tol);
        if(length(nonneg_eig)<p)
            Deig = Deig+eye(p)*abs(min(diag(Deig)));
            nonneg_eig = find(diag(Deig)>tol);
        end
        Vposdef = Veig(:,nonneg_eig);
        Dposdef = Deig(nonneg_eig,nonneg_eig);
        Vsqinv = Vposdef*diag(sqrt(1./diag(Dposdef)))*Vposdef';
    
        M_rho = Vsqinv * V_dep * Vsqinv;
    	M_ccc = eye(p) - M_rho;
        precision = mnorm(V_ind - V_dep)/2; % mccc = precision*accuracy
        accuracy = 2/mnorm(V_ind);
        varargout{1} = M_ccc;
        varargout{2} = precision;
        varargout{3} = accuracy;
        varargout{4} = mnorm(V_ind);
        varargout{5} = mnorm(V_dep);
    else
        nonneg_eig = eigs(V_ind);
        M_rho = [];
    end
    

    if(length(nonneg_eig)==p)
        mccc = 1 - mnorm(V_dep)/mnorm(V_ind);
        %mccc = 1 - mnorm(M_rho);
    else
        if(verbose)
            warning('V_ind is not postive definite');
            disp(['Low rank, df=' num2str(length(nonneg_eig))]);
        end
        mccc = 1 - mnorm(V_dep)/mnorm(V_ind);
        %mccc = 1 - mnorm(M_rho);
        
    end
end

% Generate confidence intervals    
function [mccc_ci mccc_boot] = matrix_ccc_ci(X,Y,covfun,pair_crosscov)
    
    p = size(X,2);
    
    function [mccc] = bootfun(AllData)
        XX = AllData(:,1:p); 
        YY = AllData(:,p+1:2*p); 
        mccc = matrix_ccc(XX,YY,covfun,pair_crosscov);
    end
    
    function [mccc_ci] = bootcifun(Xboot,alphalevel)
        % Xboot is a vector of bootstrapped samples
    
        nboot = length(Xboot);
        bootsort = sort(Xboot);
        low_idx = round(nboot*(alphalevel)/2) + 1;
        hi_idx = nboot - (low_idx-1);
        mccc_ci = nan(2,1);
        mccc_ci(1) = bootsort(low_idx);
        mccc_ci(2) = bootsort(hi_idx);
        
    end
    
    function [jboot] = jafterboot(bootsam)
        
        
        
    end
    
    function [jackci] = jacknifeci(Xjack)
       
        
        
    end
    
    % [mccc_ci mccc_boot] = bootci(50,{@bootfun,cat(2,X,Y)},'type','per');
    
    nboot = 50;
    alphalevel = .05;
    % bootstrap samples
    [mccc_boot mccc_bsamp] = bootstrp(nboot,@(x)(bootfun(x)),cat(2,X,Y));    
    mccc_ci = bootcifun(mccc_boot,alphalevel);
    
end