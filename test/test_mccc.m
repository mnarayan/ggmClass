function test_suite = test_standardize
	buildFunctionHandleTestSuite(localfunctions);

end

function test_null_unit_diagonal()
   
    n = 100; 
    p = 5; 
    ntrials = 100;
    
    Sigma = eye(2*p);
    rng('default');
    
    mccc = nan(1,ntrials);
    
    for trialno=1:ntrials
        AllData = randn(n,2*p)*sqrtm(Sigma);
        X = AllData(:,1:p); 
        Y = AllData(:,p+1:2*p); 
        [mccc(trialno)] = reliability.mccc(X,Y);
    end
    
    disp('Null Hyp. for Unit Diagonal')
    disp(['n=' num2str(n) ' i.i.d measurements, p=' num2str(p) ' variables']);
    disp('Mean mccc (100 trials)')
    mean(mccc)
    disp('MCMC stdev (100 trials)'); 
    std(mccc)
    
end


function test_perfect_unit_diagonal()
   
    n = 100; 
    p = 5; 
    ntrials = 100;
    
    Sigma = eye(p);
    rng('default');
    
    mccc = nan(1,ntrials);
    
    for trialno=1:ntrials
        AllData = randn(n,p)*sqrtm(Sigma);
        X = AllData(:,1:p); 
        Y = AllData(:,1:p);
        [mccc(trialno)] = reliability.mccc(X,Y);
    end
    
    disp('Null Hyp. for Unit Diagonal')
    disp(['n=' num2str(n) ' i.i.d measurements, p=' num2str(p) ' variables']);
    disp('Mean mccc (100 trials)')
    mean(mccc)
    disp('MCMC stdev (100 trials)'); 
    std(mccc)
    
end

function test_alt_unit_diagonal()
   
    n = 100; 
    p = 5; 
    ntrials = 100;
    
    Sigma = eye(2*p)*2*p;
    Sigma(1:p,p+1:2*p) = .9;
    Sigma(p+1:2*p,1:p) = .9;
    rng('default');
    
    mccc = nan(1,ntrials);
    
    for trialno=1:ntrials
        AllData = randn(n,2*p)*sqrtm(Sigma);
        X = AllData(:,1:p); 
        Y = AllData(:,p+1:2*p);
        [mccc(trialno)] = reliability.mccc(X,Y);
    end
    
    disp('Null Hyp. for Unit Diagonal')
    disp(['n=' num2str(n) ' i.i.d measurements, p=' num2str(p) ' variables']);
    disp('Mean mccc (100 trials)')
    mean(mccc)
    disp('MCMC stdev (100 trials)'); 
    std(mccc)
    
end


% function test_hiriote_example()
%
%     n = 100;
%     p = 3;
%     ntrials = 100;
%
%     Sigma = eye(2*p);
%     rng('default');
%
%     mccc = nan(1,ntrials);
%     mccc_ci = nan(2,ntrials);
%
%     for trialno=1:ntrials
%         AllData = randn(n,2*p)*sqrtm(Sigma);
%         X = AllData(:,1:p);
%         Y = AllData(:,p+1:2*p);
%         [mccc(trialno),~,mccc_ci(:,trialno)] = reliability.mccc(X,Y);
%     end
%
%     disp('Null Hyp. for Unit Diagonal')
%     disp(['n=' num2str(n) ' i.i.d measurements, p=' num2str(p) ' variables']);
%     disp('Mean mccc (100 trials)')
%     mean(mccc)
%     disp('MCMC stdev (100 trials)');
%     mean(mccc_ci,2)
%
% end