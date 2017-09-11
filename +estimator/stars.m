function [scores best_lambda ustats] = stars(grphs,varargin)	
%ESTIMATOR.STARS returns optional regularization parameter and scores. 
% USAGE: [scores] = stars(grphs)
% INPUTS
% grphs - a cell array of size (nresamples,1) or (nresamples,2). Each cell contains graphs p x p x nlambdas. Must be ordered from sparse to dense.

    if(isempty(grphs))
        error('Input is empty')
    end
    
    [nresamples npairs] = size(grphs);
    
    % Check dimensions of cell array;
    if(npairs==2)
        cp_pairs = true;
    elseif(npairs>2)
        warning(['More than 2 intersections not supported. ...' ...
                'Only using first two columns of grphs cell ...' ...
                ' array as complementary pairs'] ...
                );                
        grphs = grphs{:,1:2}; 
        cp_pairs = true;
    else
        cp_pairs = false;
    end
    
    % Check dimensions of cell data
    switch ndims(grphs{1})
    case 3
        [p,~,nlambdas] = size(grphs{1});
    case 4
        [p,~,nlambdas,nsubjects] = size(grphs{1});
        error('Only p x p x nlambdas supported');
    end   
    
    if(nargin==2)
        options=varargin{1}; 
    elseif(nargin==1)
        options = estimator.create_options(); 
    end
    
    % Get edgelet u-statistics
    ustatfun1 = @(A)(abs(A)~=0);
    ustatfun2 = @(A,B)((abs(A).*abs(B))~=0);

    ustats = zeros(p,p,nlambdas);
    for resampleno=1:nresamples
        if(cp_pairs)
            ustats = ustats + ...
                            (ustatfun2(grphs{resampleno,1},...
                                            grphs{resampleno,2} ...
                                            ));
        else
            ustats = ustats + ...
                            (ustatfun1(grphs{resampleno,1}));
        end
    end
    ustats = ustats/nresamples;
    [scores,best_lambda] = monotonic_score(ustats,options.instability_beta);
    
end



function [mt_score best_lambda] = monotonic_score(ustats,beta,varargin)

    [p,~,nlambdas] = size(ustats);
    
    % 2*pi*(1-pi)
    psistat = 2.*(ustats).*(1-ustats);    
    instability = zeros(1,nlambdas);
    sparsity = instability;
    for lambdano=1:nlambdas
        instability(lambdano) = sum(sum(triu(psistat(:,:,lambdano),1)))...
                                                        /nchoosek(p,2);
        sparsity(lambdano) = sum(sum(triu(ustats(:,:,lambdano)~=0,1)))...
                                                        /nchoosek(p,2);
    end
    
    mt_score = zeros(size(instability)); 
    mt_score(1) = instability(1);
    mt_score(2:end) = max(instability(1:end-1),instability(2:end));
    
    [min_scores,idx_mt_score] = find((mt_score <= beta));
    mt_score2 = zeros(size(mt_score)); 
    mt_score2(idx_mt_score) = min_scores;    
    min_sparsity = sparsity.*(mt_score2>0); % sparsity values at same locations
    min_sparsity(min_sparsity==0) = nan;
    
    [~,best_lambda] = max(min_sparsity(~isnan(min_sparsity))); % who is sparsest non-zero graph
    
end

