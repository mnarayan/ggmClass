function W =  kendallsW(X,varargin)
%kendallsW - returns the kendall concordance between features measured by 2 or more raters. 
%
% INPUTS
% 	- X  is a matrix of features x raters or features x subjects
%
% OUTPUTS
%   W is the kendalls W coefficient
	
    switch nargin 
    
    case 1
        % Do nothing
    otherwise
        disp('Only 1 argument supported. Use W=kendallsw(X). ')
    end
    
	options.rankfun = @(x)(tiedrank(x,0)); 
	
	[p n] = size(X); 
	ranks = zeros(p,n); 

	
	for ii=1:n	
		ranks(:,ii) 	= options.rankfun(X(:,ii));
		correction(ii) 	= tiedcorrection(ranks(:,ii));
	end
	
	sumR = sum(ranks,2);	
		
	if(sum(correction)==0)
		% no ties		
		varR  = sum((sumR-mean(sumR)).^2);
		W = 12*varR/((n^2)*p*(p^2-1));
	else
		% w. ties
		correction;
		W = (12*sum(sumR.^2) - 3*n^2*p*(p+1)^2)/...
							((n^2)*p*(p^2-1) - n*sum(correction)); 
	end
end

function correction =  tiedcorrection(ranks)
	
	u_ranks = unique(ranks); 
	n_unique = length(u_ranks);
	correction = 0;
	
	if(n_unique<length(ranks)) 
		for grp_no=1:n_unique
			tmp_idx = find(ranks==u_ranks(grp_no)); 
			n_ties = length(tmp_idx); 
			correction = correction + (n_ties^3-n_ties);
		end
	end
end


function smoothedrank(x)	
	
	
end