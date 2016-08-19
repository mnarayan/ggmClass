function [scores lrange ind3] = stars(self,obj,grphs)	
	% Applies StARS using path of resampled graphs given by grphs 
	% INPUTS
	% 	self 	- KGGM object
	% 	obj 	- structure of parameters from MONET
	% 	grphs 	- a cell array of size (1,nresamples) consisting of graphs p x p x nlambdas x obj.s
	%
	% NOTE: self is a KGGM object and obj is a struct with the following parameters
	%
	%%%% Old Code %%%%
	% obj.n 	- number of time-series observations 
	% obj.p 	- number of variables
	% obj.s 	- number of graphs or subjects for multi-subject setting. For group models set to 1. 
	% obj.dm 	- data matrix n x p x s
	% obj.nlambda - length of path 
	% obj.estMethod = 'QUIC'
	% obj.beta = stars parameter
	% 
	% Usage, QUIC replaced with call to Glasso_Helper
	
	% TO DO switch everything over to self, grphs only
	
    grphs2 = cell(1,obj.s); % ultimately setSelParams()
	grphs3 = cell(1,obj.s);
    psi_grphs = cell(1,obj.s);
    theta_grphs = cell(1,obj.s);
    Dhat = cell(1,obj.s);
    lambdas = cell(1,obj.s);
	nresamples = length(grphs); 

    disp(sprintf('Selecting Optimal Penalty ...\n'))
    for cc=1:obj.s
            theta_ustat = zeros(obj.p,obj.p,obj.nlambda); 
			%disp(sprintf('... \n'));
            %percount(max((cc-1),0)*(obj.nlambda)+2,(obj.s)*(obj.nlambda));      												                  
            parfor ii=1:length(grphs)
                    % if(mod(ii,11)==0)        
                    %         percount(max((cc-1),0)*length(grphs)+ii+1,(obj.s)*length(grphs));  
                    % end  
                    theta_ustat = theta_ustat+(abs(squeeze(grphs{ii}(:,:,:,cc)))>0);                        
            end
			
            for jj=1:obj.nlambda
                    theta_grphs{cc}(:,:,jj) = triu(squeeze(theta_ustat(:,:,jj)),1)/nresamples;
            end
            psi_grphs{cc} = 2*((theta_ustat/nresamples).*(1-(theta_ustat/nresamples)));
            
            Dhat{cc} = zeros(1,length(obj.nlambda));
            
			for jj=1:obj.nlambda
                    Dhat{cc}(jj) = sum(sum(triu(squeeze(psi_grphs{cc}(:,:,jj)),1)))/nchoosek(obj.p,2); 
                    % Total instability over all edges
            end
            
            %-- Choose final lambda parameter
            D_bl(1) = Dhat{cc}(1);
            D_bl(2:obj.nlambda) = max(Dhat{cc}(1:end-1),Dhat{cc}(2:end));
            
            [finalL ind] = find((D_bl <= obj.beta)); %% Does this change anything ?
            [finalL2 ind2] = max(D_bl(ind));
            
            % if(ind2>0)
            %         lambdas{cc} = obj.lrange(ind2); % NOTE lrange for each subject needed
            % else
            %         lambdas{cc} = obj.lrange(obj.nlambda); % NOTE lrange for each subject needed
            %         warning('Total Instability > Tolerance beta. Try increasing beta');
            % end
            
            %-- Set 'pathmode' to 0. Find the final graph estimates and chosen lambda value                         
            if(strcmp(obj.estMethod,'QUIC'))
				[theta, lrange] = self.Glasso_Helper(obj.dm,0,size(obj.dm,3),0,0); 
            else
               error('Estimation method in obj.estMethod is not supported. Please use obj.estMethod="QUIC" '); 
            end     
            
            [finalL3 ind3] = min(abs(lrange-lrange(ind2)));
            %disp(sprintf('Old Lambda: %f, New Lambda %f',old_l,obj.lrange(ind3))); 
		    disp(sprintf('Penalty Parameters Selected. \n'))			
            lambdas{cc} = lrange(ind3)                  
            grphs2{cc} = squeeze(theta(:,:,ind3));
			grphs3{cc} = theta;
               
    end