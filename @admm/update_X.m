function X = update_X(Sigma,X,Z,U,options)
	
	rho = options.rho;	
	
    [Q,L] = mexeig(rho*(Z-U) - Sigma);
	if(size(L,1)==size(L,2))
   		es = diag(L);
	else
		es = L; 
	end
    xi = (es + sqrt(es.^2 + 4*rho))./(2*rho);
    X = Q*diag(xi)*Q';   
	
end