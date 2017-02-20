function output = objective(Sigma, Theta, Z, Lambda)	
	Res = Z.*Lambda;
	Res = setdiagLS(Res,0); 
	if(exist('logdetLS'))
   		output = trace(Sigma*Theta) - logdetLS(Theta) + norm(Res(:), 1);		
	else
   		output = trace(Sigma*Theta) - log(det(Theta)) + norm(Res(:), 1);
 	end
end