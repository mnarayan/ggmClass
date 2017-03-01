classdef admm < handle

	methods(Static)
		
		output 			= admm_solver(Sigma,varargin); 
		
		options 		= admm_options(varargin);		
		
		options 		= initializeXZU(Sigma,varargin); 
   	  	
		output  		= objective(Sigma,Theta,Z,Lambda);
		
		[X,Z] 			= update_X(Sigma,X,Z,U,options); 	
		
		[Z_new,X_new] 	= update_Z(X,Z,U,options); 					
		
		U 				= update_U(X,Z,U); 
		
	end
			 		
end