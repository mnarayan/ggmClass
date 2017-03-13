function [varargout] = debias(self,varargin)
%DEBIAS returns debiased sparse inverse covariance estimate using the  Jankova-VanDeGeer Method
% 	
% 	SYNTAX: 
% 	Call on a GGM object instance self
% 		self.debias()
% 	Call using an options structure 
% 		output = ggm.debias(options)
% 

	contents = fieldnames(self); 
	if(~sum(strcmp(contents,'Theta'))|~sum(strcmp(contents,'Sigma')))
		mfile_showhelp;
		error(sprintf('self.Theta and self.Sigma needed'));
	end	

	hessianProj = @(Sigma,Theta)(Theta'*Sigma*Theta);
	debiaser 	= @(Sigma,Theta)(Theta + Theta' - hessianProj(Sigma,Theta)); 
	
	if(isobject(self))	
		sprintf('This is a %s class object',mfilename('class'));		
		self.Theta = debiaser(self.Sigma,self.Theta);	
	else		
		self.Theta_debias = debiaser(self.Sigma,self.Theta);		
		if(nargout==1)
			varargout{1} = self;	
		end	
	end


end