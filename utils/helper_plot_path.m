function [h1,h2, path] = helper_plot_path(path)
% Plot the regularization or stability path for p variables. 	
% Inputs 
% path.lambdaList
% path.G - graphs/variables along full path
% path.Lambda - chosen lambda
% path.model - chosen model
% path.score - score from the model selection procedure
% path.colormap ? 
% path.topk ?
% 

	G = path.G; 
	if(isempty(path.model))
		model = G(:,:,end-2);
	else
		model = path.model;
	end
	lambdaList = path.lambdaList;
	lambda = path.lambda;

	[p p l1] = size(path.G); 
	G = reshape(G,[p*p l1]);
	l = length(path.lambdaList);
	if(l~=l1)
		warning('Graph and LambdaList not comparable')
	end


	[maxval maxind] = sort(reshape(abs(triu(model,1)),[p*p 1]),'descend');
	topk = 50; 
	
	cmap = flipud(winter(100)); 
	roundtopk = round(topk/100)*100;
	new_cmap = repmat(cmap, [1 1 round(roundtopk/100)]);
	new_cmap = permute(new_cmap,[3 1 2]);
	new_cmap = reshape(new_cmap,[100*round(roundtopk/100) 3]);
	cmap = new_cmap; clear new_cmap;
	
	
	% Find saturation index (where coefficients/sparsity saturates along path)
	
	
	h1=figure(1); hold on;
	for ii=1:min(length(maxind),topk)
		plot(lambdaList(1:end),(abs(G(maxind(ii),1:end))),'LineWidth',2,'Color',cmap(ii,:));	
	end
	colorbar;
	hchild = get(h1,'Children')
	set(hchild(1),'TickLabels',{num2str(maxval(1)), num2str(maxval(topk))},'Ticks',[0 1]);
	
	hold off; 

	d_bl = path.score;
	
	if(length(d_bl)==l1)
		h2 = figure(2); hold on; 
		plot(lambdaList,d_bl,'LineWidth',3,'Color',[.8 .1 .1]);
		title('Model Selection Scores Vs. Regularization Parameter'); 
		hold off; 
	else
		h2 = [];
	end