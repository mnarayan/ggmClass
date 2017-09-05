%% Generate problem data
randn('seed', 0);
rand('seed', 0);

n = 100;   % number of features 
N = 10*n;  % number of samples

% generate a sparse positive definite inverse covariance matrix
Sinv      = diag(abs(ones(n,1)));
idx       = randsample(n^2, 0.001*n^2);
Sinv(idx) = ones(numel(idx), 1);
Sinv = Sinv + Sinv';   % make symmetric
if min(eig(Sinv)) < 0  % make positive definite
    Sinv = Sinv + 1.1*abs(min(eig(Sinv)))*eye(n);
end
S = inv(Sinv);


USE_SIMULATION = false;
if(USE_SIMULATION)
    % generate Gaussian samples
    D = mvnrnd(zeros(1,n), S, N);
else
    % generate Gaussian samples
    bnustudy = load('~/MATLAB/datasets/CoRR/desikan/BNU1/task-rest_bold');
    D = standardize.successive_normalize(bnustudy.output{1}.Data(:,:,1)')';
    n = size(D,2);
end

%% Solve problem
Sighat = corr(D);
output = admm.admm_solver(Sighat);
X = output.Theta;
history = output.history; 

parcor = @(Theta)(2*eye(length(Theta))-diag(1./sqrt(diag(Theta)))*Theta*diag(1./sqrt(diag(Theta))));

%% Solve Adaptive problem: Caution. Currently noisier estimates
isAdaptive = true;
if(isAdaptive)
	disp('Checking weighted ADMM')
	new_options = output.options; 
	new_options.lambda = GGM.adaptiveWeights(output.Theta,1);
	new_output = admm.admm_solver(Sighat,new_options); 
	X = new_output.Theta; 
	history = new_output.history; 
end	

isPathwise = true;
if(isPathwise)
	new_options = output.options;
	maxLambda = max(max(triu(Sighat,1)));
    minLambda = .01;
	n_lambda = 30;
    new_options.lambda = linspace(minLambda,maxLambda,n_lambda);
	new_options.maxLambda = maxLambda; 
	new_options.pathwise = true;
	new_output = admm.admm_solver(Sighat,new_options); 
end

if(~isPathwise)
	X_admm = X;
	figure; imagesc(parcor(X_admm)); axis equal image; 
	%% Reporting
	K = length(history.objval);                                                                                                        
	h = figure;
	plot(1:K, history.objval, 'k', 'MarkerSize', 10, 'LineWidth', 2); 
	ylabel('f(x^k) + g(z^k)'); xlabel('iter (k)');

	g = figure;
	subplot(2,1,1);                                                                                                                    
	semilogy(1:K, max(1e-8, history.r_norm), 'k', ...
	    1:K, history.eps_pri, 'k--',  'LineWidth', 2); 
	ylabel('||r||_2'); 

	subplot(2,1,2);                                                                                                                    
	semilogy(1:K, max(1e-8, history.s_norm), 'k', ...
	    1:K, history.eps_dual, 'k--', 'LineWidth', 2);   
	ylabel('||s||_2'); xlabel('iter (k)'); 
else
	figure('Position',[75   150   900   650]); 
	n_cols = 5; 
	n_rows = n_lambda/n_cols;
	for ii=1:2:n_lambda
		X = new_output{ii}.('Theta');
		subplot(n_rows,n_cols,ii)
		imagesc(abs(parcor(X))>minLambda); axis equal image; 
		title(sprintf('Lambda = %.2f',new_output{ii}.('options').('lambda'))); 
	end
	
    sparsity = zeros(1,n_lambda);
    lambda_vals = sparsity;
    for ii=1:n_lambda
		X = new_output{ii}.('Theta');
        sparsity(ii) = sum(sum(abs(X)>=minLambda));
        lambda_vals(ii) = new_output{ii}.('options').('lambda');
    end
    figure;
    plot(lambda_vals,...
        sparsity,...
        '-',...
        'LineWidth',2 ...
        );
     
	figure('Position',[250   150   900   650]); 
	n_cols = 5; 
	n_rows = n_lambda/n_cols;
	for ii=1:n_lambda
		subplot(n_rows,n_cols,ii)
		history = new_output{ii}.('history');
		K = length(history.objval);                                                                                                        				
		plot(1:K, history.objval, 'k', 'MarkerSize', 10, 'LineWidth', 2); 
		ylabel('f(x^k) + g(z^k)'); xlabel('iter (k)');
		title(sprintf('Lambda = %.2f',new_output{ii}.('options').('lambda'))); 
	end
	
	figure('Position',[350   150   900   900]); 
	n_cols = 5; 
	n_rows = 2*n_lambda/n_cols;
	for ii=1:n_lambda
		subplot(n_rows,n_cols,ii)
		history = new_output{ii}.('history');
		K = length(history.objval);                                                                                                        				
		semilogy(1:K, max(1e-8, history.r_norm), 'k', ...
		    1:K, history.eps_pri, 'k--',  'LineWidth', 2); 
		ylabel('||r||_2'); 
		title(sprintf('Lambda = %.2f',new_output{ii}.('options').('lambda'))); 
		
		subplot(n_rows,n_cols,n_rows*n_cols/2+ii)		
		semilogy(1:K, max(1e-8, history.s_norm), 'k', ...
		    1:K, history.eps_dual, 'k--', 'LineWidth', 2);   
		ylabel('||s||_2'); xlabel('iter (k)'); 
		title(sprintf('Lambda = %.2f',new_output{ii}.('options').('lambda'))); 
	end
	
end
	