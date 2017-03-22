function test_suite = test_standardize
	buildFunctionHandleTestSuite(localfunctions);

end


function test_unit_diagonal_null
	
	n = 1000; 
	p = 10;
	
	mu = zeros(1,p); 
	sig = ones(1,p);
	unit_sig = ones(1,p);
	errtol = 1e-1;
	tol = 1e-5;
	
	rng('default');
	X = repmat(mu,[n 1]) + randn(n,p)*diag(sig);
	
	[Xnew mu_hat sig_hat] = standardize.standardize_cols(X); 
	
	assert(any(abs(mean(Xnew))>tol)==0,'Not centered');
	assert(any(abs(std(Xnew,1)-unit_sig)>tol)==0,'Not scaled');
	
	assert(any(abs(mu_hat-mu)>errtol)==0,'Mean wrong');
	assert(any(abs(sig_hat-sig)>errtol)==0,'Std wrong');
	
	[Xnew mu_hat sig_hat] = standardize.standardize_cols(X,'center'); 
	
	assert(any(abs(mean(Xnew))>tol)==0,'Not centered');
	
	assert(any(abs(mu_hat-mu)>errtol)==0,'Mean wrong');
	assert(any(abs(sig_hat-sig)>errtol)==0,'Std wrong');
	
	[Xnew mu_hat sig_hat] = standardize.standardize_cols(X,'scale'); 
	
	assert(any(abs(std(Xnew,1)-unit_sig)>tol)==0,'Not scaled');
	
	assert(any(abs(mu_hat-mu)>errtol)==0,'Mean wrong');
	assert(any(abs(sig_hat-sig)>errtol)==0,'Std wrong');
end

function test_unit_diagonal_mean
	
	n = 400; 
	p = 10;
	
	mu = zeros(1,p);
	mu(2) = 10;	 
	sig = ones(1,p);
	unit_sig = ones(1,p);
	errtol = 1e-1;
	tol = 1e-5;
	
	rng('default');
	X = repmat(mu,[n 1]) + randn(n,p)*diag(sig);

	[Xnew mu_hat sig_hat] = standardize.standardize_cols(X); 
	
	assert(any(abs(mean(Xnew))>tol)==0,'Not centered');
	assert(any(abs(std(Xnew,1)-unit_sig)>tol)==0,'Not scaled');
	
	assert(any(abs(mu_hat(2)-mu(2))>errtol)==0,'Mean wrong');
	assert(any(abs(sig_hat-sig)>errtol)==0,'Std wrong');

	[Xnew mu_hat sig_hat] = standardize.standardize_cols(X,'center'); 
	
	assert(any(abs(mean(Xnew))>tol)==0,'Not centered');
	
	assert(any(abs(mu_hat(2)-mu(2))>errtol)==0,'Mean wrong');
	assert(any(abs(sig_hat-sig)>errtol)==0,'Std wrong');
	
	[Xnew mu_hat sig_hat] = standardize.standardize_cols(X,'scale'); 
	
	assert(any(abs(std(Xnew,1)-unit_sig)>errtol/2)==0,'Not scaled');
	
	assert(any(abs(mu_hat(2)-mu(2))>errtol)==0,'Mean wrong');
	assert(any(abs(sig_hat-sig)>errtol)==0,'Std wrong');
	
	[Xnew results] = standardize.successive_normalize(X); 
	
	assert(any(abs(mean(Xnew))>tol)==0,'Not centered');
	assert(any(abs(std(Xnew,1)-unit_sig)>errtol/2)==0,'Not scaled');
	
	results

	
end


function test_diagonal_null 


	
end


function test_diagonal_mean


	
end