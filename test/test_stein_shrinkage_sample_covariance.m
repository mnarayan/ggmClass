%% Generate problem data
randn('seed', 0);
rand('seed', 0);

n = 100;   % number of features 
N = .8*n;  % number of samples

% generate a sparse positive definite inverse covariance matrix
Sinv      = eye(n);
idx       = randsample(n^2, 0.001*n^2);
Sinv(idx) = -ones(numel(idx), 1);
Sinv(find(eye(n))) = 1;
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
Sighat = cov(D);
opts = estimator.create_options();
output = estimator.stein_shrinkage_sample_covariance(Sighat,opts);