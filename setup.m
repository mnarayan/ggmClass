addpath('@GGM');
addpath('@admm');
addpath('solvers');
addpath('solvers/QUIC');

addpath('external');
addpath('external/kendalltau'); 
USE_LIGHTSPEED = false;
if(USE_LIGHTSPEED)
	addpath('external/lightspeed'); 
end

% Tests
USE_TESTS=false;
if(USE_TESTS)
	addpath('test');
	addpath(genpath('../pggm-sims')); 
end