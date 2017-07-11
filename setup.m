%%%%%%%%%%%%%%%%%%%%
USE_LIGHTSPEED = true;
USE_TESTS=true;
USE_MEX=true;
%%%%%%%%%%%%%%%%%%%%
MATLABDIR =[getenv('HOME') filesep 'MATLAB'];
MATLIBPATH=[getenv('HOME') filesep 'MATLAB' filesep 'matlab-library'];
%%%%%%%%%%%%%%%%%%%%

% External packages specific to this class
addpath('external');
addpath(fullfile('external','kendalltau')); 
if(USE_LIGHTSPEED)
	addpath(fullfile('external','lightspeed')); 
end
if(USE_MEX)
	addpath(fullfile('external','mex')); 
end
addpath(genpath(fullfile('external',...
                'matlab-bgl')));
addpath(genpath(fullfile(MATLIBPATH,...
				'lib','KPMtools')));
addpath(genpath(fullfile(MATLIBPATH,...
				'lib','MATLAB-ParseArgs')));

% Necessary for GGM estimation
addpath('solvers');
addpath(fullfile('solvers','QUIC'));

% Tests
if(USE_TESTS)
	addpath(genpath(fullfile(MATLIBPATH,...
					'test','matlab-xunit')));
	addpath('test');
	if(exist(fullfile(MATLABDIR,'pggm-sims')))
		addpath(genpath(fullfile(MATLABDIR,'pggm-sims'))); 
	else
		warning('simulation package unavailable. All tests will not run');
	end
end