addpath('external')
if(exist('external/mpm','dir'))
	addpath('external/mpm')
end
if(exist('mpm'))
        disp('Invoking Matlab Package Manager')
else
        disp('Installing Matlab Package Manager');
        system('git clone https://github.com/mobeets/mpm.git external/mpm');		
end

disp('Installing package dependencies to external/.')
PYTHON_EXE = '/usr/bin/python';
MPM_INSTALL_DIR = fullfile(pwd,'external');
MPM_PATH = fullfile(pwd,'external','mpm');
try

	mpm_args =  ['  -r  ' which('requirements.txt') '  --installdir ' MPM_INSTALL_DIR];
	cmd = [PYTHON_EXE ' ' fullfile(MPM_PATH, 'mpm.py') ' ' mpm_args];
	[~, output] = system(cmd);   
    mpm_paths(MPM_INSTALL_DIR);
	disp(output);
catch me
        disp(me)
        warning('Could not install packages in requirements.txt. Please install them manually to external/ or add pre-existing packages to your path');
end
