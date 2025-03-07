% This function made the installation of QuaLiKiz for METIS.
% A configuration must be created for each type of linux system.
% Files are directely downloaded from Git sever.
%
% syntax:
%
% path_to_NN_data = install_QLKNN10D(server,flag)
%
% with 
%      server is the server type name (default = 'intra' for IRFM server).
%      flag = compilation flag use by Tubs 
%
function path_to_QLKZ_data = install_QLKZ(server,flag)

% for Jenkins use
try
  addpath(pwd);
catch
  disp('path is already setted');
end

% somme information
disp('Start of installation of QuaLiKiz for METIS ');

% check input
if nargin == 0
    server = '';
    flag   = '';
elseif nargin < 2
    flag = '';
end

if isempty(server)
    [s,t] = unix('hostname -a');
    if s ~= 0
       error(t); 
    end
    if length(t) < 5
        [s,t] = unix('hostname -f');
        if s ~= 0
            error(t);
        end
    end
    if ~isempty(strfind(t,'.intra.'))
        server = 'irfm_intra';
        disp('Server detected = irfm_intra');
        disp('test of compilation has been done with: intel/2018 mpi/2018');
    end
end

if isempty(flag)
    switch server
        case 'irfm_intra'
            flag = 'TOOLCHAIN=intel FC=mpiifort LINK=mpiifort VERBOSE=1 QLK_HAVE_NAG=0';
    end
end

if ~isempty(flag)
    fprintf('Flag use to build the project are: %s\n',flag);
end

% memorise current directory
disp('Creating directory for project')
cdmem = pwd;
% make clean directory for the installation
cd(fileparts(which('install_QLKZ')));
fprintf('changing directory to %s\n',pwd);
if isdir('QLKZ')
  s = movefile('QLKZ',sprintf('%s_%s','QLKZ',datestr(now,30)));
  if s <= 0
    cd(cdmem);
    error('Not able to move previous directory');
   end
end

s = mkdir('QLKZ');
if s <= 0
    cd(cdmem);
    error('Not able to create QLKZ directory');
end

cd('QLKZ');

% download file from Git servers
disp('donwloading files from Git server');
s = unix(sprintf('%s','git clone https://gitlab.com/qualikiz-group/QuaLiKiz.git'));
if s ~= 0
    cd(cdmem);
    error('Git error during cloning of https://gitlab.com/qualikiz-group/QuaLiKiz.git');
end

cd('QuaLiKiz');
fprintf('changing directory to %s\n',pwd);
s = unix(sprintf('%s','git submodule update --init'));
if s ~= 0
    cd(cdmem);
    error('Git error during cloning submodule');
end
% checkout branch master
disp('Switch to branch main :')
!git checkout main
!git branch

%
disp('Cleaning the project');
[s,t] = unix(sprintf('make %s VERBOSE=1 clean',flag));
if s ~= 0
    cd(cdmem);
    if ~isempty(strfind(t,'make version'))
        error(sprintf('Wrong version of tool make in use; check which module you are using.\n%s\n',t));
    else
        error(t); 
    end
end
!rm -rf build
!rm -f bin/qlknn_mex.mexa64

disp('Start project compilation')
s = unix(sprintf('module purge \n module load intel/2018 mpi/2018 \n make %s VERBOSE=1',flag));
if s ~= 0
    cd(cdmem);
    error('Error during project compilation');
else
   disp('Installation of QuaLiKiz for METIS finished successfully')
end






