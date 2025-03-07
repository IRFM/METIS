% This function made the installation of QLK-NN 10D for METIS.
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
function path_to_NN_data = install_QLKNN10D(server,flag)

% for Jenkins use
try
  addpath(pwd);
catch
  disp('path is already setted');
end

% somme information
disp('Start of installation of QLK-NN 10D for METIS (target = Matlab mexfile)');
disp('QLK10D for Matlab has been tested with Matlab 2017B');

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
        disp('test of compilation has been done with: gcc/8.2.0, matlab/2017B & make/4.2');
    end
end
if isempty(flag)
    switch server
        case 'irfm_intra'
            flag = 'TOOLCHAIN=gcc TUBSCFG_MPI=0 TUBSCFG_MKL=0 TUBSCFG_OPT="-mno-avx -mno-bmi2"';
    end
end
if ~isempty(flag)
    fprintf('Flag use to build the project are: %s\n',flag);
end

% memorise current directory
disp('Creating directory for project')
cdmem = pwd;
% make clean directory for the installation
cd(fileparts(which('install_QLKNN10D')));
fprintf('changing directory to %s\n',pwd);
if isdir('QLKNN10D')
  s = movefile('QLKNN10D',sprintf('%s_%s','QLKNN10D',datestr(now,30)));
  if s <= 0
    cd(cdmem);
    error('Not able to move previous directory');
   end
end
s = mkdir('QLKNN10D');
if s <= 0
    cd(cdmem);
    error('Not able to create QLKNN10D directory');
end
cd('QLKNN10D');

% change LD_LIBRARY_PATH to prenvent Curl problem
ldp = getenv('LD_LIBRARY_PATH');
% parse it to remove matlab
ldp = strsplit(ldp,':');
new_ldp = '';
for k = 1:length(ldp)
    if isempty(strfind(lower(ldp{k}),'matlab'))
        new_ldp = sprintf('%s%s:',new_ldp,ldp{k});
    end
end
new_ldp = new_ldp(1:end-1);
if ~isempty(strfind(getenv('SHELL'),'csh'))
    new_ldp = sprintf('setenv LD_LIBRARY_PATH %s\n',new_ldp);
else
    new_ldp = sprintf('export LD_LIBRARY_PATH=%s\n',new_ldp);    
end


% download file from Git servers
disp('donwloading files from Git server');
s = unix(sprintf('%s %s',new_ldp,'git clone https://gitlab.com/qualikiz-group/QLKNN-fortran.git'));
if s ~= 0
    cd(cdmem);
    error('Git error during cloning of https://gitlab.com/qualikiz-group/QLKNN-fortran.git');
end
s = unix(sprintf('%s %s',new_ldp,'git clone https://gitlab.com/qualikiz-group/qlknn-hyper-namelists.git'));
if s ~= 0
    cd(cdmem);
    error('Git error during cloning of https://gitlab.com/qualikiz-group/qlknn-hyper-namelists.git');
end
s = unix(sprintf('%s %s',new_ldp,'git clone https://gitlab.com/qualikiz-group/qlknn-hyper.git'));
if s ~= 0    
    cd(cdmem);
    error('Git error during cloning of https://gitlab.com/qualikiz-group/qlknn-hyper.git');
end
cd('QLKNN-fortran');
fprintf('changing directory to %s\n',pwd);
s = unix(sprintf('%s %s',new_ldp,'git submodule update --init'));
if s ~= 0
    cd(cdmem);
    error('Git error during cloning submodule');
end
% checkout branch master
disp('Switch to branch master :')
!git checkout master
!git branch

%
disp('Cleaning the project');
[s,t] = unix(sprintf('%s make %s VERBOSE=1 clean',new_ldp,flag));
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
% make TOOLCHAIN=gcc  TUBSCFG_MPI=0 TUBSCFG_MKL=0 VERBOSE=1 MEX=/Applications/Matlab_2013B/glnxa64/mex /ZONE_TRAVAIL/JA132999/METIS4ITM/metis4itm/trunk/coef/QLKNN10D/QLKNN-fortran/bin/qlknn_mex.mexa64
s = unix(sprintf('%s make %s VERBOSE=1 %s/bin/glnxa64/mex %s/bin/qlknn_mex.mexa64',new_ldp,flag,matlabroot,pwd));
if s ~= 0
    cd(cdmem);
    error('Error during project compilation');
else
    % copy file in coef of METIS
    sc = copyfile(sprintf('%s/bin/qlknn_mex.mexa64',pwd),fileparts(which('install_QLKNN10D')));    
    if sc <= 0
       error('unable to copy mexfile in METIS coef directory'); 
    end
end

% test 
cd(fileparts(which('install_QLKNN10D')));
load data4testingQLKNN10D.mat
path_to_NN_data = fullfile(fileparts(which('install_QLKNN10D')),'QLKNN10D','qlknn-hyper-namelists');
[Nout, dNout_dNin] = qlknn_mex(path_to_NN_data, data.Nin, 1, data.opts);
figure;
x = linspace(0,1,size(data.Nin,1));
subplot(2,2,1)
plot(x,data.Nout,'.r',x,Nout,'b');
subplot(2,2,2)
for k = 1:size(Nin,2)
    plot(x,data.dNout_dNin(:,:,k),'.r',x,dNout_dNin(:,:,k),'ob');
    hold on
end
subplot(2,2,3)
plot(x,data.Nout - Nout);
subplot(2,2,4)
for k = 1:size(Nin,2)
    plot(x,data.dNout_dNin(:,:,k) - dNout_dNin(:,:,k));
    hold on
end
cd(cdmem);
disp('Installation of QLK-NN 10D for METIS finished successfully')






