% Install FACIT from sources files in metis/coef/FACIT_TOOLS/runs
% Currently only for IRFM zone partenaires or intranet
function path_to_QLKZ_data = install_FACIT(server)

% somme information
disp('Start of installation of FACIT for METIS ');

% check input
if nargin == 0
    server = '';
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
        disp('test of compilation has been done with:  inteloneapi/2025.1.0 lapack/3.9.1-inteloneapi-2025.1.0');
    end
end

% memorise current directory
disp('Creating directory for project')
cdmem = pwd;
% make clean directory for the installation
cd(fileparts(which('install_FACIT')));
fprintf('changing directory to %s\n',pwd);

cd('FACIT_TOOLS/runs');
fprintf('changing directory to %s\n',pwd);

disp('Start project compilation')
s = unix(sprintf('module purge \n module load inteloneapi/2025.1.0 lapack/3.9.1-inteloneapi-2025.1.0 \nifx -c constants.f90 -o constants.o \nifx -c facit.f90 -o facit.o -llapack \nifx -c FACIT_interface.f90 -o FACIT_interface.o \nifx -O2 -r8 -mcmodel=large facit.o FACIT_interface.o constants.o -o FACIT_interface -llapack'));

if s ~= 0
    cd(cdmem);
    error('Error during project compilation');
else
   disp('Installation of FACIT for METIS finished successfully')
end






