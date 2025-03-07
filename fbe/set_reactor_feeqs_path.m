function set_reactor_feeqs_path(directoryname)

% get env variable (isenv does not work in older Matlab version)
try
    feeqs_path_env = getenv('FEEQS_PATH');
catch
    feeqs_path_env = '';
end

% test if FEEQS is in the path
if ~exist('reactor_inverse4metis_fast') || ~exist('start_up_Unix') || ...
        ~exist('extract_metisdata4feeqs')
    if nargin > 0
        if ~exist(directoryname,'dir')
            error('Invalid path to FEEQS.M code');
        end
    elseif ~isempty(feeqs_path_env)
        directoryname = feeqs_path_env;
    elseif ~isempty(which('isfolder')) && isfolder('/Applications/software/FEEQS.M')
        directoryname = '/Applications/software/FEEQS.M';
    elseif isempty(which('isfolder'))
        error('You ar currently used a too old Matlab version for runing FEEQS.M: please switch to newer than 2017b MAtlab version');
    else
        directoryname = uigetdir('*', 'select path to FEEQS root directory');
    end
    if isnumeric(directoryname)
        error('Path to FEEQS code must be set');
    else
        % added root directory
        addpath(directoryname)
        % try to compile if needed
        cdmem = pwd;
        try
            cd(directoryname);
            evalin('base','start_up_Unix');
            cd(cdmem);
            
        catch
            cd(cdmem);
            error('Error during FEEQS initialisation')
        end
        % added others directories
        addpath(fullfile(directoryname,'Projects','JT60_SA'));
        addpath(fullfile(directoryname,'Projects','JT60_SA','Lib'));
        addpath(fullfile(directoryname,'Projects','JT60_SA','Data'));
        % add path to fast code
        addpath(fullfile(directoryname,'Projects','FastCoilCurrentIdentification'));
        addpath(fullfile(directoryname,'Projects','FastCoilCurrentIdentification','Lib'));
        % add path to reactor project
        addpath(fullfile(directoryname,'Projects','Reactors'));
    end
end