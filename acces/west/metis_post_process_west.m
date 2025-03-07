% compute diagnostic synthetic
function metis_post_process_west(shot,run,occurrence,user,tokamak,dataversion,backend)

%  pre proccess equilibrium for Tofu
disp('Pre proccessing of IDS equilibrium for Tofu')
if isempty(occurrence)
    occurrence = '0';
elseif isnumeric(occurrence)
    occurrence = sprintf('%d',occurrence);
end
%[s,t] = unix(sprintf('module purge\nmodule load tools_dc\nTEQUIshiftfast %d %d %s %s %s %s',shot,run,occurrence,user,tokamak));
[s,t] = unix(sprintf('TEQUIshiftfast %d %d %s %s %s %s',shot,run,occurrence,user,tokamak));
if s == 0
    disp('Pre proceccessing of IDS equilibrium for Tofu done');
else
    disp('error during pPre proceccessing of IDS equilibrium for Tofu');
    error(t);
end


if nargin < 3
    occurrence = '';
end
if nargin < 4
    user = getenv('USER');
elseif isempty(user)
    user = getenv('USER');
end
if nargin < 5
    tokamak = 'west';
elseif isempty(tokamak)
    tokamak = 'west';
end

% select diagnostic to be computed
list_diag = {'bolo','sxr','interf','brem'};
[s,v] = listdlg('PromptString',{'Select diagnostic synthetics','to be computed:'},...
    'SelectionMode','multiple',...
    'ListString',list_diag);
if (v == 0) || isempty(s)
    return
end


% preprocessing of equilibrium
% TEQUIshiftfast 53838 2 0 JM240179 west
% usage: equilibrium_imas_shift_time.py [-h] [--dry-run]
%                                      [shot] [run] [occ] [user] [machine]
% 1- test if derectory exist otherwise create it
%    ~/public/equilibrium_data_shift_time/$machine/0
if ~isdir(fullfile(getenv('HOME'),'public'))
    [se,m] = mkdir(ullfile(getenv('HOME'),'public'));
    if se~= 1
        error(m);
    end
end
if ~isdir(fullfile(getenv('HOME'),'public','equilibrium_data_shift_time'))
    [se,m] = mkdir(fullfile(getenv('HOME'),'public','equilibrium_data_shift_time'));
    if se~= 1
        error(m);
    end
end
if ~isdir(fullfile(getenv('HOME'),'public','equilibrium_data_shift_time',strtrim(tokamak)))
    [se,m] = mkdir(fullfile(getenv('HOME'),'public','equilibrium_data_shift_time',strtrim(tokamak)));
    if se~= 1
        error(m);
    end
end
if ~isdir(fullfile(getenv('HOME'),'public','equilibrium_data_shift_time',strtrim(tokamak),'0'))
    [se,m] = mkdir(fullfile(getenv('HOME'),'public','equilibrium_data_shift_time',strtrim(tokamak),'0'));
    if se~= 1
        error(m);
    end
end

% % 2- call convertisseur
% occurrence = strtrim(occurrence);
% if ~isempty(occurrence)
%     conv_cmd = sprintf('TEQUIshiftfast %d %d %s %s %s ',shot,run,occurrence,user,tokamak);
% else
%     conv_cmd = sprintf('TEQUIshiftfast %d %d 0 %s %s ',shot,run,user,tokamak);
% end
% cmd = sprintf('module purge\nmodule load tools_dc\nwhich TEQUIshiftfast\n%s',conv_cmd);
% fprintf('Conversion command:\n%s\n\n',cmd);
% se   = unix(cmd);
% if se ~= 0
% 	  warndlg('Error in equilibrium preproccessing','TEQUIshiftfast');
% 	  return
% end

%2- delete cache file if exist
%data_IMAS_equi_Shot53643_Run0001_Occ0_JA132999_west_temp.mat
occurrence = strtrim(occurrence);
if isempty(occurrence)
    file_name = sprintf('data_IMAS_equi_Shot%d_Run%4.4d_Occ0_%s_%s.mat', ...
        shot,run,user,tokamak);
else
    file_name = sprintf('data_IMAS_equi_Shot%d_Run%4.4d_Occ%s_%s_%s.mat', ...
        shot,run,occurrence,user,tokamak);
end
file_name = fullfile(getenv('HOME'),'public','equilibrium_data_shift_time',strtrim(tokamak),'0',file_name);
if exist(file_name)
    delete(file_name)
end

% create the command
% tfw_METIS_SynthDiags.py [-h] [-diag DIAG1,DIAG2, ...] [-shot SHOT] [-run RUN]
%                               [-occ OCC] [-usr USR] [-machine MACHINE]
log_file = sprintf('post_process_west_%d_%d.log',shot,run);
%exe_cmd = fullfile(fileparts(which('metis_post_process_west')),'py','tfw_METIS_SynthDiags.py');
exe_cmd    = '/Applications/tools/tools_dc/tofu_tools_dc/tfw_METIS_SynthDiags.py';

occurrence = strtrim(occurrence);
if length(s)  == 1
    if ~isempty(occurrence)
        tofu_cmd = sprintf('python %s -shot %d -run %d -occ %s -usr %s -database %s -diag %s',exe_cmd,shot,run,occurrence,user,tokamak,list_diag{s});
    else
        tofu_cmd = sprintf('python %s -shot %d -run %d -usr %s -database %s -diag %s',exe_cmd,shot,run,user,tokamak,list_diag{s});
    end
% elseif length(s)  == length(list_diag)
%     if ~isempty(occurrence)
%         tofu_cmd = sprintf('python %s -shot %d -run %d -occ %s -usr %s -database %s',exe_cmd,shot,run,occurrence,user,tokamak);
%     else
%         tofu_cmd = sprintf('python %s -shot %d -run %d -usr %s -database %s' ,exe_cmd,shot,run,user,tokamak);
%     end
%     
else
    if ~isempty(occurrence)
        tofu_cmd = sprintf('python %s -shot %d -run %d -occ %s -usr %s -database %s -diag %s',exe_cmd,shot,run,occurrence,user,tokamak,list_diag{s(1)});
    else
        tofu_cmd = sprintf('python %s -shot %d -run %d -usr %s -database %s -diag %s',exe_cmd,shot,run,user,tokamak,list_diag{s(1)});
    end
    for k=2:length(s)
        if ~isempty(occurrence)
            tofu_cmd = sprintf('%s ; python %s -shot %d -run %d -occ %s -usr %s -database %s -diag %s',tofu_cmd,exe_cmd,shot,run,occurrence,user,tokamak,list_diag{s(k)});
        else
            tofu_cmd = sprintf('%s ; python %s -shot %d -run %d -usr %s -database %s -diag %s',tofu_cmd,exe_cmd,shot,run,user,tokamak,list_diag{s(k)});
        end
    end
end

%cmd = sprintf('module purge\nmodule load tools_dc\n%s ',tofu_cmd);
cmd = tofu_cmd;
fprintf('Tofu command:\n%s\n\n',cmd);
%cmd = sprintf('module purge\nmodule load tools_dc\n%s | tee %s',tofu_cmd,log_file);
se   = unix(cmd);
if se ~= 0
    t = sprintf('Error in Tofu with command: %s',cmd);
    warndlg(t,'Error launching post treatement');
end
