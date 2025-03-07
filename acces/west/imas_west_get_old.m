function [ imas_out ] = imas_west_get(shot, ids_name, varargin)
% ========================
% IMAS WEST get IDS script
% ========================
% Short usage:
% -----------
% [imas_out] = imas_west_get(shot, ids_name)
% To see all ids_name, get data dictionary, type in Linux Terminal: dd_doc
%
% Long usage:
% ----------
% [imas_out] = imas_west_get(shot, ids_name, imas_run, imas_occurrence, imas_user, imas_machine)
%
% Inputs:
% ------
% With defaults: shot, ids_name, imas_run=0, imas_occurrence=0,
%                imas_user='imas_public', imas_machine='west'
%
% Output:
% ------
% imas_out: requested IDS Matlab structure


% Inputs
% ------
% Input variables check:
% only want 3 optional inputs at most

numvarargs = length(varargin);
if (numvarargs > 4)
    error('ERROR:TINTimas_read:TooManyInputs', ...
        'requires at most 6 inputs, for help type >>help imas_west_get');
end
% If shot number empty or 0
if (nargin < 2)
    ME = MException('wrapper_TINTlidcalc:NoInputs', ...
        'ERROR: script requires at least shot number and IDS name as inputs');
    throw(ME);
end

% set defaults for optional inputs
optargs = {0 0 'imas_public' 'west'};

% now put these defaults into the valuesToUse cell array,
% and overwrite the ones specified in varargin.
optargs(1:numvarargs) = varargin;

% Place optional args in memorable variable names
[imas_run, imas_occurrence, imas_user, imas_machine] = optargs{:};

% Print
disp(' ')
fprintf('shot            = %d \n', shot);
fprintf('ids_name        = %s \n', ids_name);
fprintf('imas_run        = %d \n', imas_run);
fprintf('imas_occurrence = %d \n', imas_occurrence);
fprintf('imas_user       = %s \n', imas_user);
fprintf('imas_machine    = %s \n', imas_machine);
fprintf('Reading data... \n\n');

% Open IMAS file
% --------------
% Check if file exists
ids_file = ['~', imas_user, '/public/imasdb/', imas_machine, '/3/0/ids_', ...
    int2str(shot), sprintf('%04d', imas_run), '.datafile'];
status = system(['test -e ', ids_file]);
% Open file if exists
if (~status)
    % Function  imas_open_env('ids', shot, imas_run, imas_user, imas_machine, imas_version_number);
    idx = imas_open_env('ids', shot, imas_run, imas_user, imas_machine, '3');
else
    disp(' ')
    ME = MException('imas_west_get:IDSNotExist', ...
        'ERROR IDS file does not exist for this shot, run, user or machine');
    throw(ME);
end

% For interfero IDS
% -----------------
% Check IMAS version and transform it to integer
if (strcmp(ids_name, 'interferometer'))
    imas_version = str2double(strrep(getenv('IMAS_VERSION'), '.', ''));
    if (imas_version >= 3150)
        ids_name = 'interferometer';
    else
        ids_name = 'interfero_polarimeter';
    end
end

% Get IDS data with:
% function ids_get(id_from_imas_open, ids_name) for different occurences
if (~strcmp(ids_name, 'equilibrium'))
    if (imas_occurrence == 0)
        imas_out = ids_get(idx, ids_name);
    else
        ids_name_occ = [ids_name, '/', num2str(imas_occurrence)];
        imas_out = ids_get(idx, ids_name_occ);
    end
    fprintf('Data provider: %s \n\n', ...
        imas_out.ids_properties.provider)
    fprintf('Time of creation of the IDS: %s \n\n', ...
        imas_out.ids_properties.creation_date)
end

% For equilibrium IDS
% -------------------
% run python script if equilibrium mat file exists does not exist
if (strcmp(ids_name, 'equilibrium'))
    current_file = mfilename('fullpath');
    fprintf('current_file = %s', current_file)
    [filepath, ~] = fileparts(current_file);
    fprintf('filepath = %s', filepath)
    equi_mat = sprintf([filepath ...
        '/../tmp/matlab_equilibria/data_vactheqx_Shot%d_Run%d_Occ%d_%s_%s.mat'], ...
        shot, imas_run, imas_occurrence, imas_user, imas_machine);
    disp(' ')
    fprintf('\n equi mat file: %s \n \n', equi_mat)
    if ~exist(equi_mat, 'file')
        % read data
        ret_unix = unix(sprintf(['python ' filepath ...
            '/equilibrium_python_to_matlab.py %d %d %d %s %s'], ...
            shot, imas_run, imas_occurrence, imas_user, imas_machine));
        if (ret_unix ~= 0)
            ME = MException('imas_west_get:NoFileFound', ...
                'unable to read EQUINOX data');
            throw(ME);
        end
    end
    if exist(equi_mat, 'file')
        load(equi_mat);
        imas_out.timeEquiIDS = timeEquiIDS;
        imas_out.wall = wall;
        imas_out.ip = ip;
        imas_out.q_95 = q_95;
        imas_out.q_axis = q_axis;
        imas_out.li_3 = li_3;
        imas_out.beta_pol = beta_pol;
        imas_out.w_mhd = w_mhd;
        imas_out.mag_ax_R = mag_ax_R;
        imas_out.mag_ax_Z = mag_ax_Z;
        imas_out.RNodes = RNodes;
        imas_out.ZNodes = ZNodes;
        imas_out.Psi_val = Psi_val;
        imas_out.boundPlasma = boundPlasma;
        imas_out.xPoint = xPoint;
        imas_out.RInterp = RInterp;
        imas_out.ZInterp = ZInterp;
        imas_out.Psi_interp = Psi_interp;
        imas_out.rho_tor = rho_tor;
        imas_out.q_safFac = q_safFac;
        imas_out.elong = elong;
        imas_out.triang_up = triang_up;
        imas_out.triang_low = triang_low;
        imas_out.j_tor = j_tor;
        imas_out.pressure = pressure;
        imas_out.f_df_dpsi = f_df_dpsi;
        imas_out.dpress_dpsi = dpress_dpsi;
        imas_out.psi_prof = psi_prof;
        imas_out.limPoint = limPoint;
    else
        ME = MException('imas_west_get:NoFileFound', ...
            'ERROR IMAS mat file not found');
        throw(ME);
    end
end

disp('End imas west get')
end
