% in new version of IMAS, imas_create is no more available
function idx = imas_create( name, shot, run, refshot, refrun )
% idx  = imas_create(name, shot, run, refshot, refrun )
% Create a database.
%
% idx  : database index
% name : name of the database (by convention ids).
% shot : shot number.
% run  : run number.
% refshot : currently not use
% refrun  : currently not use

import imasjava.*

if (nargin < 3 || nargin >5)
    error('Bad number of input arguments. must be 5');
end

if (nargin<5)
	refrun=0;
end

if (nargin<4)
	refshot=0;
end

if ~ischar(name);
    error('First input argument must be a string (BD name)')
end
if ~isnumeric(shot)
        error('Second input argument must be a numeric (shot num)')
end
if ~isnumeric(run)
        error('Third input argument must be a numeric (run number)');
end

% test if environment is defined
envmustbesetted = 0;
if ~isappdata(0,'UAL_USER')
  envmustbesetted = 1;
end
if ~isappdata(0,'UAL_TOKAMAK')
  envmustbesetted = 1;
end
if ~isappdata(0,'UAL_DATAVERSION')
  envmustbesetted = 1;
end
if envmustbesetted
  imasdb;
end

% get backend
id_backend = get_imas_backend_id;
if isempty(id_backend)
    %idx = imas.create(name, shot, run, refshot, refrun);
    idx = imas_create_env( name, shot, run, refshot, refrun, strtrim(getappdata(0,'UAL_USER')),strtrim(getappdata(0,'UAL_TOKAMAK')),strtrim(getappdata(0,'UAL_DATAVERSION')));
else
    % use now imas_create_env_backend
    % in backup switch to other backend
    try
        idx = imas_create_env_backend(shot, run, refshot, refrun, strtrim(getappdata(0,'UAL_USER')), ...
              strtrim(getappdata(0,'UAL_TOKAMAK')),strtrim(getappdata(0,'UAL_DATAVERSION')),id_backend);    
    catch
        idx = imas_create_env( name, shot, run, refshot, refrun, strtrim(getappdata(0,'UAL_USER')),strtrim(getappdata(0,'UAL_TOKAMAK')),strtrim(getappdata(0,'UAL_DATAVERSION')));
        fprintf('selected backend (%d)is nos available switching back to default backend: %d\n',id_backend,imas_get_backendID(idx));
    end
end
