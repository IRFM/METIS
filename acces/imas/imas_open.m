% in new version of IMAS, imas_open is no more available
function expIdx = imas_open(name,shot,run)

% test input
if nargin < 3
  error('syntax: expIdx = imas_open(name,shot,run)')
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
    % use now imas_open_env
    expIdx = imas_open_env(name,shot,run,strtrim(getappdata(0,'UAL_USER')),strtrim(getappdata(0,'UAL_TOKAMAK')),strtrim(getappdata(0,'UAL_DATAVERSION')));
else    
    % use now imas_open_env_backend
    % in backup switch to other backend
    try
        expIdx = imas_open_env_backend(shot,run,strtrim(getappdata(0,'UAL_USER')),strtrim(getappdata(0,'UAL_TOKAMAK')),strtrim(getappdata(0,'UAL_DATAVERSION')),id_backend);
    catch
        % go back to MDS+ if not available
        expIdx = imas_open_env(name,shot,run,strtrim(getappdata(0,'UAL_USER')),strtrim(getappdata(0,'UAL_TOKAMAK')),strtrim(getappdata(0,'UAL_DATAVERSION')));
        fprintf('selected backend (%d)is nos available switching back to default backend: %d\n',id_backend,imas_get_backendID(expIdx));
    end
end

