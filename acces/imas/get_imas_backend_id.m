function id_backend = get_imas_backend_id(ual_backend)

if isempty(which('imas_open_env_backend'))
    % no switch to HDF5 as defaut
    id_backend = [];
    return
end

if (nargin == 0) || isempty(ual_backend);
    if isappdata(0,'UAL_BACKEND')
        ual_backend = getappdata(0,'UAL_BACKEND');
    elseif ~isempty(getenv('IMAS_AL_BACKEND'))
        ual_backend = getenv('IMAS_AL_BACKEND');
        setappdata(0,'UAL_BACKEND',ual_backend)
    else
        % switch to defaut define by ITER I/O if available
        if isempty(which('imas_open_env_backend'))
            ual_backend = '';
        else
            ual_backend = 'HDF5';           
        end
    end
end
if isnumeric(ual_backend)
    ual_backend = sprintf('%d',ual_backend);
end
switch ual_backend
    case {'12','MDS+','mds+'}
        id_backend = 12;
    case {'13','HDF5','hdf5'}
        id_backend = 13;        
    otherwise
        id_backend = [];
end