function IMAS_VERSION = get_imas_version

% get IMAS dictionnary major version
if isappdata(0,'IMAS_VERSION')
    IMAS_VERSION  = getappdata(0,'IMAS_VERSION');
else
    IMAS_VERSION  = [];
end
if isempty(IMAS_VERSION)
    IMAS_VERSION_ENV = getenv('IMAS_VERSION');
    if isempty(IMAS_VERSION_ENV)
        IMAS_VERSION = 3;
        setappdata(0,'IMAS_VERSION',IMAS_VERSION);
    else
        IMAS_VERSION_ENV = strtrim(IMAS_VERSION_ENV);
        IMAS_VERSION_ENV(abs(IMAS_VERSION_ENV) < abs(' ')) = [];
        mask = (abs(IMAS_VERSION_ENV) >= abs('0')) & (abs(IMAS_VERSION_ENV) <= abs('9'));
        IMAS_VERSION_ENV = IMAS_VERSION_ENV(mask);
        IMAS_VERSION = str2double(IMAS_VERSION_ENV(1));
        if isempty(IMAS_VERSION)
            IMAS_VERSION = 3;
        end
        setappdata(0,'IMAS_VERSION',IMAS_VERSION);
    end
end
if IMAS_VERSION > 3
        setappdata(0,'UAL_DATAVERSION',sprintf('%d',IMAS_VERSION));
end
