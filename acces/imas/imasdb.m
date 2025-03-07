% in new version of IMAS, imas_open is no more available
% replace imasdb script under shell
function [ual_tokamak,ual_user,ual_dataversion,ual_backend] = imasdb(ual_tokamak,ual_user,ual_dataversion,ual_backend)

% test input
if nargin < 1
    ual_tokamak = '';
end
if nargin < 2
    ual_user = '';
end
if nargin < 3
    ual_dataversion = '';
end
if nargin < 4
    ual_backend = '';
end
% test if environnement is defined
if ~isempty(ual_user)
    setappdata(0,'UAL_USER',ual_user);
elseif isappdata(0,'UAL_USER')
    ual_user = getappdata(0,'UAL_USER');
else
    ual_user = getenv('USER');
    setappdata(0,'UAL_USER',ual_user);
end
if ~isempty(ual_tokamak)
    setappdata(0,'UAL_TOKAMAK',ual_tokamak);
elseif isappdata(0,'UAL_TOKAMAK')
    ual_tokamak = getappdata(0,'UAL_TOKAMAK');
end
if ~isempty(ual_dataversion)
    setappdata(0,'UAL_DATAVERSION',ual_dataversion);
elseif isappdata(0,'UAL_DATAVERSION')
    ual_dataversion = getappdata(0,'UAL_DATAVERSION');
else
    ual_dataversion = getenv('IMAS_VERSION');
    setappdata(0,'UAL_DATAVERSION',ual_dataversion);
end
if ~isempty(ual_backend)
    setappdata(0,'UAL_BACKEND',ual_backend);
elseif isappdata(0,'UAL_BACKEND')
    ual_backend = getappdata(0,'UAL_BACKEND');
elseif ~isempty(getenv('IMAS_AL_BACKEND'))
    ual_backend = getenv('IMAS_AL_BACKEND');
end

% mode interactif or not
if (nargin > 0) && ~isempty(ual_tokamak)
    return
end

% with output
if (nargout > 0) && ~isempty(ual_tokamak) && ~isempty(ual_user) && ...
    ~isempty(ual_dataversion) && ~isempty(ual_backend)
    return
end

% information to open WEST public database
disp(' ');
disp('===================================================================')
disp('To access to WEST public database, fields must be:')
disp('    user    = imas_public');
disp('    tokamak = west');
disp(' ');
fprintf('Your user name is %s\n',getenv('USER'));
rep = dir('~/public/imasdb/');
if ~isempty(rep)
    list = {};
    for k=1:length(rep)
        if isempty(strmatch(rep(k).name,{'.','..'},'exact'))
            list{end+1} = rep(k).name;
        end
    end
    if ~isempty(list)
        disp(' ');
        disp('Your private available databases are:')
        for k = 1:length(list)
            fprintf('%s\t',list{k});
        end
        fprintf('\n');
    end
end
% search for available IMAS version
p = which('imas_open_env');
while ~isempty(strfind(p,'ual'))
    [p,v] = fileparts(p);
end
while ~isempty(strfind(p,'matlab'))
    [p,v] = fileparts(p);
end
while ~isempty(strfind(p,'mex'))
    [p,v] = fileparts(p);
end
%[p,v] = fileparts(p); % change in IMAS installation
rep = dir(p);
if ~isempty(rep)
    list = {};
    for k=1:length(rep)
        if isempty(strmatch(rep(k).name,{'.','..'},'exact'))
            list{end+1} = rep(k).name;
        end
    end
    if ~isempty(list)
        disp(' ');
        disp('Available IMAS versions are:')
        for k = 1:length(list)
            fprintf('%s\t',list{k});
        end
        fprintf('\n');
    end
end
disp(' ')
disp('List of available backend (string): ');
disp('     12 = MDS+')
disp('     13 = HDF5')
disp('     Inf or NaN = default systemp setting')
disp('===================================================================')


%-------------------------------------------------------------------------
% Check computer type to see if we can use a GUI
%-------------------------------------------------------------------------
useGUI   = 1; % Assume we can use a GUI

if isunix,     % Unix?
    useGUI = length(getenv('DISPLAY')) > 0;
    if useGUI
      [s,t] = unix('xhost');
      if s~= 0
        useGUI = false;
      end
    end
end % if


if useGUI
    prompt={'         Tokamak:               ','User:','Version:','Backend:'};
    def={ual_tokamak,ual_user,ual_dataversion,ual_backend};
    dlgTitle='Set IMAS environment:';
    lineNo=1;
    answer=zinputdlg(prompt,dlgTitle,lineNo,def);
    drawnow
    if isempty(answer)
        
        warning('IMAS environment setting cancelled');
        return
    end
    ual_tokamak        = strtrim(answer{1});
    if isempty(ual_tokamak)
        ual_tokamak = 'test';
        disp('Using default value "test" as machine name');
    end
    ual_user           = strtrim(answer{2});
    if isempty(ual_user)
        ual_user = getenv('USER');
        fprintf('Using default user %s\n',ual_user);
    end
    ual_dataversion    = strtrim(answer{3});
    if isempty(ual_dataversion)
        ual_dataversion = list{end};
        fprintf('Using default ual version %s',ual_dataversion);
    end
    ual_backend    = strtrim(answer{4});
    if isempty(ual_backend)
        ual_backend = NaN;
        fprintf('Using default ual backend');
    end
    setappdata(0,'UAL_TOKAMAK',ual_tokamak);
    setappdata(0,'UAL_USER',ual_user);
    setappdata(0,'UAL_DATAVERSION',ual_dataversion);
    setappdata(0,'UAL_BACKEND',ual_backend);
else
    ual_tokamak = input(sprintf('IMAS environment Tokamak (%s) ? ',ual_tokamak),'s');
    if ~isempty(ual_tokamak)
        setappdata(0,'UAL_TOKAMAK',strtrim(ual_tokamak));
    else
        disp('Using default value "test" as machine name');
        setappdata(0,'UAL_TOKAMAK','test');
    end
    ual_user = input(sprintf('IMAS environment User (%s) ? ',ual_user),'s');
    if ~isempty(ual_user)
        setappdata(0,'UAL_USER',strtrim(ual_user));
    else
        setappdata(0,'UAL_USER',strtrim(getenv('USER')));
        fprintf('Using default user %s\n',getenv('USER'));
    end
    ual_dataversion = input(sprintf('IMAS environment Version (%s) ? ',ual_dataversion),'s');
    if ~isempty(ual_dataversion)
        setappdata(0,'UAL_DATAVERSION',strtrim(ual_dataversion));
    else
        ual_dataversion = list{end};
        setappdata(0,'UAL_DATAVERSION',ual_dataversion);
        fprintf('Using default ual version %s',ual_dataversion);        
    end
    ual_backend = input(sprintf('IMAS backend (%s) ? ',ual_backend),'s');
    if ~isempty(ual_backend)
        setappdata(0,'UAL_BACKEND',strtrim(ual_backend));
    else
        disp('Using default system setting');
        setappdata(0,'UAL_BACKEND','NaN');
    end

end
