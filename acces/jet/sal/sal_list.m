% function providing data tree information from SAL server.
% 
% syntax:
%
%    data = sal_list(path,[authorisation,revision])
%
% input:
%    path                 = path to the data  (see SAL documentation)
%    authorisation        = authorisation data structure as return by sal_auth 
%                           (optional, useful only if multi servers connexion)
%    revision             = data version (default 0)
%
% output:  tree structure  (see SAL documentation for  details)
%
% reference: https://simple-access-layer.github.io/documentation/server/rest_api.html
% contact:   jean-francois.artaud@cea.fr
%
% test:
% data = sal_list('data/pulse/87737/ppf/signal/jetppf/magn');
%
function data = sal_list(data_path,authorisation_machine,revision)

if nargin < 1
    error('syntax: data = sal_list(data_path,[authorisation_machine,revision])');
elseif (nargin < 2) || isempty(authorisation_machine)
    authorisation_machine = 'JET';
end
if ~isstruct(authorisation_machine)
    authorisation_machine = getappdata(0,sprintf('SAL_AUTHORISATION_%s',authorisation_machine));
end

% test token validity
if authorisation_machine.validity < now
    try
        authorisation_machine = sal_auth;
    catch
        error('Token validity has expired: You must renew you authentication');
    end
end

% see details on https://simple-access-layer.github.io/documentation/server.html
if (nargin < 3) || isempty(revision)
    revision = 0; % default provided by the server
end
    

% parse data path for escape caracters
% Converts a node name from native JPF format to SAL format.
%SAL nodes have a limited character set so escape special characters.
data_path  = strrep(data_path,'>', '_out_');
data_path  = strrep(data_path,'<', '_in_');
data_path  = strrep(data_path,':', '_sq_');
data_path  = strrep(data_path,'$', '_do_');
data_path  = strrep(data_path,'&', '_et_');
data_path  = strrep(data_path,' ', '_sp_');
if data_path(1) ~= '/'
    data_path = strcat('/',data_path);
end

header{1,1} = 'Host';
header{1,2} = authorisation_machine.server;
header{2,1} = 'Authorization';
header{2,2} = sprintf('Bearer %s',authorisation_machine.token);
opt = weboptions('RequestMethod','get','CertificateFilename','','ContentType','text','HeaderFields',header,'Timeout',60);
url = sprintf('https://%s%s?revision=%d&auth=%s', ...
      authorisation_machine.server,data_path,revision,authorisation_machine.token);

try
    t = webread(url,opt);
catch    
    if ~isempty(strfind(lasterr,'NOT FOUND'))
        t = [];
    else
        % built wegt querry
        %
        disp('using backup solution with wget to retreive the data')
        [s,t] = unix(sprintf('wget -q -O - --header="Host: %s" --header="Authorization: Bearer %s" "https://%s%s?revision=%d&auth=%s"', ...
            authorisation_machine.server,authorisation_machine.token,authorisation_machine.server,data_path,revision,authorisation_machine.token));
        
        % error handling
        if s ~= 0
            error(t);
        end
    end
end

% decode json structure
if ~isempty(t)
    data = jsondecode(t);
else
    data = t;
    return;
end

% convert fields endcode in base64 to matlab native type
data = convert_from_base64(data);





         
         





