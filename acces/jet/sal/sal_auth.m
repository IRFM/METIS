% function negociating connexion autorisation to SAL server.
% return a session token valid for 4 hours.
% syntax:
%
%    authorisation = sal_auth(user,passwd,[machine])
%
% input:
%    user    = user name on SAL server
%    passwd  = password or SecurID token for JET server
%    machine = machine name or server adresse (default JET: sal.jet.uk).
%
% output: autorisation data structure with fields
%       user      = username
%       token     = token  of session
%       validity  = limit of validity in time (to be compared with date return by now)
%       server    = address of the SAL server
%
% reference: https://simple-access-layer.github.io/documentation/server/rest_api.html
% contact:   jean-francois.artaud@cea.fr
%
function authorisation = sal_auth(user,passwd,machine)

% default 
authorisation = [];


% do not reask for authorisation if there is still one valide
if nargin < 3
     machine = 'JET';
end   
switch machine
    case 'JET'
        server = 'sal.jet.uk';
    otherwise
        server = machine;
end

% test if is a existing authorisation for the machine
if isappdata(0,sprintf('SAL_AUTHORISATION_%s',machine))
    authorisation = getappdata(0,sprintf('SAL_AUTHORISATION_%s',machine));
    if ~isempty(authorisation) && isfield(authorisation,'validity') && (authorisation.validity >= now)
           try
                data = sal_get('data/pulse/87737/ppf/signal/jetppf/magn/ipla');
                if ~isempty(data)
                    return
                else
                    authorisation = [];
                end    
           catch
               authorisation = [];
           end
    else
        authorisation = [];
    end
end

if nargin == 0
    % GUI mode
	prompt={'User name','SecurID key'};
	def={'    prebut    ','1234567890'};
	dlgTitle='Identification for JET data access through SAL';

    lineNo=1;
    answer=zinputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(answer)
	    return
    end
    user  = strtrim(answer{1});
    securidkey = strtrim(answer{2});
    % try SAL authetification if securidkey is defined
    if (length(securidkey) >= 10) 
        try
            authorisation = sal_auth(user,securidkey);
        catch
            disp('SAL identification failed !');
        end
        if ~isempty(authorisation)
            disp('Connected to JET SAL server');
        end
    else
            disp('Wrong passwd format !');
    end
    return
elseif nargin < 2
    error('syntax: uthorisation = sal_auth(user,passwd,[machine]);');
% elseif nargin == 2
%     machine = 'JET';
end


% wget command
if isnumeric(passwd)
    passwd= sprintf('%d',passwd);
end

% import matlab.net.http.*
% creds = Credentials('Username',user,'Password',passwd);
% options = HTTPOptions('Credentials', creds,'CertificateFilename','');
% url = sprintf('https://%s/auth',server);
% [response, request] = RequestMessage().send(url,options);
% authorizationField = request.getFields('Authorization');
% authInfo = authorizationField.convert;
% disp(string(authInfo));
%auth_sal = base64encode(strcat(user,':',passwd),'');
% for debug
%!wget -v -d  --server-response --save-headers --http-user=username --http-password=pinc123456 https://sal.jet.uk/auth

%opt = weboptions('RequestMethod','get','Username',user,'Password',passwd,'CertificateFilename','','ContentType','text','debug',true);
opt = weboptions('RequestMethod','get','Username',user,'Password',passwd,'CertificateFilename','','ContentType','text','Timeout',60);
url = sprintf('https://%s/auth',server);

%try
    t = webread(url,opt);
    if isempty(t)
        error('t is empty');
    end
% catch
%     disp(lasterr)
%     disp('using backup solution with wget to connect to the server')
%     [s,t] = unix(sprintf('wget -q -O - --http-user=%s --http-password=%s https://%s/auth ',user,passwd,server));
%     if (s~= 0) && isempty(strfind(t,'{"authorisation": {'))
%         error(t);
%     end
% end
r=jsondecode(t);
if isfield(r,'authorisation') && isfield(r.authorisation,'token') && ~isempty(r.authorisation.token)
    authorisation = r.authorisation;
    authorisation.validity = now + 3.90 / 24;
    authorisation.server   = server;
    setappdata(0,sprintf('SAL_AUTHORISATION_%s',machine),authorisation)
else
    error('No vaild token !');
end
