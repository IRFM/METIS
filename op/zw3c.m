% ZW3C acces au web
%-------------------------------------
% fichier zw3c.m ->  zw3c
%
%
% fonction Matlab 5 :
% 
%  Cette fonction permet de rammener des donnees via un serveur http.
%  Elle soumet un formulaire a l'adresse donnees par l'url et recupere 
%  le resultat. Elle utilise la focntion w3c du consotium W3C.
%
% syntaxe  :
%  
%  [status,reponse,cmd] = zw3c(url,method,name1,data1,...nameN,dataN);
%  
% entrees :
%   
%    url    = URL du serveur (complete)
%    method = methode pour envoyer les donnees :
%                 method = 'get' ou vide  -> methode get
%                 method = 'post'         -> methode post
%                 
%    nameX  = nom du Xieme champ du formulaire
%    dataX  = valeur du Xieme champ du formulaire
%                 
%    
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.6, du 19/07/2001.
% 
% 
% liste des modifications :
%
%--------------------------------------------------------------
%
function [status,reponse,cmd] = zw3c(url,method,varargin)


if nargin  < 2
    method = '';
elseif isempty(method)
    method = '';
end
form ='';
formp ='';
if nargin >=4
    for k = 1:2:(nargin-3)
       if ischar(varargin{k+1})
         chaine=varargin{k+1};
	      % translation
         chainec=strrep(chaine,'=','\=');			       	 
         form = [form,sprintf(' "%s=%s"',varargin{k},varargin{k+1})];
       else
         form = [form,sprintf(' "%s=%g"',varargin{k},varargin{k+1})];
       end
    end
end

% localisation executable
exe='/usr/drfc/cgc/matlab5/zineb/tar/w3c-libwww-5.3.2/ComLine/src/w3c';
%exe='/usr/drfc/cgc/matlab5/zineb/tar/w3c-libwww-5.3.2/ComLine/src/.libs/lt-w3c'; 
temp ='';
if strcmp(lower(method),'post') & ~isempty(form)
    cmd = sprintf('%s -n -post %s -form %s',exe,url,form);
elseif ~isempty(form)
    cmd = sprintf('%s -n -get %s  -form %s',exe,url,form);
else
    cmd = sprintf('%s -n -source -get %s',exe,url);
end

% translation
cmd=strrep(cmd,'\','\\');			       
cmd=strrep(cmd,'!','\!');			       

% appel de w3c
[status,reponse]= unix(cmd);
  
  
