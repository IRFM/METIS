% ZINTERFACE_CROSSREF  creation des references croisees de l'interface
%----------------------------------------------------------------------
% fichier : zinterface_crossref.m -> zinterface_crossref
% 
% 
% fonction Matlab 5 : 
% 
% fonction pour la creation des references croisees de l'interface
%  
% syntaxe :  
%  function   zinterface_crossref(racine,parent,hfig)
%  
% entrées :  
%  racine  = get(0) avant creation de la nouvelle fenetre
%  parent  = fenetre appelante
%  hfig    = figure courante (qui vient d'etre creee) 
%  
% sorties :  
%   
%  
% fonction écrite par J-F Artaud , poste 46-78
% version  1.7  du  29/09/2001  
%  
%  
% liste des modifications :  
%   * 19/02/2002 -> test sur valeur de uicrossref
%   * 22/03/2002 -> rajout du style de chaque boutons  (Ch. Passeron)
%	 * 26/03/2002 -> set(0,'DefaultFigurePaperPositionMode','auto') avant appel zcapture
%   * 16/05/2002 -> relit toujours les donnees
%-------------------------------------------------------------------------------  
%  
function zinterface_crossref(racine,parent,hfig)

uicrossref = getappdata(0,'uicrossref') ;
if isempty(uicrossref)
   return
elseif uicrossref~=1
	return
end
%if isempty(getappdata(0,'uicrossref'))
%   return
%end

% recherche si la figure a deja ete memorisee
filefig = fullfile(getappdata(0,'root'),'image',strcat(get(hfig,'Tag'),'.png'));
if ~exist(filefig,'file')
	if strcmp(getenv('USER'),'cgc')
		% sauvegarde de la figure
		DefFigMode = get(0,'DefaultFigurePaperPositionMode') ;
		set(0,'DefaultFigurePaperPositionMode','auto') ;
		zcapture(hfig);
		set(0,'DefaultFigurePaperPositionMode',DefFigMode) ;
	end
end

% lecture des donnees
redo = 0;
try 
     interface = evalin('base','interface_crossref');
catch
     redo = 1;
end
if redo == 1    
	ww       = which('zinfo');
	filename = fullfile(fileparts(ww),'zinterface_crossref_data');
	try
	    evalin('base',strcat('load(''',filename,''')'));
	catch
	   zassignin('base','interface_crossref',[]);
	end
	interface = evalin('base','interface_crossref');
end


% info de la figure
data                     = [];
data.name                = get(hfig,'Name');
data.tag                 = get(hfig,'Tag');
if isstruct(parent)
	if strcmp(parent.Type,'figure')
		data.parent.name         = parent.Name;
		data.parent.tag          = parent.Tag;
	else
		data.parent.name         = '';
		data.parent.tag          = '';
	end
else
		data.parent.name         = '';
		data.parent.tag          = '';
end
hc                       = racine.CallbackObject;
if ishandle(hc)
	obj                      = [];
	obj.tag                  = get(hc,'tag');
	obj.value                = get(hc,'value');
	obj.string               = get(hc,'string');
	obj.tooltip              = get(hc,'tooltip');
	obj.userdata             = get(hc,'userdata');
	obj.variable             = getappdata(hc,'variable');
	obj.callback             = get(hc,'callback');
	data.callback_obj        = obj;
else
	obj                      = [];
	obj.tag                  = '';
	obj.value                = [];
	obj.string               = '';
	obj.tooltip              = '';
	obj.userdata             = [];
	obj.variable             = '';
	obj.callback             = '';
	data.callback_obj        = obj;
end

% les variables
hch    = get(hfig,'children');
nosave = 0;
for k = 1:length(hch)
	hc           = hch(k);
	obj          = [];
	varinfo.var  = [];
	obj.tag      = get(hc,'tag');
	obj.value    = get(hc,'value');
	obj.style    = get(hc,'style');
	obj.string   = get(hc,'string');
	obj.tooltip  = get(hc,'tooltip');
	obj.userdata = get(hc,'userdata');
	obj.variable = getappdata(hc,'variable');
	obj.callback = get(hc,'callback');
	
	if ~isempty(obj.variable) | ~isempty(obj.callback)
		data = setfield(data,obj.tag,obj);
	end	
end

% info de la figure
interface = setfield(interface,get(hfig,'tag'),data);
zassignin('base','interface_crossref',interface);
if strcmp(getenv('USER'),'cgc')
	ww       = which('zinfo');
	filename = fullfile(fileparts(ww),'zinterface_crossref_data');
%	evalin('base',strcat('save(''',filename,''',''interface_crossref'',''variable_interface'')'));
	evalin('base',strcat('save(''',filename,''',''interface_crossref'')'));
end
