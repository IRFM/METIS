% ZUICHGFCT	formulaire de saisie nom fonction des modules externes
%--------------------------------------------------------------
% fichier  ->  zuichgfct
%					zuicreeform, gmarge
%
%
% fonction Matlab 5 :
%	fonction de creation de GUI pour le formulaire
% 	de saisie du nom de la fonction
%	des modules externes , parametres  g��aux
% 	sous le mode edition du formulaire principal
% 
% syntaxe  :
%	hout = zuichgfct(module,liste_fct,hfct,hmodule) 
%
% entrees :
%	module : nom du module
%	liste_fct : liste des fonctions associ�s
%	hfct : handle du bouton "nom du module" pour la mise a jour en retour
%	hmodule : handles du bouton "texte de la fonction" pour la mise a jour en retour
%
% sorties :
%     hout  =  handle de la fenetre cree.
% 
% fonction �rite par formulaire  , poste XX-XX  
% version  1.7  du  29/09/2001
%
% liste des modifications : 
%
%--------------------------------------------------------------
function hout = zuichgfct(module,liste_fct,hfct,hmodule)

if nargin < 1
	return
end
if nargin < 2 | isempty(liste_fct)
	liste_fct = ' ' ;
end
% 
var = strcat('param.fonction.',module) ;
fct = evalin('base',var) ;
ind = strmatch(fct,liste_fct,'exact') ;
if isempty(ind)
	if ~exist(fct)
		fct = '';

	end
	ind = 1;
elseif isempty(fct)
	fct = '';
	ind = 1;
end
ind =ind(1);

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle([module 'chgfct']) ;
if ishandle(hform)
	%zuiformvisible(hform) ;
 	%zuidata(hui.popup_nom_fct,zuidata(hui.edit_nom_fct),'zuimax') ;
	%zuicloseone(hform) ;
	%hout = zuichgfct(module,liste_fct,hfct,hmodule) ;
	%return
	close(hform);
end

% le formulaire
form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

% 1ere ligne nom de la machine
col1 = {'titre','text@full',sprintf(' module %s ',module),[],''};
form{length(form)+1} = {col1};
col1 = {'text_nom_fct','text','function name',[],''};
col2 = {'popup_nom_fct','popup',liste_fct,ind,'associated function name',liste_fct,'void'};
col3 = {'edit_nom_fct','edit',fct,40,'Ohter function',[],var};
form{length(form)+1} = {col1,col2,col3};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
%celldisp(form)

hout=zuicreeform('Modules externes',[module 'chgfct'],'zuichgfct_fct','zuichgfct_ctrl',form) ;
setappdata(hout,'hfct',hfct) ;
setappdata(hout,'hmodule',hmodule) ;
setappdata(hout,'module',module) ;
