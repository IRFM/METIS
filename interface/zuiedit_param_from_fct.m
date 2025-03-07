% ZUIDEDIT_PARAM_FROM_FCT  gestion des callbacks du formulaire d'informations (from) 
%--------------------------------------------------------------
% fichier zuiedit_param_from_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	d'informations (from) , parametres  généraux
%	sous le mode edition du formulaire principal
%
% syntaxe :
% 
% entrees :
%  action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 3.0 du 17/01/2005.
% 
% liste des modifications : 
%  * 28/08/2001 -> ajout du callback d'affichage des infos du createur  (J-F Artaud)
%  * 27/09/2001 -> rajout de "zuidownloadform" dans la partie "Validation"
%  * 17/01/2005 -> remplace rm par rm -f
%--------------------------------------------------------------
function zuiedit_param_from_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('from');
% information pour l'assistant
zuicr(hfig,action) ;
hoc = getfield(h,action) ;

% selon ation
switch lower(action)

% fonctions equilibre
case {'radio_info'}
	disp('Ouverture du fichier sous nedit') ;
	tmp = zfrominfo(1) ;
	[s,t] = unix(sprintf('%s %s &',getappdata(0,'editeur'),tmp));
	if s~=0
		herror = errordlg('à l''ouverture du fichier ??','PROBLEME') ;
		zwaitfor(herror) ;
		zuicloseone(hform) ;
	end
	% eval(['[s,t] = unix(''rm -f ' tmp ''');']) ;		
	zuireset(h.radio_info)

case {'radio_creainfo'}
	disp('Ouverture du fichier sous nedit') ;
	tmp = zfrominfo(2) ;
	[s,t] = unix(sprintf('%s %s &',getappdata(0,'editeur'),tmp));
	if s~=0
		herror = errordlg('à l''ouverture du fichier ??','PROBLEME') ;
		zwaitfor(herror) ;
		zuicloseone(hform) ;
	end
	eval(['[s,t] = unix(''rm -f ' tmp ''');']) ;		
	zuireset(h.radio_creainfo)

case {'annulation','close'}
	zuicloseone(hfig);	
	
case {'init','raz'}
	zuiuploadform(hfig)
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuireset(h.raz);
	
case 'radio_createur'
	% ajout du 28/08/2001 (J-F Artaud) 
	hcr=zuicreefunform('ztsacces','param.from.option',1,1);
	zuireset(h.radio_createur);

case 'validation'
	zuidownloadform(hfig) ;
	zuireset(h.validation) ;
	zuicloseone(hfig) ;
	
otherwise
	warning('ation non prise en compte')
	   
end

