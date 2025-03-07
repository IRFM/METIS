% ZUIDEDIT_PARAM_SPLIT_FCT  gestion des callbacks du formulaire de split 
%--------------------------------------------------------------
% fichier zuiedit_param_split_ctrl.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de split , parametres  généraux
% 	sous le mode edition du formulaire principal
%
% syntaxe :
%	zuiedit_param_split_fct(action)
%
% entrees :
%  action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 1.3, du 10/04/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
function zuiedit_param_split_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% disp('callback : ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_param_split');

% information pour l'assistant
zuicr(hfig,action) ;
hoc = getfield(h,action) ;

% selon ation
switch lower(action)

% fonctions equilibre
% 		
	case {'init','raz'}
		zuiformvisible(hfig) ;
		zuiformreset(hfig) ;
		zuiuploadform(hfig) ;
		zuireset(h.raz) ;
	
	case {'annulation','close'}
		zuicloseone(hfig);	
	
	case 'validation'
		zuidownloadform(hfig);
		zuiformcache(hfig) ;
		zuisavenonok;
		zuireset(h.validation);
	
	otherwise
		warning('ation non prise en compte')
	   
end

