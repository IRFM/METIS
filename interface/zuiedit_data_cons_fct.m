% ZUIDEDIT_DATA_CONS_FCT  gestion des callbacksdu formulaire d'edition de consignes
%--------------------------------------------------------------
% fichier zuiedit_data_cons_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	d'edition de consignes
%
% syntaxe :
%	zuiedit_data_cons_fct(action)
%
% entrees : 
%	action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 31/07/2001.
% 
% liste des modifications : 
% * 14/03/2002 -> ajout de la stabilite MHD
%
%--------------------------------------------------------------
function zuiedit_data_cons_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info = zinfo ;
% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_cons') ;

% information pour l'assistant
zuicr(hfig,action) ;

% selon ation
switch lower(action)

case 'limites'
	zuiedit_data_cons_limites ;
	zuireset(h.limites) ;

case 'chauffages'
	zuiedit_data_cons_chauf ;
	zuireset(h.chauffages) ;

case 'inject'
	zuiedit_data_cons_inj ;
	zuireset(h.inject) ;

case 'mhd'
	% zuiedit_data_cons_mhd ;
	nom         = 'data.cons.stab' ;
	x           = evalin('base','data.gene.temps') ;
	y           = evalin('base',nom) ;
	texte_x     = 'temps' ;
	texte_y     = 'stab'
	var_x       = 'void';
	var_y       = nom ;
	canal       = 1 ;
	code_retour = '' ;
	liste_ref   = {} ;
	var_ref     = {} ;
	texte_prop  = '' ;
	var_prop    = '' ;
 	hout = zuieditcons(nom,info.data.cons.stab,x,y,texte_x,texte_y,var_x,var_y, ...
 	                   canal,code_retour,liste_ref,var_ref,texte_prop,var_prop)  ;
	zuireset(h.mhd) ;

case {'btn_quit','close'}
	zuicloseone(hfig);	
	
case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
otherwise
	warning('ation non prise en compte')
	   
end

