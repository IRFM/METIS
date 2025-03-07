% ZUIDEDIT_DATA_SRCIDN_FCT gestion callbacks du formulaire d'edition de profils de source-idn
%--------------------------------------------------------------------
% fichier zuiedit_data_srcidn_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de profils de source-idn
%
% syntaxe :
% 
% entrees :
%	action       =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 02/09/2003.
% 
% liste des modifications : 
%
% * 27/08/2003 -> ajout de kini pour edit  
% * 02/09/2003 -> ajout des editions de source.??.q 
%
%--------------------------------------------------------------
function zuiedit_data_srcidn_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_srcidn') ;

% information pour l'assistant
zuicr(hfig,action) ;

% variables d'entrée de l'éditeur de consignes zuieditcons
x           = evalin('base','param.gene.x') ;
kini        = evalin('base','param.gene.kmin') ;
texte_x     = 'x' ;
var_x       = 'void';
canal       = 1 ;
code_retour = '' ;
liste_ref   = {'Pe','Pion','Ne','Psi','Jmoy','\alphae','Zeff','vide'} ;
var_ref     = {  {'param.gene.x','data.prof.pe','-'},
                 {'param.gene.x','data.prof.pion','-'},
                 {'param.gene.x','data.prof.ne','-'},
                 {'param.gene.x','data.prof.psi','-'},
                 {'param.gene.x','data.prof.jmoy','-'},
                 {'param.gene.x','data.prof.ae','-'},
                 {'param.gene.x','data.prof.zeff','-'},
		 {'[]','[]',''}  } ;
texte_prop = 'idn' ;
var_prop   = 'data.cons.idn' ;

% selon ation
switch lower(action)

case 'edit_el'
	nom        = 'data.source.idn.el' ;
	y          = evalin('base','data.source.idn.el') ;
	texte_y    = 'idn.el' ;
	var_y      = 'data.source.idn.el';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_el) ;

case 'import_el'
	nom     = 'data.source.idn.el' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_el) ;
	
case 'dessin_el'
	hout = zuiedit_profcmplx('data.source.idn.el','data.gene.temps','param.gene.x')
	zuireset(h.dessin_el) ;
	

case 'edit_ion'
	nom        = 'data.source.idn.ion' ;
	y          = evalin('base','data.source.idn.ion') ;
	texte_y    = 'idn.ion' ;
	var_y      = 'data.source.idn.ion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ion) ;

case 'import_ion'
	nom     = 'data.source.idn.ion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ion) ;

case 'dessin_ion'
	hout = zuiedit_profcmplx('data.source.idn.ion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ion) ;
	
case 'edit_ne'
	nom        = 'data.source.idn.ne' ;
	y          = evalin('base','data.source.idn.ne') ;
	texte_y    = 'idn.ne' ;
	var_y      = 'data.source.idn.ne';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ne) ;

case 'import_ne'
	nom     = 'data.source.idn.ne' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ne) ;

case 'dessin_ne'
	hout = zuiedit_profcmplx('data.source.idn.ne','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ne) ;
	
case 'edit_j'
	nom        = 'data.source.idn.j' ;
	y          = evalin('base','data.source.idn.j') ;
	texte_y    = 'idn.j' ;
	var_y      = 'data.source.idn.j';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_j) ;

case 'import_j'
	nom     = 'data.source.idn.j' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_j) ;

case 'dessin_j'
	hout = zuiedit_profcmplx('data.source.idn.j','data.gene.temps','param.gene.x')
	zuireset(h.dessin_j) ;
	
case 'edit_w'
	nom        = 'data.source.idn.w' ;
	y          = evalin('base','data.source.idn.w') ;
	texte_y    = 'idn.w' ;
	var_y      = 'data.source.idn.w';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_w) ;

case 'import_w'
	nom     = 'data.source.idn.w' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_w) ;

case 'dessin_w'
	hout = zuiedit_profcmplx('data.source.idn.w','data.gene.temps','param.gene.x')
	zuireset(h.dessin_w) ;
	
case 'edit_wb'
	nom        = 'data.source.idn.wb' ;
	y          = evalin('base','data.source.idn.wb') ;
	texte_y    = 'idn.wb' ;
	var_y      = 'data.source.idn.wb';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_wb) ;

case 'import_wb'
	nom     = 'data.source.idn.wb' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_wb) ;

case 'dessin_wb'
	hout = zuiedit_profcmplx('data.source.idn.wb','data.gene.temps','param.gene.x')
	zuireset(h.dessin_wb) ;
	
case 'edit_q'
	nom        = 'data.source.idn.q' ;
	y          = evalin('base','data.source.idn.q') ;
	texte_y    = 'idn.q' ;
	var_y      = 'data.source.idn.q';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_q) ;

case 'import_q'
	nom     = 'data.source.idn.q' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_q) ;

case 'dessin_q'
	hout = zuiedit_profcmplx('data.source.idn.q','data.gene.temps','param.gene.x')
	zuireset(h.dessin_q) ;
	
case 'edit_fluce'
	nom        = 'data.source.idn.fluce' ;
	y          = evalin('base','data.source.idn.fluce') ;
	texte_y    = 'idn.fluce' ;
	var_y      = 'data.source.idn.fluce';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_fluce) ;

case 'import_fluce'
	nom     = 'data.source.idn.fluce' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_fluce) ;

case 'dessin_fluce'
	hout = zuiedit_profcmplx('data.source.idn.fluce','data.gene.temps','param.gene.x')
	zuireset(h.dessin_fluce) ;


case 'edit_flucion'
	nom        = 'data.source.idn.flucion' ;
	y          = evalin('base','data.source.idn.flucion') ;
	texte_y    = 'idn.flucion' ;
	var_y      = 'data.source.idn.flucion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_flucion) ;

case 'import_flucion'
	nom     = 'data.source.idn.flucion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_flucion) ;

case 'dessin_flucion'
	hout = zuiedit_profcmplx('data.source.idn.flucion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_flucion) ;

case 'edit_psupra'
	nom        = 'data.source.idn.psupra' ;
	y          = evalin('base','data.source.idn.psupra') ;
	texte_y    = 'idn.psupra' ;
	var_y      = 'data.source.idn.psupra' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_psupra) ;

case 'import_psupra'
	nom     = 'data.source.idn.psupra' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_psupra) ;

case 'dessin_psupra'
	hout = zuiedit_profcmplx('data.source.idn.psupra','data.gene.temps','param.gene.x')
	zuireset(h.dessin_psupra) ;


case 'edit_paniso'
	nom        = 'data.source.idn.paniso' ;
	y          = evalin('base','data.source.idn.paniso') ;
	texte_y    = 'idn.paniso' ;
	var_y      = 'data.source.idn.paniso';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_paniso) ;

case 'import_paniso'
	nom     = 'data.source.idn.paniso' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_paniso) ;

case 'dessin_paniso'
	hout = zuiedit_profcmplx('data.source.idn.paniso','data.gene.temps','param.gene.x')
	zuireset(h.dessin_paniso) ;


case {'btn_quit','close'}
	zuicloseone(hfig);	
	
case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
otherwise
	warning('ation non prise en compte')
	   
end

