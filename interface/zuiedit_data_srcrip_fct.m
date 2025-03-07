% ZUIDEDIT_DATA_SRCRIP_FCT gestion callbacks du formulaired'edition de profils de source-ripple
%--------------------------------------------------------------------
% fichier zuiedit_data_srcrip_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de profils de source-rip
%
% syntaxe :
% 
% entrees :
%	action       =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par J-F Artaud, poste 62-15
% version 2.2, du 02/09/2003.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
function zuiedit_data_srcrip_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_srcrip') ;

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
texte_prop = 'ip' ;
var_prop   = 'data.cons.ip' ;

% selon ation
switch lower(action)

case 'edit_el'
	nom        = 'data.source.rip.el' ;
	y          = evalin('base','data.source.rip.el') ;
	texte_y    = 'rip.el' ;
	var_y      = 'data.source.rip.el';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_el) ;

case 'import_el'
	nom     = 'data.source.rip.el' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_el) ;
	
case 'dessin_el'
	hout = zuiedit_profcmplx('data.source.rip.el','data.gene.temps','param.gene.x')
	zuireset(h.dessin_el) ;
	

case 'edit_ion'
	nom        = 'data.source.rip.ion' ;
	y          = evalin('base','data.source.rip.ion') ;
	texte_y    = 'rip.ion' ;
	var_y      = 'data.source.rip.ion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ion) ;

case 'import_ion'
	nom     = 'data.source.rip.ion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ion) ;

case 'dessin_ion'
	hout = zuiedit_profcmplx('data.source.rip.ion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ion) ;
	
case 'edit_ne'
	nom        = 'data.source.rip.ne' ;
	y          = evalin('base','data.source.rip.ne') ;
	texte_y    = 'rip.ne' ;
	var_y      = 'data.source.rip.ne';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ne) ;

case 'import_ne'
	nom     = 'data.source.rip.ne' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ne) ;

case 'dessin_ne'
	hout = zuiedit_profcmplx('data.source.rip.ne','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ne) ;
	
case 'edit_j'
	nom        = 'data.source.rip.j' ;
	y          = evalin('base','data.source.rip.j') ;
	texte_y    = 'rip.j' ;
	var_y      = 'data.source.rip.j';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_j) ;

case 'import_j'
	nom     = 'data.source.rip.j' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_j) ;

case 'dessin_j'
	hout = zuiedit_profcmplx('data.source.rip.j','data.gene.temps','param.gene.x')
	zuireset(h.dessin_j) ;
	
case 'edit_w'
	nom        = 'data.source.rip.w' ;
	y          = evalin('base','data.source.rip.w') ;
	texte_y    = 'rip.w' ;
	var_y      = 'data.source.rip.w';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_w) ;

case 'import_w'
	nom     = 'data.source.rip.w' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_w) ;

case 'dessin_w'
	hout = zuiedit_profcmplx('data.source.rip.w','data.gene.temps','param.gene.x')
	zuireset(h.dessin_w) ;
	
case 'edit_wb'
	nom        = 'data.source.rip.wb' ;
	y          = evalin('base','data.source.rip.wb') ;
	texte_y    = 'rip.wb' ;
	var_y      = 'data.source.rip.wb';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_wb) ;

case 'import_wb'
	nom     = 'data.source.rip.wb' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_wb) ;

case 'dessin_wb'
	hout = zuiedit_profcmplx('data.source.rip.wb','data.gene.temps','param.gene.x')
	zuireset(h.dessin_wb) ;
	
case 'edit_q'
	nom        = 'data.source.rip.q' ;
	y          = evalin('base','data.source.rip.q') ;
	texte_y    = 'rip.q' ;
	var_y      = 'data.source.rip.q';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_q) ;

case 'import_q'
	nom     = 'data.source.rip.q' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_q) ;

case 'dessin_q'
	hout = zuiedit_profcmplx('data.source.rip.q','data.gene.temps','param.gene.x')
	zuireset(h.dessin_q) ;
	
case 'edit_fluce'
	nom        = 'data.source.rip.fluce' ;
	y          = evalin('base','data.source.rip.fluce') ;
	texte_y    = 'rip.fluce' ;
	var_y      = 'data.source.rip.fluce';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_fluce) ;

case 'import_fluce'
	nom     = 'data.source.rip.fluce' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_fluce) ;

case 'dessin_fluce'
	hout = zuiedit_profcmplx('data.source.rip.fluce','data.gene.temps','param.gene.x')
	zuireset(h.dessin_fluce) ;


case 'edit_flucion'
	nom        = 'data.source.rip.flucion' ;
	y          = evalin('base','data.source.rip.flucion') ;
	texte_y    = 'rip.flucion' ;
	var_y      = 'data.source.rip.flucion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_flucion) ;

case 'import_flucion'
	nom     = 'data.source.rip.flucion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_flucion) ;

case 'dessin_flucion'
	hout = zuiedit_profcmplx('data.source.rip.flucion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_flucion) ;


case 'edit_psupra'
	nom        = 'data.source.rip.psupra' ;
	y          = evalin('base','data.source.rip.psupra') ;
	texte_y    = 'rip.psupra' ;
	var_y      = 'data.source.rip.psupra' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_psupra) ;

case 'import_psupra'
	nom     = 'data.source.rip.psupra' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_psupra) ;

case 'dessin_psupra'
	hout = zuiedit_profcmplx('data.source.rip.psupra','data.gene.temps','param.gene.x')
	zuireset(h.dessin_psupra) ;


case 'edit_paniso'
	nom        = 'data.source.rip.paniso' ;
	y          = evalin('base','data.source.rip.paniso') ;
	texte_y    = 'rip.paniso' ;
	var_y      = 'data.source.rip.paniso';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_paniso) ;

case 'import_paniso'
	nom     = 'data.source.rip.paniso' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_paniso) ;

case 'dessin_paniso'
	hout = zuiedit_profcmplx('data.source.rip.paniso','data.gene.temps','param.gene.x')
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

