% ZUIDEDIT_DATA_SRCFCE_FCT gestion callbacks du formulaired'edition de profils de source-fce
%--------------------------------------------------------------------
% fichier zuiedit_data_srcfce_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de profils de source-fce
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
% * 02/09/2003 -> ajout des editions de source.??.q %
%--------------------------------------------------------------
function zuiedit_data_srcfce_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_srcfce') ;

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
texte_prop = 'fce' ;
var_prop   = 'data.cons.fce' ;

% selon ation
switch lower(action)

case 'edit_el'
	nom        = 'data.source.fce.el' ;
	y          = evalin('base','data.source.fce.el') ;
	texte_y    = 'fce.el' ;
	var_y      = 'data.source.fce.el';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_el) ;

case 'import_el'
	nom     = 'data.source.fce.el' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_el) ;
	
case 'dessin_el'
	hout = zuiedit_profcmplx('data.source.fce.el','data.gene.temps','param.gene.x')
	zuireset(h.dessin_el) ;
	
case 'edit_ion'
	nom        = 'data.source.fce.ion' ;
	y          = evalin('base','data.source.fce.ion') ;
	texte_y    = 'fce.ion' ;
	var_y      = 'data.source.fce.ion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ion) ;

case 'import_ion'
	nom     = 'data.source.fce.ion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ion) ;

case 'dessin_ion'
	hout = zuiedit_profcmplx('data.source.fce.ion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ion) ;
	
case 'edit_ne'
	nom        = 'data.source.fce.ne' ;
	y          = evalin('base','data.source.fce.ne') ;
	texte_y    = 'fce.ne' ;
	var_y      = 'data.source.fce.ne';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ne) ;

case 'import_ne'
	nom     = 'data.source.fce.ne' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ne) ;

case 'dessin_ne'
	hout = zuiedit_profcmplx('data.source.fce.ne','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ne) ;
	
case 'edit_j'
	nom        = 'data.source.fce.j' ;
	y          = evalin('base','data.source.fce.j') ;
	texte_y    = 'fce.j' ;
	var_y      = 'data.source.fce.j';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_j) ;

case 'import_j'
	nom     = 'data.source.fce.j' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_j) ;

case 'dessin_j'
	hout = zuiedit_profcmplx('data.source.fce.j','data.gene.temps','param.gene.x')
	zuireset(h.dessin_j) ;
	
case 'edit_w'
	nom        = 'data.source.fce.w' ;
	y          = evalin('base','data.source.fce.w') ;
	texte_y    = 'fce.w' ;
	var_y      = 'data.source.fce.w';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_w) ;

case 'import_w'
	nom     = 'data.source.fce.w' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_w) ;

case 'dessin_w'
	hout = zuiedit_profcmplx('data.source.fce.w','data.gene.temps','param.gene.x')
	zuireset(h.dessin_w) ;
	
case 'edit_wb'
	nom        = 'data.source.fce.wb' ;
	y          = evalin('base','data.source.fce.wb') ;
	texte_y    = 'fce.wb' ;
	var_y      = 'data.source.fce.wb';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_wb) ;

case 'import_wb'
	nom     = 'data.source.fce.wb' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_wb) ;

case 'dessin_wb'
	hout = zuiedit_profcmplx('data.source.fce.wb','data.gene.temps','param.gene.x')
	zuireset(h.dessin_wb) ;
	
case 'edit_q'
	nom        = 'data.source.fce.q' ;
	y          = evalin('base','data.source.fce.q') ;
	texte_y    = 'fce.q' ;
	var_y      = 'data.source.fce.q';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_q) ;

case 'import_q'
	nom     = 'data.source.fce.q' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_q) ;

case 'dessin_q'
	hout = zuiedit_profcmplx('data.source.fce.q','data.gene.temps','param.gene.x')
	zuireset(h.dessin_q) ;
	
case 'edit_fluce'
	nom        = 'data.source.fce.fluce' ;
	y          = evalin('base','data.source.fce.fluce') ;
	texte_y    = 'fce.fluce' ;
	var_y      = 'data.source.fce.fluce';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_fluce) ;

case 'import_fluce'
	nom     = 'data.source.fce.fluce' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_fluce) ;

case 'dessin_fluce'
	hout = zuiedit_profcmplx('data.source.fce.fluce','data.gene.temps','param.gene.x')
	zuireset(h.dessin_fluce) ;
	
case 'edit_flucion'
	nom        = 'data.source.fce.flucion' ;
	y          = evalin('base','data.source.fce.flucion') ;
	texte_y    = 'fce.flucion' ;
	var_y      = 'data.source.fce.flucion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_flucion) ;

case 'import_flucion'
	nom     = 'data.source.fce.flucion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_flucion) ;

case 'dessin_flucion'
	hout = zuiedit_profcmplx('data.source.fce.flucion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_flucion) ;
	
case 'edit_psupra'
	nom        = 'data.source.fce.psupra' ;
	y          = evalin('base','data.source.fce.psupra') ;
	texte_y    = 'fce.psupra' ;
	var_y      = 'data.source.fce.psupra' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_psupra) ;

case 'import_psupra'
	nom     = 'data.source.fce.psupra' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_psupra) ;

case 'dessin_psupra'
	hout = zuiedit_profcmplx('data.source.fce.psupra','data.gene.temps','param.gene.x')
	zuireset(h.dessin_psupra) ;
	
case 'edit_paniso'
	nom        = 'data.source.fce.paniso' ;
	y          = evalin('base','data.source.fce.paniso') ;
	texte_y    = 'fce.paniso' ;
	var_y      = 'data.source.fce.paniso';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_paniso) ;

case 'import_paniso'
	nom     = 'data.source.fce.paniso' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_paniso) ;

case 'dessin_paniso'
	hout = zuiedit_profcmplx('data.source.fce.paniso','data.gene.temps','param.gene.x')
	zuireset(h.dessin_paniso') ;
	
case {'btn_quit','close'}
	zuicloseone(hfig);	
	
case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
otherwise
	warning('ation non prise en compte')
	   
end

