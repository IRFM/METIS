% ZUIDEDIT_DATA_SRCEXT_FCT  gestion  callbacks du formulaire d'edition de profils de source-ext
%--------------------------------------------------------------------
% fichier zuiedit_data_srcext_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de profils de source-ext
%
% syntaxe :
% 
% entrees :
%	action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 09/03/2003.
% 
% liste des modifications : 
%
% * 27/08/2003 -> ajout de kini pour edit  
% * 02/09/2003 -> ajout des editions de source.??.q 
%--------------------------------------------------------------
function zuiedit_data_srcext_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_srcext') ;

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
texte_prop = 'ext' ;
var_prop   = 'data.cons.ext' ;

% selon ation
switch lower(action)

case 'edit_el'
	nom        = 'data.source.ext.el' ;
	y          = evalin('base','data.source.ext.el') ;
	texte_y    = 'ext.el' ;
	var_y      = 'data.source.ext.el';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_el) ;

case 'import_el'
	nom     = 'data.source.ext.el' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_el) ;
	
case 'dessin_el'
	hout = zuiedit_profcmplx('data.source.ext.el','data.gene.temps','param.gene.x')
	zuireset(h.dessin_el) ;
	

case 'edit_ion'
	nom        = 'data.source.ext.ion' ;
	y          = evalin('base','data.source.ext.ion') ;
	texte_y    = 'ext.ion' ;
	var_y      = 'data.source.ext.ion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ion) ;

case 'import_ion'
	nom     = 'data.source.ext.ion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ion) ;

case 'dessin_ion'
	hout = zuiedit_profcmplx('data.source.ext.ion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ion) ;
	
case 'edit_ne'
	nom        = 'data.source.ext.ne' ;
	y          = evalin('base','data.source.ext.ne') ;
	texte_y    = 'ext.ne' ;
	var_y      = 'data.source.ext.ne';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ne) ;

case 'import_ne'
	nom     = 'data.source.ext.ne' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ne) ;

case 'dessin_ne'
	hout = zuiedit_profcmplx('data.source.ext.ne','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ne) ;

case 'edit_j'
	nom        = 'data.source.ext.j' ;
	y          = evalin('base','data.source.ext.j') ;
	texte_y    = 'ext.j' ;
	var_y      = 'data.source.ext.j';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_j) ;
	
case 'import_j'
	nom     = 'data.source.ext.j' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_j) ;

case 'dessin_j'
	hout = zuiedit_profcmplx('data.source.ext.j','data.gene.temps','param.gene.x')
	zuireset(h.dessin_j) ;
	
case 'edit_w'
	nom        = 'data.source.ext.w' ;
	y          = evalin('base','data.source.ext.w') ;
	texte_y    = 'ext.w' ;
	var_y      = 'data.source.ext.w';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_w) ;

case 'import_w'
	nom     = 'data.source.ext.w' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_w) ;

case 'dessin_w'
	hout = zuiedit_profcmplx('data.source.ext.w','data.gene.temps','param.gene.x')
	zuireset(h.dessin_w) ;

case 'edit_wb'
	nom        = 'data.source.ext.wb' ;
	y          = evalin('base','data.source.ext.wb') ;
	texte_y    = 'ext.wb' ;
	var_y      = 'data.source.ext.wb';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_wb) ;

case 'import_wb'
	nom     = 'data.source.ext.wb' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_wb) ;

case 'dessin_wb'
	hout = zuiedit_profcmplx('data.source.ext.wb','data.gene.temps','param.gene.x')
	zuireset(h.dessin_wb) ;

case 'edit_q'
	nom        = 'data.source.ext.q' ;
	y          = evalin('base','data.source.ext.q') ;
	texte_y    = 'ext.q' ;
	var_y      = 'data.source.ext.q';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_q) ;

case 'import_q'
	nom     = 'data.source.ext.q' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_q) ;

case 'dessin_q'
	hout = zuiedit_profcmplx('data.source.ext.q','data.gene.temps','param.gene.x')
	zuireset(h.dessin_q) ;

case 'edit_fluce'
	nom        = 'data.source.ext.fluce' ;
	y          = evalin('base','data.source.ext.fluce') ;
	texte_y    = 'ext.fluce' ;
	var_y      = 'data.source.ext.fluce';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_fluce) ;

case 'import_fluce'
	nom     = 'data.source.ext.fluce' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_fluce) ;

case 'dessin_fluce'
	hout = zuiedit_profcmplx('data.source.ext.fluce','data.gene.temps','param.gene.x')
	zuireset(h.dessin_fluce) ;


case 'edit_flucion'
	nom        = 'data.source.ext.flucion' ;
	y          = evalin('base','data.source.ext.flucion') ;
	texte_y    = 'ext.flucion' ;
	var_y      = 'data.source.ext.flucion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_flucion) ;

case 'import_flucion'
	nom     = 'data.source.ext.flucion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_flucion) ;

case 'dessin_flucion'
	hout = zuiedit_profcmplx('data.source.ext.flucion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_flucion) ;

case 'edit_psupra'
	nom        = 'data.source.ext.psupra' ;
	y          = evalin('base','data.source.ext.psupra') ;
	texte_y    = 'ext.psupra' ;
	var_y      = 'data.source.ext.psupra' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_psupra) ;

case 'import_psupra'
	nom     = 'data.source.ext.psupra' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_psupra) ;

case 'dessin_psupra'
	hout = zuiedit_profcmplx('data.source.ext.psupra','data.gene.temps','param.gene.x')
	zuireset(h.dessin_psupra) ;


case 'edit_paniso'
	nom        = 'data.source.ext.paniso' ;
	y          = evalin('base','data.source.ext.paniso') ;
	texte_y    = 'ext.paniso' ;
	var_y      = 'data.source.ext.paniso';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_paniso) ;

case 'import_paniso'
	nom     = 'data.source.ext.paniso' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_paniso) ;

case 'dessin_paniso'
	hout = zuiedit_profcmplx('data.source.ext.paniso','data.gene.temps','param.gene.x')
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

