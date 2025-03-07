% ZUIDEDIT_DATA_SRCN0_FCT  gestion callbacks du formulaired'edition de profils de source-n0
%--------------------------------------------------------------------
% fichier zuiedit_data_srcn0_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de profils de source-n0
%
% syntaxe :
% 
% entrees :
%  action =  tag du uicontrol active
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
%
%--------------------------------------------------------------
function zuiedit_data_srcn0_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_srcn0') ;

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
texte_prop = 'n0' ;
var_prop   = 'data.cons.n0' ;

% selon ation
switch lower(action)

case 'edit_el'
	nom        = 'data.source.n0.el' ;
	y          = evalin('base','data.source.n0.el') ;
	texte_y    = 'n0.el' ;
	var_y      = 'data.source.n0.el';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_el) ;

case 'import_el'
	nom     = 'data.source.n0.el' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_el) ;
	
case 'dessin_el'
	hout = zuiedit_profcmplx('data.source.n0.el','data.gene.temps','param.gene.x')
	zuireset(h.dessin_el) ;
	

case 'edit_ion'
	nom        = 'data.source.n0.ion' ;
	y          = evalin('base','data.source.n0.ion') ;
	texte_y    = 'n0.ion' ;
	var_y      = 'data.source.n0.ion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ion) ;

case 'import_ion'
	nom     = 'data.source.n0.ion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ion) ;

case 'dessin_ion'
	hout = zuiedit_profcmplx('data.source.n0.ion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ion) ;
	
case 'edit_ne'
	nom        = 'data.source.n0.ne' ;
	y          = evalin('base','data.source.n0.ne') ;
	texte_y    = 'n0.ne' ;
	var_y      = 'data.source.n0.ne';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ne) ;

case 'import_ne'
	nom     = 'data.source.n0.ne' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ne) ;

case 'dessin_ne'
	hout = zuiedit_profcmplx('data.source.n0.ne','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ne) ;
	
case 'edit_j'
	nom        = 'data.source.n0.j' ;
	y          = evalin('base','data.source.n0.j') ;
	texte_y    = 'n0.j' ;
	var_y      = 'data.source.n0.j';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_j) ;

case 'import_j'
	nom     = 'data.source.n0.j' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_j) ;

case 'dessin_j'
	hout = zuiedit_profcmplx('data.source.n0.j','data.gene.temps','param.gene.x')
	zuireset(h.dessin_j) ;
	
case 'edit_w'
	nom        = 'data.source.n0.w' ;
	y          = evalin('base','data.source.n0.w') ;
	texte_y    = 'n0.w' ;
	var_y      = 'data.source.n0.w';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_w) ;

case 'import_w'
	nom     = 'data.source.n0.w' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_w) ;

case 'dessin_w'
	hout = zuiedit_profcmplx('data.source.n0.w','data.gene.temps','param.gene.x')
	zuireset(h.dessin_w) ;
	
case 'edit_wb'
	nom        = 'data.source.n0.wb' ;
	y          = evalin('base','data.source.n0.wb') ;
	texte_y    = 'n0.wb' ;
	var_y      = 'data.source.n0.wb';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_wb) ;

case 'import_wb'
	nom     = 'data.source.n0.wb' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_wb) ;

case 'dessin_wb'
	hout = zuiedit_profcmplx('data.source.n0.wb','data.gene.temps','param.gene.x')
	zuireset(h.dessin_wb) ;
	
case 'edit_q'
	nom        = 'data.source.n0.q' ;
	y          = evalin('base','data.source.n0.q') ;
	texte_y    = 'n0.q' ;
	var_y      = 'data.source.n0.q';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_q) ;

case 'import_q'
	nom     = 'data.source.n0.q' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_q) ;

case 'dessin_q'
	hout = zuiedit_profcmplx('data.source.n0.q','data.gene.temps','param.gene.x')
	zuireset(h.dessin_q) ;
	
case 'edit_fluce'
	nom        = 'data.source.n0.fluce' ;
	y          = evalin('base','data.source.n0.fluce') ;
	texte_y    = 'n0.fluce' ;
	var_y      = 'data.source.n0.fluce';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_fluce) ;

case 'import_fluce'
	nom     = 'data.source.n0.fluce' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_fluce) ;

case 'dessin_fluce'
	hout = zuiedit_profcmplx('data.source.n0.fluce','data.gene.temps','param.gene.x')
	zuireset(h.dessin_fluce) ;


case 'edit_flucion'
	nom        = 'data.source.n0.flucion' ;
	y          = evalin('base','data.source.n0.flucion') ;
	texte_y    = 'n0.flucion' ;
	var_y      = 'data.source.n0.flucion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_flucion) ;

case 'import_flucion'
	nom     = 'data.source.n0.flucion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_flucion) ;

case 'dessin_flucion'
	hout = zuiedit_profcmplx('data.source.n0.flucion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_flucion) ;

case 'edit_psupra'
	nom        = 'data.source.n0.psupra' ;
	y          = evalin('base','data.source.n0.psupra') ;
	texte_y    = 'n0.psupra' ;
	var_y      = 'data.source.n0.psupra' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_psupra) ;

case 'import_psupra'
	nom     = 'data.source.n0.psupra' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_psupra) ;

case 'dessin_psupra'
	hout = zuiedit_profcmplx('data.source.n0.psupra','data.gene.temps','param.gene.x')
	zuireset(h.dessin_psupra) ;


case 'edit_paniso'
	nom        = 'data.source.n0.paniso' ;
	y          = evalin('base','data.source.n0.paniso') ;
	texte_y    = 'n0.paniso' ;
	var_y      = 'data.source.n0.paniso';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_paniso) ;

case 'import_paniso'
	nom     = 'data.source.n0.paniso' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_paniso) ;

case 'dessin_paniso'
	hout = zuiedit_profcmplx('data.source.n0.paniso','data.gene.temps','param.gene.x')
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

