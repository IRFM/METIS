% ZUIDEDIT_DATA_SRCFCI_FCT gestion callbacks du formulaired'edition de profils de source-fci
%--------------------------------------------------------------------
% fichier zuiedit_data_srcfci_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de profils de source-fci
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
function zuiedit_data_srcfci_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_srcfci') ;

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
texte_prop = 'fci' ;
var_prop   = 'data.cons.fci' ;

% selon ation
switch lower(action)

case 'edit_el'
	nom        = 'data.source.fci.el' ;
	y          = evalin('base','data.source.fci.el') ;
	texte_y    = 'fci.el' ;
	var_y      = 'data.source.fci.el';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_el) ;

case 'import_el'
	nom     = 'data.source.fci.el' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_el) ;
	
case 'dessin_el'
	hout = zuiedit_profcmplx('data.source.fci.el','data.gene.temps','param.gene.x')
	zuireset(h.dessin_el) ;
	

case 'edit_ion'
	nom        = 'data.source.fci.ion' ;
	y          = evalin('base','data.source.fci.ion') ;
	texte_y    = 'fci.ion' ;
	var_y      = 'data.source.fci.ion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ion) ;

case 'import_ion'
	nom     = 'data.source.fci.ion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ion) ;

case 'dessin_ion'
	hout = zuiedit_profcmplx('data.source.fci.ion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ion) ;
	
case 'edit_ne'
	nom        = 'data.source.fci.ne' ;
	y          = evalin('base','data.source.fci.ne') ;
	texte_y    = 'fci.ne' ;
	var_y      = 'data.source.fci.ne';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ne) ;

case 'import_ne'
	nom     = 'data.source.fci.ne' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ne) ;

case 'dessin_ne'
	hout = zuiedit_profcmplx('data.source.fci.ne','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ne) ;
	
case 'edit_j'
	nom        = 'data.source.fci.j' ;
	y          = evalin('base','data.source.fci.j') ;
	texte_y    = 'fci.j' ;
	var_y      = 'data.source.fci.j';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_j) ;

case 'import_j'
	nom     = 'data.source.fci.j' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_j) ;

case 'dessin_j'
	hout = zuiedit_profcmplx('data.source.fci.j','data.gene.temps','param.gene.x')
	zuireset(h.dessin_j) ;
	
case 'edit_w'
	nom        = 'data.source.fci.w' ;
	y          = evalin('base','data.source.fci.w') ;
	texte_y    = 'fci.w' ;
	var_y      = 'data.source.fci.w';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_w) ;

case 'import_w'
	nom     = 'data.source.fci.w' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_w) ;

case 'dessin_w'
	hout = zuiedit_profcmplx('data.source.fci.w','data.gene.temps','param.gene.x')
	zuireset(h.dessin_w) ;
	
case 'edit_wb'
	nom        = 'data.source.fci.wb' ;
	y          = evalin('base','data.source.fci.wb') ;
	texte_y    = 'fci.wb' ;
	var_y      = 'data.source.fci.wb';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_wb) ;

case 'import_wb'
	nom     = 'data.source.fci.wb' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_wb) ;

case 'dessin_wb'
	hout = zuiedit_profcmplx('data.source.fci.wb','data.gene.temps','param.gene.x')
	zuireset(h.dessin_wb) ;
	
case 'edit_q'
	nom        = 'data.source.fci.q' ;
	y          = evalin('base','data.source.fci.q') ;
	texte_y    = 'fci.q' ;
	var_y      = 'data.source.fci.q';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_q) ;

case 'import_q'
	nom     = 'data.source.fci.q' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_q) ;

case 'dessin_q'
	hout = zuiedit_profcmplx('data.source.fci.q','data.gene.temps','param.gene.x')
	zuireset(h.dessin_q) ;
	
case 'edit_fluce'
	nom        = 'data.source.fci.fluce' ;
	y          = evalin('base','data.source.fci.fluce') ;
	texte_y    = 'fci.fluce' ;
	var_y      = 'data.source.fci.fluce';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_fluce) ;

case 'import_fluce'
	nom     = 'data.source.fci.fluce' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_fluce) ;

case 'dessin_fluce'
	hout = zuiedit_profcmplx('data.source.fci.fluce','data.gene.temps','param.gene.x')
	zuireset(h.dessin_fluce) ;


case 'edit_flucion'
	nom        = 'data.source.fci.flucion' ;
	y          = evalin('base','data.source.fci.flucion') ;
	texte_y    = 'fci.flucion' ;
	var_y      = 'data.source.fci.flucion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_flucion) ;

case 'import_flucion'
	nom     = 'data.source.fci.flucion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_flucion) ;

case 'dessin_flucion'
	hout = zuiedit_profcmplx('data.source.fci.flucion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_flucion) ;


case 'edit_psupra'
	nom        = 'data.source.fci.psupra' ;
	y          = evalin('base','data.source.fci.psupra') ;
	texte_y    = 'fci.psupra' ;
	var_y      = 'data.source.fci.psupra' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_psupra) ;

case 'import_psupra'
	nom     = 'data.source.fci.psupra' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_psupra) ;

case 'dessin_psupra'
	hout = zuiedit_profcmplx('data.source.fci.psupra','data.gene.temps','param.gene.x')
	zuireset(h.dessin_psupra) ;


case 'edit_paniso'
	nom        = 'data.source.fci.paniso' ;
	y          = evalin('base','data.source.fci.paniso') ;
	texte_y    = 'fci.paniso' ;
	var_y      = 'data.source.fci.paniso';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_paniso) ;

case 'import_paniso'
	nom     = 'data.source.fci.paniso' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_paniso) ;

case 'dessin_paniso'
	hout = zuiedit_profcmplx('data.source.fci.paniso','data.gene.temps','param.gene.x')
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

