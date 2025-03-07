% ZUIDEDIT_DATA_SRChyb_FCT  gestion callbacks du formulaire d'edition de profils de source-hyb
%--------------------------------------------------------------------
% fichier zuiedit_data_srchyb_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de profils de source-hyb
%
% syntaxe :
%	zuiedit_data_srchyb_fct(action)
%
% entrees :
%	action =  tag du uicontrol active
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
function zuiedit_data_srchyb_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_srchyb') ;

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
texte_prop = 'hyb' ;
var_prop   = 'data.cons.hyb' ;

% selon ation
switch lower(action)

case 'edit_el'
	nom        = 'data.source.hyb.el' ;
	y          = evalin('base','data.source.hyb.el') ;
	texte_y    = 'hyb.el' ;
	var_y      = 'data.source.hyb.el';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_el) ;

case 'import_el'
	nom     = 'data.source.hyb.el' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_el) ;
	
case 'dessin_el'
	hout = zuiedit_profcmplx('data.source.hyb.el','data.gene.temps','param.gene.x')
	zuireset(h.dessin_el) ;
	

case 'edit_ion'
	nom        = 'data.source.hyb.ion' ;
	y          = evalin('base','data.source.hyb.ion') ;
	texte_y    = 'hyb.ion' ;
	var_y      = 'data.source.hyb.ion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ion) ;

case 'import_ion'
	nom     = 'data.source.hyb.ion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ion) ;

case 'dessin_ion'
	hout = zuiedit_profcmplx('data.source.hyb.ion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ion) ;
	

case 'edit_ne'
	nom        = 'data.source.hyb.ne' ;
	y          = evalin('base','data.source.hyb.ne') ;
	texte_y    = 'hyb.ne' ;
	var_y      = 'data.source.hyb.ne';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ne) ;

case 'import_ne'
	nom     = 'data.source.hyb.ne' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ne) ;

case 'dessin_ne'
	hout = zuiedit_profcmplx('data.source.hyb.ne','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ne) ;
	

case 'edit_j'
	nom        = 'data.source.hyb.j' ;
	y          = evalin('base','data.source.hyb.j') ;
	texte_y    = 'hyb.j' ;
	var_y      = 'data.source.hyb.j';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_j) ;

case 'import_j'
	nom     = 'data.source.hyb.j' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_j) ;

case 'dessin_j'
	hout = zuiedit_profcmplx('data.source.hyb.j','data.gene.temps','param.gene.x')
	zuireset(h.dessin_j) ;
	

case 'edit_w'
	nom        = 'data.source.hyb.w' ;
	y          = evalin('base','data.source.hyb.w') ;
	texte_y    = 'hyb.w' ;
	var_y      = 'data.source.hyb.w';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_w) ;

case 'import_w'
	nom     = 'data.source.hyb.w' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_w) ;

case 'dessin_w'
	hout = zuiedit_profcmplx('data.source.hyb.w','data.gene.temps','param.gene.x')
	zuireset(h.dessin_w) ;
	

case 'edit_wb'
	nom        = 'data.source.hyb.wb' ;
	y          = evalin('base','data.source.hyb.wb') ;
	texte_y    = 'hyb.wb' ;
	var_y      = 'data.source.hyb.wb';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_wb) ;

case 'import_wb'
	nom     = 'data.source.hyb.wb' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_wb) ;

case 'dessin_wb'
	hout = zuiedit_profcmplx('data.source.hyb.wb','data.gene.temps','param.gene.x')
	zuireset(h.dessin_wb) ;
	

case 'edit_q'
	nom        = 'data.source.hyb.q' ;
	y          = evalin('base','data.source.hyb.q') ;
	texte_y    = 'hyb.q' ;
	var_y      = 'data.source.hyb.q';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_q) ;

case 'import_q'
	nom     = 'data.source.hyb.q' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_q) ;

case 'dessin_q'
	hout = zuiedit_profcmplx('data.source.hyb.q','data.gene.temps','param.gene.x')
	zuireset(h.dessin_q) ;
	

case 'edit_fluce'
	nom        = 'data.source.hyb.fluce' ;
	y          = evalin('base','data.source.hyb.fluce') ;
	texte_y    = 'hyb.fluce' ;
	var_y      = 'data.source.hyb.fluce';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_fluce) ;

case 'import_fluce'
	nom     = 'data.source.hyb.fluce' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_fluce) ;

case 'dessin_fluce'
	hout = zuiedit_profcmplx('data.source.hyb.fluce','data.gene.temps','param.gene.x')
	zuireset(h.dessin_fluce) ;
	

case 'edit_flucion'
	nom        = 'data.source.hyb.flucion' ;
	y          = evalin('base','data.source.hyb.flucion') ;
	texte_y    = 'hyb.flucion' ;
	var_y      = 'data.source.hyb.flucion';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_flucion) ;

case 'import_flucion'
	nom     = 'data.source.hyb.flucion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_flucion) ;

case 'dessin_flucion'
	hout = zuiedit_profcmplx('data.source.hyb.flucion','data.gene.temps','param.gene.x')
	zuireset(h.dessin_flucion) ;
	

case 'edit_psupra'
	nom        = 'data.source.hyb.psupra' ;
	y          = evalin('base','data.source.hyb.psupra') ;
	texte_y    = 'hyb.psupra' ;
	var_y      = 'data.source.hyb.psupra' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_psupra) ;

case 'import_psupra'
	nom     = 'data.source.hyb.psupra' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_psupra) ;

case 'dessin_psupra'
	hout = zuiedit_profcmplx('data.source.hyb.psupra','data.gene.temps','param.gene.x')
	zuireset(h.dessin_psupra) ;
	

case 'edit_paniso'
	nom        = 'data.source.hyb.paniso' ;
	y          = evalin('base','data.source.hyb.paniso') ;
	texte_y    = 'hyb.paniso' ;
	var_y      = 'data.source.hyb.paniso';
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_paniso) ;

case 'import_paniso'
	nom     = 'data.source.hyb.paniso' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_paniso) ;

case 'dessin_paniso'
	hout = zuiedit_profcmplx('data.source.hyb.paniso','data.gene.temps','param.gene.x')
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

