% ZUIDEDIT_DATA_COEF_FCT gestion callbacksdu formulaire d'edition de profils des coefficients de transport
%--------------------------------------------------------------------
% fichier zuiedit_data_coef_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de profils des coefficients de transport
%
% syntaxe :
% 
% entrees :
%	action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 02/08/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
function zuiedit_data_coef_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_coef') ;

% information pour l'assistant
zuicr(hfig,action) ;

% variables d'entrée de l'éditeur de consignes zuieditcons
x           = evalin('base','param.gene.x') ;
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

% selon ation
switch lower(action)

case 'edit_eta'
	nom        = 'data.coef.eta' ;
	y          = evalin('base','data.coef.eta') ;
	texte_y    = 'eta' ;
	var_y      = 'data.coef.eta';
	texte_prop = 'te0' ;
	var_prop   = 'data.cons.asser.te0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_eta) ;

case 'import_eta'
	nom     = 'data.coef.eta' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_eta) ;
	
case 'dessin_eta'
	hout = zuiedit_profcmplx('data.coef.eta','data.gene.temps','param.gene.x')
	zuireset(h.dessin_eta) ;

case 'edit_ee'
	nom        = 'data.coef.ee' ;
	y          = evalin('base','data.coef.ee') ;
	texte_y    = 'ee' ;
	var_y      = 'data.coef.ee';
	texte_prop = 'te0' ;
	var_prop   = 'data.cons.asser.te0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ee) ;

case 'import_ee'
	nom     = 'data.coef.ee' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ee) ;

case 'dessin_ee'
	hout = zuiedit_profcmplx('data.coef.ee','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ee) ;

case 'edit_ei'
	nom        = 'data.coef.ei' ;
	y          = evalin('base','data.coef.ei') ;
	texte_y    = 'ei' ;
	var_y      = 'data.coef.ei';
	texte_prop = 'ti0' ;
	var_prop   = 'data.cons.asser.ti0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ei) ;

case 'import_ei'
	nom     = 'data.coef.ei' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ei) ;

case 'dessin_ei'
	hout = zuiedit_profcmplx('data.coef.ei','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ei) ;

case 'edit_en'
	nom        = 'data.coef.en' ;
	y          = evalin('base','data.coef.en') ;
	texte_y    = 'en' ;
	var_y      = 'data.coef.en';
	texte_prop = 'ne0' ;
	var_prop   = 'data.cons.asser.ne0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_en) ;

case 'import_en'
	nom     = 'data.coef.en' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_en) ;

case 'dessin_en'
	hout = zuiedit_profcmplx('data.coef.en','data.gene.temps','param.gene.x')
	zuireset(h.dessin_en) ;

case 'edit_ej'
	nom        = 'data.coef.ej' ;
	y          = evalin('base','data.coef.ej') ;
	texte_y    = 'ej' ;
	var_y      = 'data.coef.ej';
	texte_prop = 'q0' ;
	var_prop   = 'data.cons.asser.q0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ej) ;

case 'import_ej'
	nom     = 'data.coef.ej' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ej) ;

case 'dessin_ej'
	hout = zuiedit_profcmplx('data.coef.ej','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ej) ;

case 'edit_ve'
	nom        = 'data.coef.ve' ;
	y          = evalin('base','data.coef.ve') ;
	texte_y    = 've' ;
	var_y      = 'data.coef.ve';
	texte_prop = 'vloop' ;
	var_prop   = 'data.cons.asser.vloop' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ve) ;

case 'import_ve'
	nom     = 'data.coef.ve' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ve) ;

case 'dessin_ve'
	hout = zuiedit_profcmplx('data.coef.ve','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ve) ;

case 'edit_ep'
	nom        = 'data.coef.ep' ;
	y          = evalin('base','data.coef.ep') ;
	texte_y    = 'ep' ;
	var_y      = 'data.coef.ep';
	texte_prop = 'te0' ;
	var_prop   = 'data.cons.asser.te0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ep) ;

case 'import_ep'
	nom     = 'data.coef.ep' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ep) ;

case 'dessin_ep'
	hout = zuiedit_profcmplx('data.coef.ep','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ep) ;

case 'edit_ie'
	nom        = 'data.coef.ie' ;
	y          = evalin('base','data.coef.ie') ;
	texte_y    = 'ie' ;
	var_y      = 'data.coef.ie' ;
	texte_prop = 'te0' ;
	var_prop   = 'data.cons.asser.te0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ie) ;

case 'import_ie'
	nom     = 'data.coef.ie' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ie) ;

case 'dessin_ie'
	hout = zuiedit_profcmplx('data.coef.ie','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ie) ;

case 'edit_ii'
	nom        = 'data.coef.ii' ;
	y          = evalin('base','data.coef.ii') ;
	texte_y    = 'ii' ;
	var_y      = 'data.coef.ii';
	texte_prop = 'ti0' ;
	var_prop   = 'data.cons.asser.ti0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ii) ;

case 'import_ii'
	nom     = 'data.coef.ii' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ii) ;

case 'dessin_ii'
	hout = zuiedit_profcmplx('data.coef.ii','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ii) ;

case 'edit_in'
	nom        = 'data.coef.in' ;
	y          = evalin('base','data.coef.in') ;
	texte_y    = 'in' ;
	var_y      = 'data.coef.in';
	texte_prop = 'ne0' ;
	var_prop   = 'data.cons.asser.ne0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_in) ;

case 'import_in'
	nom     = 'data.coef.in' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_in) ;

case 'dessin_in'
	hout = zuiedit_profcmplx('data.coef.in','data.gene.temps','param.gene.x')
	zuireset(h.dessin_in) ;

case 'edit_ij'
	nom        = 'data.coef.ij' ;
	y          = evalin('base','data.coef.ij') ;
	texte_y    = 'ij' ;
	var_y      = 'data.coef.ij';
	texte_prop = 'q0' ;
	var_prop   = 'data.cons.asser.q0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ij) ;

case 'import_ij'
	nom     = 'data.coef.ij' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ij) ;

case 'dessin_ij'
	hout = zuiedit_profcmplx('data.coef.ij','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ij) ;

case 'edit_vi'
	nom        = 'data.coef.vi' ;
	y          = evalin('base','data.coef.vi') ;
	texte_y    = 'vi' ;
	var_y      = 'data.coef.vi';
	texte_prop = 'vloop' ;
	var_prop   = 'data.cons.asser.vloop' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_vi) ;

case 'import_vi'
	nom     = 'data.coef.vi' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_vi) ;

case 'dessin_vi'
	hout = zuiedit_profcmplx('data.coef.vi','data.gene.temps','param.gene.x')
	zuireset(h.dessin_vi) ;

case 'edit_ip'
	nom        = 'data.coef.ip' ;
	y          = evalin('base','data.coef.ip') ;
	texte_y    = 'ip' ;
	var_y      = 'data.coef.ip';
	texte_prop = 'ti0' ;
	var_prop   = 'data.cons.asser.ti0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ip) ;

case 'import_ip'
	nom     = 'data.coef.ip' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ip) ;

case 'dessin_ip'
	hout = zuiedit_profcmplx('data.coef.ip','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ip) ;

case 'edit_ne'
	nom        = 'data.coef.ne' ;
	y          = evalin('base','data.coef.ne') ;
	texte_y    = 'ne' ;
	var_y      = 'data.coef.ne';
	texte_prop = 'te0' ;
	var_prop   = 'data.cons.asser.te0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ne) ;

case 'import_ne'
	nom     = 'data.coef.ne' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ne) ;

case 'dessin_ne'
	hout = zuiedit_profcmplx('data.coef.ne','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ne) ;

case 'edit_ni'
	nom        = 'data.coef.ni' ;
	y          = evalin('base','data.coef.ni') ;
	texte_y    = 'ni' ;
	var_y      = 'data.coef.ni';
	texte_prop = 'ti0' ;
	var_prop   = 'data.cons.asser.ti0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ni) ;

case 'import_ni'
	nom     = 'data.coef.ni' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ni) ;

case 'dessin_ni'
	hout = zuiedit_profcmplx('data.coef.ni','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ni) ;

case 'edit_nn'
	nom        = 'data.coef.nn' ;
	y          = evalin('base','data.coef.nn') ;
	texte_y    = 'nn' ;
	var_y      = 'data.coef.nn';
	texte_prop = 'ne0' ;
	var_prop   = 'data.cons.asser.ne0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_nn) ;

case 'import_nn'
	nom     = 'data.coef.nn' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_nn) ;

case 'dessin_nn'
	hout = zuiedit_profcmplx('data.coef.nn','data.gene.temps','param.gene.x')
	zuireset(h.dessin_nn) ;

case 'edit_nj'
	nom        = 'data.coef.nj' ;
	y          = evalin('base','data.coef.nj') ;
	texte_y    = 'nj' ;
	var_y      = 'data.coef.nj';
	texte_prop = 'q0' ;
	var_prop   = 'data.cons.asser.q0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_nj) ;

case 'import_nj'
	nom     = 'data.coef.nj' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_nj) ;

case 'dessin_nj'
	hout = zuiedit_profcmplx('data.coef.nj','data.gene.temps','param.gene.x')
	zuireset(h.dessin_nj) ;

case 'edit_vn'
	nom        = 'data.coef.vn' ;
	y          = evalin('base','data.coef.vn') ;
	texte_y    = 'vn' ;
	var_y      = 'data.coef.vn';
	texte_prop = 'vloop' ;
	var_prop   = 'data.cons.vloop' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_vn) ;

case 'import_vn'
	nom     = 'data.coef.vn' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_vn) ;

case 'dessin_vn'
	hout = zuiedit_profcmplx('data.coef.vn','data.gene.temps','param.gene.x')
	zuireset(h.dessin_vn) ;

case 'edit_fev'
	nom        = 'data.coef.fev' ;
	y          = evalin('base','data.coef.fev') ;
	texte_y    = 'fev' ;
	var_y      = 'data.coef.fev';
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_fev) ;

case 'import_fev'
	nom     = 'data.coef.fev' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_fev) ;

case 'dessin_fev'
	hout = zuiedit_profcmplx('data.coef.fev','data.gene.temps','param.gene.x')
	zuireset(h.dessin_fev) ;

case 'edit_fefe'
	nom        = 'data.coef.fefe' ;
	y          = evalin('base','data.coef.fefe') ;
	texte_y    = 'fefe' ;
	var_y      = 'data.coef.fefe';
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_fefe) ;

case 'import_fefe'
	nom     = 'data.coef.fefe' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_fefe) ;

case 'dessin_fefe'
	hout = zuiedit_profcmplx('data.coef.fefe','data.gene.temps','param.gene.x')
	zuireset(h.dessin_fefe) ;

case 'edit_fiv'
	nom        = 'data.coef.fiv' ;
	y          = evalin('base','data.coef.fiv') ;
	texte_y    = 'fiv' ;
	var_y      = 'data.coef.fiv' ;
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_fiv) ;

case 'import_fiv'
	nom     = 'data.coef.fiv' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_fiv) ;

case 'dessin_fiv'
	hout = zuiedit_profcmplx('data.coef.fiv','data.gene.temps','param.gene.x')
	zuireset(h.dessin_fiv) ;

case 'edit_fifi'
	nom        = 'data.coef.fifi' ;
	y          = evalin('base','data.coef.fifi') ;
	texte_y    = 'fifi' ;
	var_y      = 'data.coef.fifi';
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_fifi) ;

case 'import_fifi'
	nom     = 'data.coef.fifi' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_fifi) ;

case 'dessin_fifi'
	hout = zuiedit_profcmplx('data.coef.fifi','data.gene.temps','param.gene.x')
	zuireset(h.dessin_fifi) ;

case 'edit_rotv'
	nom        = 'data.coef.rotv' ;
	y          = evalin('base','data.coef.rotv') ;
	texte_y    = 'rotv' ;
	var_y      = 'data.coef.rotv';
	texte_prop = 'vloop' ;
	var_prop   = 'data.cons.vloop' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_rotv) ;

case 'import_rotv'
	nom     = 'data.coef.rotv' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_rotv) ;

case 'dessin_rotv'
	hout = zuiedit_profcmplx('data.coef.rotv','data.gene.temps','param.gene.x')
	zuireset(h.dessin_rotv) ;

case 'edit_rot'
	nom        = 'data.coef.rot' ;
	y          = evalin('base','data.coef.rot') ;
	texte_y    = 'rot' ;
	var_y      = 'data.coef.rot';
	texte_prop = 'ti0' ;
	var_prop   = 'data.cons.asser.ti0' ;
	zuieditcons(nom,info.data.coef.eta,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_rot) ;

case 'import_rot'
	nom     = 'data.coef.rot' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_rot) ;

case 'dessin_rot'
	hout = zuiedit_profcmplx('data.coef.rot','data.gene.temps','param.gene.x')
	zuireset(h.dessin_rot) ;

case {'btn_quit','close'}
	zuicloseone(hfig);	
	
case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
otherwise
	warning('ation non prise en compte')
	   
end

