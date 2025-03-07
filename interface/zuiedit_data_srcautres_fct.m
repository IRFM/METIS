% ZUIDEDIT_DATA_SRCAUTRES_FCT  gestion  callbacks du formulaired'edition de profils de source-autres
%--------------------------------------------------------------------
% fichier zuiedit_data_srcautres_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de profils de source-autres
%
% syntaxe :
% 
% entrees :
%  action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 2.1, du 27/08/2003.
% 
% liste des modifications : 
%
% * 27/08/2003 -> ajout de kini pour edit  
%
%--------------------------------------------------------------
function zuiedit_data_srcautres_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_srcautres') ;

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

% selon ation
switch lower(action)

case 'edit_prad'
	nom        = 'data.source.prad' ;
	y          = evalin('base','data.source.prad') ;
	texte_y    = 'prad' ;
	var_y      = 'data.source.prad';
	texte_prop = 'zeffm' ;
	var_prop   = 'data.cons.zeffm' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_prad) ;

case 'import_prad'
	nom     = 'data.source.prad' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_prad) ;
	
case 'dessin_prad'
	hout = zuiedit_profcmplx('data.source.prad','data.gene.temps','param.gene.x')
	zuireset(h.dessin_prad) ;


case 'edit_brem'
	nom        = 'data.source.brem' ;
	y          = evalin('base','data.source.brem') ;
	texte_y    = 'brem' ;
	var_y      = 'data.source.brem';
	texte_prop = 'asser.te0' ;
	var_prop   = 'data.cons.asser.te0' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_brem) ;

case 'import_brem'
	nom     = 'data.source.brem' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_brem) ;

case 'dessin_brem'
	hout = zuiedit_profcmplx('data.source.brem','data.gene.temps','param.gene.x')
	zuireset(h.dessin_brem) ;


case 'edit_cyclo'
	nom        = 'data.source.cyclo' ;
	y          = evalin('base','data.source.cyclo') ;
	texte_y    = 'cyclo' ;
	var_y      = 'data.source.cyclo';
	texte_prop = 'asser.te0' ;
	var_prop   = 'data.cons.asser.te0' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_cyclo) ;

case 'import_cyclo'
	nom     = 'data.source.cyclo' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_cyclo) ;

case 'dessin_cyclo'
	hout = zuiedit_profcmplx('data.source.cyclo','data.gene.temps','param.gene.x')
	zuireset(h.dessin_cyclo) ;


case 'edit_ohm'
	nom        = 'data.source.ohm' ;
	y          = evalin('base','data.source.ohm') ;
	texte_y    = 'ohm' ;
	var_y      = 'data.source.ohm';
	texte_prop = 'vloop' ;
	var_prop   = 'data.cons.vloop' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ohm) ;

case 'import_ohm'
	nom     = 'data.source.ohm' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ohm) ;

case 'dessin_ohm'
	hout = zuiedit_profcmplx('data.source.ohm','data.gene.temps','param.gene.x')
	zuireset(h.dessin_ohm) ;

case 'edit_qei'
	nom        = 'data.source.qei' ;
	y          = evalin('base','data.source.qei') ;
	texte_y    = 'qei' ;
	var_y      = 'data.source.qei';
	texte_prop = 'asser.ti0' ;
	var_prop   = 'data.cons.asser.ti0' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_qei) ;

case 'import_qei'
	nom     = 'data.source.qei' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_qei) ;

case 'dessin_qei'
	hout = zuiedit_profcmplx('data.source.qei','data.gene.temps','param.gene.x')
	zuireset(h.dessin_qei) ;

case 'edit_jboot'
	nom        = 'data.source.jboot' ;
	y          = evalin('base','data.source.jboot') ;
	texte_y    = 'jboot' ;
	var_y      = 'data.source.jboot';
	texte_prop = 'asser.ti0' ;
	var_prop   = 'data.cons.asser.ti0' ;
	zuieditcons(nom,info.data.coef.eta,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_jboot) ;

case 'import_jboot'
	nom     = 'data.source.jboot' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_jboot) ;

case 'dessin_jboot'
	hout = zuiedit_profcmplx('data.source.jboot','data.gene.temps','param.gene.x')
	zuireset(h.dessin_jboot) ;

case {'btn_quit','close'}
	zuicloseone(hfig);	
	
case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
otherwise
	warning('ation non prise en compte')
	   
end

