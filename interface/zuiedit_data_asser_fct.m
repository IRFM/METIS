% ZUIDEDIT_DATA_ASSER_FCT gestion callbacks du formulaire d'edition de consignes des asservissements
%--------------------------------------------------------------
% fichier zuiedit_data_asser_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de consignes des asservissements
%
% syntaxe :
% 
% entrees :
%  action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 02/08/2001.
% 
% liste des modifications : 
%
%  * 27/09/2001 -> modification du libelle en x : texte_x = 'temps'
%--------------------------------------------------------------
function zuiedit_data_asser_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_asser') ;

% information pour l'assistant
zuicr(hfig,action) ;

% variables d'entrée de l'éditeur de consignes zuieditcons
x           = evalin('base','data.gene.temps') ;
texte_x     = 'temps' ;
var_x       = 'void';
code_retour = 'abs' ;
liste_ref   = {} ;
var_ref     = {} ;
texte_prop  = '' ;
var_prop    = '' ;

% selon ation
switch lower(action)

case 'edit_nl0'
	nom     = 'data.cons.asser.nl0' ;
	y       = evalin('base','data.cons.asser.nl0') ;
	texte_y = 'nl0' ;
	var_y   = 'data.cons.asser.nl0';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.asser.nl0,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_nl0) ;

case 'import_nl0'
	hout = zuiedit_import_mode('data.cons.asser.nl0','consigne') ;
	zuireset(h.import_nl0) ;

case 'edit_ne0'
	nom     = 'data.cons.asser.ne0' ;
	y       = evalin('base','data.cons.asser.ne0') ;
	texte_y = 'ne0' ;
	var_y   = 'data.cons.asser.ne0';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.asser.ne0,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_ne0) ;

case 'import_ne0'
	hout = zuiedit_import_mode('data.cons.asser.ne0','consigne') ;
	zuireset(h.import_ne0) ;

case 'edit_ne1'
	nom     = 'data.cons.ne1' ;
	y       = evalin('base','data.cons.asser.ne1') ;
	texte_y = 'ne1' ;
	var_y   = 'data.cons.asser.ne1';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.asser.ne1,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_ne1) ;

case 'import_ne1'
	hout = zuiedit_import_mode('data.cons.asser.ne1','consigne') ;
	zuireset(h.import_ne1) ;

case 'edit_nemoy'
	nom     = 'data.cons.asser.nemoy' ;
	y       = evalin('base','data.cons.asser.nemoy') ;
	texte_y = 'nemoy' ;
	var_y   = 'data.cons.asser.nemoy';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.asser.nemoy,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_nemoy) ;

case 'import_nemoy'
	hout = zuiedit_import_mode('data.cons.asser.nemoy','consigne') ;
	zuireset(h.import_nemoy) ;

case 'edit_te0'
	nom     = 'data.cons.asser.te0' ;
	y       = evalin('base','data.cons.asser.te0') ;
	texte_y = 'te0' ;
	var_y   = 'data.cons.asser.te0';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.asser.te0,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_te0) ;

case 'import_te0'
	hout = zuiedit_import_mode('data.cons.asser.te0','consigne') ;
	zuireset(h.import_te0) ;

case 'edit_te1'
	nom     = 'data.cons.asser.te1' ;
	y       = evalin('base','data.cons.asser.te1') ;
	texte_y = 'te1' ;
	var_y   = 'data.cons.asser.te1';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.asser.te1,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_te1) ;

case 'import_te1'
	hout = zuiedit_import_mode('data.cons.asser.te1','consigne') ;
	zuireset(h.import_te1) ;

case 'edit_ti0'
	nom     = 'data.cons.asser.ti0' ;
	y       = evalin('base','data.cons.asser.ti0') ;
	texte_y = 'ti0' ;
	var_y   = 'data.cons.asser.ti0';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.asser.ti0,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_ti0) ;

case 'import_ti0'
	hout = zuiedit_import_mode('data.cons.asser.ti0','consigne') ;
	zuireset(h.import_ti0) ;

case 'edit_ti1'
	nom     = 'data.cons.asser.ti1' ;
	y       = evalin('base','data.cons.asser.ti1') ;
	texte_y = 'ti1' ;
	var_y   = 'data.cons.asser.ti1';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.asser.ti1,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_ti1) ;

case 'import_ti1'
	hout = zuiedit_import_mode('data.cons.asser.ti1','consigne') ;
	zuireset(h.import_ti1) ;

case 'edit_li'
	nom     = 'data.cons.asser.li' ;
	y       = evalin('base','data.cons.asser.li') ;
	texte_y = 'li' ;
	var_y   = 'data.cons.asser.li';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.asser.li,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_li) ;

case 'import_li'
	hout = zuiedit_import_mode('data.cons.asser.li','consigne') ;
	zuireset(h.import_li) ;

case 'edit_q0'
	nom     = 'data.cons.asser.q0' ;
	y       = evalin('base','data.cons.asser.q0') ;
	texte_y = 'q0' ;
	var_y   = 'data.cons.asser.q0';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.asser.q0,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_q0) ;

case 'import_q0'
	hout = zuiedit_import_mode('data.cons.asser.q0','consigne') ;
	zuireset(h.import_q0) ;

case 'edit_module_c1'
	nom     = 'data.cons.asser.c(:,1)' ;
	y       = evalin('base','data.cons.asser.c(:,1)') ;
	texte_y = 'c(:,1)' ;
	var_y   = 'data.cons.asser.c';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.asser.c,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_module_c1) ;

case 'import_c1'
	hout = zuiedit_import_mode('data.cons.asser.c(:,1)','consigne') ;
	zuireset(h.import_c1) ;

case 'edit_module_c2'
	nom     = 'data.cons.asser.c(:,2)' ;
	y       = evalin('base','data.cons.asser.c(:,2)') ;
	texte_y = 'c(:,2)' ;
	var_y   = 'data.cons.asser.c';
	canal   = 2 ;
	hout = zuieditcons(nom,info.data.cons.asser.c,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_c2) ;

case 'prec_c2'
	rep =questdlg('On recopie les valeurs du canal précédent ?', ...
	              'Confirmation', ...
	              'Oui','Non','Oui');
	switch rep
	case 'Oui'
		y = evalin('base','data.cons.asser.c(:,1)') ;
		zassignin('base','data.cons.asser.c(:,2)',y) ;
	case 'Non'
	end
	zuireset(h.prec_c2)

case 'import_c2'
	hout = zuiedit_import_mode('data.cons.asser.c(:,2)','consigne') ;
	zuireset(h.import_c2) ;

case 'edit_module_c3'
	nom     = 'data.cons.asser.c(:,3)' ;
	y       = evalin('base','data.cons.asser.c(:,3)') ;
	texte_y = 'c(:,3)' ;
	var_y   = 'data.cons.asser.c';
	canal   = 3 ;
	hout = zuieditcons(nom,info.data.cons.asser.c,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_c3) ;

case 'prec_c3'
	rep =questdlg('On recopie les valeurs du canal précédent ?', ...
	              'Confirmation', ...
	              'Oui','Non','Oui');
	switch rep
	case 'Oui'
		y = evalin('base','data.cons.asser.c(:,2)') ;
		zassignin('base','data.cons.asser.c(:,3)',y) ;
	case 'Non'
	end
	zuireset(h.prec_c3)

case 'import_c3'
	hout = zuiedit_import_mode('data.cons.asser.c(:,3)','consigne') ;
	zuireset(h.import_c3) ;

case 'edit_module_c4'
	nom     = 'data.cons.asser.c(:,4)' ;
	y       = evalin('base','data.cons.asser.c(:,4)') ;
	texte_y = 'c(:,4)' ;
	var_y   = 'data.cons.asser.c';
	canal   = 4 ;
	hout = zuieditcons(nom,info.data.cons.c,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_c4) ;

case 'prec_c4'
	rep =questdlg('On recopie les valeurs du canal précédent ?', ...
	              'Confirmation', ...
	              'Oui','Non','Oui');
	switch rep
	case 'Oui'
		y = evalin('base','data.cons.asser.c(:,3)') ;
		zassignin('base','data.cons.asser.c(:,4)',y) ;
	case 'Non'
	end
	zuireset(h.prec_c4)

case 'import_c4'
	hout = zuiedit_import_mode('data.cons.asser.c(:,4)','consigne') ;
	zuireset(h.import_c4) ;

case 'edit_module_c5'
	nom     = 'data.cons.asser.c(:,5)' ;
	y       = evalin('base','data.cons.asser.c(:,5)') ;
	texte_y = 'c(:,5)' ;
	var_y   = 'data.cons.asser.c';
	canal   = 5 ;
	hout = zuieditcons(nom,info.data.cons.c,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_c5) ;

case 'prec_c5'
	rep =questdlg('On recopie les valeurs du canal précédent ?', ...
	              'Confirmation', ...
	              'Oui','Non','Oui');
	switch rep
	case 'Oui'
		y = evalin('base','data.cons.asser.c(:,4)') ;
		zassignin('base','data.cons.asser.c(:,5)',y) ;
	case 'Non'
	end
	zuireset(h.prec_c5)

case 'import_c5'
	hout = zuiedit_import_mode('data.cons.asser.c(:,5)','consigne') ;
	zuireset(h.import_c5) ;

case {'btn_quit','close'}
	zuicloseone(hfig);	
	
case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
otherwise
	warning('ation non prise en compte')
	   
end

