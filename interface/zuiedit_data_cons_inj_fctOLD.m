%
% ZUIDEDIT_DATA_CONS_INJ_FCT
%      	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de consignes
%--------------------------------------------------------------
% fonction Matlab 5 :
%
% fichier zuiedit_data_cons_inj_fct.m  
%
% syntaxe :
% 
% entrees :
% 
%  action       =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 2/08/2001.
% 
% liste des modifications : 
%	* 13/09/2001 -> Pour les glacons, appel de l'editeur zuiedit_mode
%
%--------------------------------------------------------------
function zuiedit_data_cons_inj_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_cons_inj') ;

% information pour l'assistant
zuicr(hfig,action) ;

% variables d'entrée de l'éditeur de consignes zuieditcons
x           = evalin('base','data.gene.temps') ;
texte_x     = 'x' ;
var_x       = 'void';
code_retour = '' ;
liste_ref   = {} ;
var_ref     = {} ;
texte_prop  = '' ;
var_prop    = '' ;

% selon ation
switch lower(action)

% data.cons.c
% -----------
case 'edit_module_c1'
	nom     = 'data.cons.c(:,1)' ;
	y       = evalin('base','data.cons.c(:,1)') ;
	texte_y = 'c(:,1)' ;
	var_y   = 'data.cons.c';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.c,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_module_c1) ;

case 'import_c1'
	hout = zuiedit_import_mode('data.cons.c(:,1)','consigne') ;
	zuireset(h.import_c1) ;

case 'edit_module_c2'
	nom     = 'data.cons.c(:,2)' ;
	y       = evalin('base','data.cons.c(:,2)') ;
	texte_y = 'c(:,2)' ;
	var_y   = 'data.cons.c';
	canal   = 2 ;
	hout = zuieditcons(nom,info.data.cons.c,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_c2) ;

case 'prec_c2'
	rep =questdlg('On recopie les valeurs du canal précédent ?', ...
	              'Confirmation', ...
	              'Oui','Non','Oui');
	switch rep
	case 'Oui'
		y = evalin('base','data.cons.c(:,1)') ;
		zassignin('base','data.cons.c(:,2)',y) ;
	case 'Non'
	end
	zuireset(h.prec_c2)

case 'import_c2'
	hout = zuiedit_import_mode('data.cons.c(:,2)','consigne') ;
	zuireset(h.import_c2) ;

case 'edit_module_c3'
	nom     = 'data.cons.c(:,3)' ;
	y       = evalin('base','data.cons.c(:,3)') ;
	texte_y = 'c(:,3)' ;
	var_y   = 'data.cons.c';
	canal   = 3 ;
	hout = zuieditcons(nom,info.data.cons.c,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_c3) ;

case 'prec_c3'
	rep =questdlg('On recopie les valeurs du canal précédent ?', ...
	              'Confirmation', ...
	              'Oui','Non','Oui');
	switch rep
	case 'Oui'
		y = evalin('base','data.cons.c(:,2)') ;
		zassignin('base','data.cons.c(:,3)',y) ;
	case 'Non'
	end
	zuireset(h.prec_c3)

case 'import_c3'
	hout = zuiedit_import_mode('data.cons.c(:,3)','consigne') ;
	zuireset(h.import_c3) ;

case 'edit_module_c4'
	nom     = 'data.cons.c(:,4)' ;
	y       = evalin('base','data.cons.c(:,4)') ;
	texte_y = 'c(:,4)' ;
	var_y   = 'data.cons.c';
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
		y = evalin('base','data.cons.c(:,3)') ;
		zassignin('base','data.cons.c(:,4)',y) ;
	case 'Non'
	end
	zuireset(h.prec_c4)

case 'import_c4'
	hout = zuiedit_import_mode('data.cons.c(:,4)','consigne') ;
	zuireset(h.import_c4) ;

case 'edit_module_c5'
	nom     = 'data.cons.c(:,5)' ;
	y       = evalin('base','data.cons.c(:,5)') ;
	texte_y = 'c(:,5)' ;
	var_y   = 'data.cons.c';
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
		y = evalin('base','data.cons.c(:,4)') ;
		zassignin('base','data.cons.c(:,5)',y) ;
	case 'Non'
	end
	zuireset(h.prec_c5)

case 'import_c5'
	hout = zuiedit_import_mode('data.cons.c(:,5)','consigne') ;
	zuireset(h.import_c5) ;

% data.cons.pomp
% --------------
case 'edit_module_pomp'
	nom     = 'data.cons.pomp' ;
	y       = evalin('base','data.cons.pomp') ;
	texte_y = 'pomp' ;
	var_y   = 'data.cons.pomp';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.pomp,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_pomp) ;

case 'import_pomp'
	hout = zuiedit_import_mode('data.cons.pomp','consigne') ;
	zuireset(h.import_pomp) ;

% data.cons.glacon
% ----------------
case 'edit_module_glacon1'
	zuiedit_mode('data.cons.glacon(:,1)',{0,1},{'attente','lancement'}) ;
%	nom     = 'data.cons.glacon(:,1)' ;
%	y       = evalin('base','data.cons.glacon(:,1)') ;
%	texte_y = 'glacon(:,1)' ;
%	var_y   = 'data.cons.glacon';
%	canal   = 1 ;
%	hout = zuieditcons(nom,info.data.cons.glacon,x,y,texte_x,texte_y,var_x,var_y,canal, ...
%	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_glacon1) ;

case 'import_glacon1'
	hout = zuiedit_import_mode('data.cons.glacon(:,1)','consigne') ;
	zuireset(h.import_glacon1) ;

case 'edit_module_glacon2'
	zuiedit_mode('data.cons.glacon(:,2)',{0,1},{'attente','lancement'}) ;
	zuireset(h.edit_module_glacon2) ;

case 'prec_glacon2'
	rep =questdlg('On recopie les valeurs du canal précédent ?', ...
	              'Confirmation', ...
	              'Oui','Non','Oui');
	switch rep
	case 'Oui'
		y = evalin('base','data.cons.glacon(:,1)') ;
		zassignin('base','data.cons.glacon(:,2)',y) ;
	case 'Non'
	end
	zuireset(h.prec_glacon2)

case 'import_glacon2'
	hout = zuiedit_import_mode('data.cons.glacon(:,2)','consigne') ;
	zuireset(h.import_glacon2) ;

case 'edit_module_glacon3'
	zuiedit_mode('data.cons.glacon(:,3)',{0,1},{'attente','lancement'}) ;
	zuireset(h.edit_module_glacon3) ;

case 'prec_glacon3'
	rep =questdlg('On recopie les valeurs du canal précédent ?', ...
	              'Confirmation', ...
	              'Oui','Non','Oui');
	switch rep
	case 'Oui'
		y = evalin('base','data.cons.glacon(:,2)') ;
		zassignin('base','data.cons.glacon(:,3)',y) ;
	case 'Non'
	end
	zuireset(h.prec_glacon3)

case 'import_glacon3'
	hout = zuiedit_import_mode('data.cons.glacon(:,31)','consigne') ;
	zuireset(h.import_glacon3) ;

% data.cons.zeffm
% ---------------
case 'edit_module_zeffm'
	nom     = 'data.cons.zeffm' ;
	y       = evalin('base','data.cons.zeffm') ;
	texte_y = 'zeffm' ;
	var_y   = 'data.cons.zeffm';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.zeffm,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_zeffm) ;

case 'import_zeffm'
	hout = zuiedit_import_mode('data.cons.zeffm','consigne') ;
	zuireset(h.import_zeffm) ;

% data.cons.nhnd
% ---------------
case 'edit_module_nhnd'
	nom     = 'data.cons.nhnd' ;
	y       = evalin('base','data.cons.nhnd') ;
	texte_y = 'nhnd' ;
	var_y   = 'data.cons.nhnd';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.cons.nhnd,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_nhnd) ;

case 'import_nhnd'
	hout = zuiedit_import_mode('data.cons.nhnd','consigne') ;
	zuireset(h.import_nhnd) ;

% autres
% ------
case {'btn_quit','close'}
	zuicloseone(hfig);	
	
case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
otherwise
	warning('ation non prise en compte')
	   
end

