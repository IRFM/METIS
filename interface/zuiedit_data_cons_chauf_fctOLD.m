%
% ZUIDEDIT_DATA_CONS_CHAUF_FCT
%      	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de consignes
%--------------------------------------------------------------
% fonction Matlab 5 :
%
% fichier zuiedit_data_cons_chauf_fct.m  
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
% version 1.6, du 31/07/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
function zuiedit_data_cons_chauf_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_cons_chauf') ;

% information pour l'assistant
zuicr(hfig,action) ;

% variables d'entrée de l'éditeur de consignes zuieditcons
x          = evalin('base','data.gene.temps') ;
texte_x    = 'x' ;
var_x      = 'void';
liste_ref  = {} ;
var_ref    = {} ;
texte_prop = '' ;
var_prop   = '' ;

% selon ation
switch lower(action)

case 'edit_module_fci1'
	nom     = 'data.cons.fci(:,1)' ;
	y       = evalin('base','data.cons.fci(:,1)') ;
	texte_y = 'fci(:,1)' ;
	var_y   = 'data.cons.fci';
	canal   = 1 ;
	code_retour = 'abs' ;
	hout = zuieditcons(nom,info.data.cons.fci,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_module_fci1) ;

case 'edit_phase_fci1'
	nom     = 'data.cons.fci(:,1)' ;
	y       = evalin('base','data.cons.fci(:,1)') ;
	texte_y = 'fci(:,1)' ;
	var_y   = 'data.cons.fci';
	canal   = 1 ;
	code_retour = 'degres' ;
	hout = zuieditcons(nom,info.data.cons.fci,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_phase_fci1) ;

case 'import_fci1'
	hout = zuiedit_import_mode('data.cons.fci(:,1)','consigne')
	zuireset(h.import_fci1) ;

case 'edit_module_fci2'
	nom     = 'data.cons.fci(:,2)' ;
	y       = evalin('base','data.cons.fci(:,2)') ;
	texte_y = 'fci(:,2)' ;
	var_y   = 'data.cons.fci';
	canal   = 2 ;
	code_retour = 'abs' ;
	hout = zuieditcons(nom,info.data.cons.fci,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_fci2) ;

case 'edit_phase_fci2'
	nom     = 'data.cons.fci(:,2)' ;
	y       = evalin('base','data.cons.fci(:,2)') ;
	texte_y = 'fci(:,2)' ;
	var_y   = 'data.cons.fci';
	canal   = 2 ;
	code_retour = 'degres' ;
	hout = zuieditcons(nom,info.data.cons.fci,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_phase_fci2) ;

case 'prec_fci2'
	y = evalin('base','data.cons.fci(:,1)') ;
	zassignin('base','data.cons.fci(:,2)',y) ;
	zuireset(h.prec_fci2)

case 'import_fci2'
	hout = zuiedit_import_mode('data.cons.fci(:,2)','consigne')
	zuireset(h.import_fci2) ;

case 'edit_module_fci3'
	nom     = 'data.cons.fci(:,3)' ;
	y       = evalin('base','data.cons.fci(:,3)') ;
	texte_y = 'fci(:,3)' ;
	var_y   = 'data.cons.fci';
	canal   = 3 ;
	code_retour = 'abs' ;
	hout = zuieditcons(nom,info.data.cons.fci,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_fci3) ;

case 'edit_phase_fci3'
	nom     = 'data.cons.fci(:,3)' ;
	y       = evalin('base','data.cons.fci(:,3)') ;
	texte_y = 'fci(:,3)' ;
	var_y   = 'data.cons.fci';
	canal   = 3 ;
	code_retour = 'degres' ;
	hout = zuieditcons(nom,info.data.cons.fci,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_phase_fci3) ;

case 'prec_fci3'
	y = evalin('base','data.cons.fci(:,2)') ;
	zassignin('base','data.cons.fci(:,3)',y) ;
	zuireset(h.prec_fci3) ;

case 'import_fci3'
	hout = zuiedit_import_mode('data.cons.fci(:,3)','consigne')
	zuireset(h.import_fci3) ;

case 'edit_module_fce1'
	nom     = 'data.cons.fce(:,1)' ;
	y       = evalin('base','data.cons.fce(:,1)') ;
	texte_y = 'fce(:,1)' ;
	var_y   = 'data.cons.fce';
	canal   = 1 ;
	code_retour = 'abs' ;
	hout = zuieditcons(nom,info.data.cons.fce,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_fce1) ;

case 'edit_phase_fce1'
	nom     = 'data.cons.fce(:,1)' ;
	y       = evalin('base','data.cons.fce(:,1)') ;
	texte_y = 'fce(:,1)' ;
	var_y   = 'data.cons.fce';
	canal   = 1 ;
	code_retour = 'degres' ;
	hout = zuieditcons(nom,info.data.cons.fce,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_phase_fce1) ;

case 'import_fce1'
	hout = zuiedit_import_mode('data.cons.fce(:,1)','consigne')
	zuireset(h.import_fce1) ;

case 'edit_module_hyb1'
	nom     = 'data.cons.hyb(:,1)' ;
	y       = evalin('base','data.cons.hyb(:,1)') ;
	texte_y = 'hyb(:,1)' ;
	var_y   = 'data.cons.hyb';
	canal   = 1 ;
	code_retour = 'abs' ;
	hout = zuieditcons(nom,info.data.cons.hyb,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_hyb1) ;

case 'edit_phase_hyb1'
	nom     = 'data.cons.hyb(:,1)' ;
	y       = evalin('base','data.cons.hyb(:,1)') ;
	texte_y = 'hyb(:,1)' ;
	var_y   = 'data.cons.hyb';
	canal   = 1 ;
	code_retour = 'angle' ;
	hout = zuieditcons(nom,info.data.cons.hyb,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_phase_hyb1) ;

case 'import_hyb1'
	hout = zuiedit_import_mode('data.cons.hyb(:,1)','consigne')
	zuireset(h.import_hyb1) ;

case 'edit_module_hyb2'
	nom     = 'data.cons.hyb(:,2)' ;
	y       = evalin('base','data.cons.hyb(:,2)') ;
	texte_y = 'hyb(:,2)' ;
	var_y   = 'data.cons.hyb';
	canal   = 2 ;
	code_retour = 'abs' ;
	hout = zuieditcons(nom,info.data.cons.hyb,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_hyb2) ;

case 'edit_phase_hyb2'
	nom     = 'data.cons.hyb(:,2)' ;
	y       = evalin('base','data.cons.hyb(:,2)') ;
	texte_y = 'hyb(:,2)' ;
	var_y   = 'data.cons.hyb';
	canal   = 2 ;
	code_retour = 'angle' ;
	hout = zuieditcons(nom,info.data.cons.hyb,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_phase_hyb2) ;

case 'prec_hyb2'
	y = evalin('base','data.cons.hyb(:,1)') ;
	zassignin('base','data.cons.hyb(:,2)',y) ;
	zuireset(h.prec_hyb2) ;

case 'import_hyb2'
	hout = zuiedit_import_mode('data.cons.hyb(:,2)','consigne')
	zuireset(h.import_hyb2) ;

% %%%%%%
case {'edit_module_idn'}
	action
	nom     = 'data.cons.idn(:,1)' ;
	y       = evalin('base','data.cons.idn(:,1)') ;
	texte_y = 'idn(:,1)' ;1
	var_y   = 'data.cons.idn';
	canal   = 1 ;
	code_retour = 'abs' ;
	hout = zuieditcons(nom,info.data.cons.idn,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_module_idn1) ;

case 'edit_phase_idn1'
	nom     = 'data.cons.idn(:,1)' ;
	y       = evalin('base','data.cons.idn(:,1)') ;
	texte_y = 'idn(:,1)' ;
	var_y   = 'data.cons.idn';
	canal   = 2 ;
	code_retour = 'angle' ;
	hout = zuieditcons(nom,info.data.cons.idn,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(h.edit_phase_idn1) ;

case 'import_idn1'
	hout = zuiedit_import_mode('data.cons.idn(:,1)','consigne')
	zuireset(h.import_idn1) ;

case {'btn_quit','close'}
	zuicloseone(hfig);	
	
case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
otherwise
	warning('ation non prise en compte')
	   
end

