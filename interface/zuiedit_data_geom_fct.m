% ZUIDEDIT_DATA_GEOM_FCT gestion des callbacks du formulaire d'edition des consignes de g�m�rie
%---------------------------------------------
% fichier zuiedit_data_geom_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	d'edition des consignes de geomerie
%
% syntaxe :
%	zuiedit_data_geom_fct(action)

% entrees :
%  action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 17/092001.
% 
% liste des modifications : 
%  * 27/09/2001 -> modification du libelle en x : texte_x = 'temps'
%
%--------------------------------------------------------------
function zuiedit_data_geom_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_geom') ;

% information pour l'assistant
zuicr(hfig,action) ;

% variables d'entree de l'editeur de consignes zuieditcons
x           = evalin('base','data.gene.temps') ;
texte_x     = 'time' ;
var_x       = 'void';
code_retour = 'abs' ;
liste_ref   = {} ;
var_ref     = {} ;
texte_prop  = '' ;
var_prop    = '' ;

% selon ation
switch lower(action)

case 'edit_r0'
	nom     = 'data.geo.r0' ;
	y       = evalin('base','data.geo.r0') ;
	texte_y = 'r0' ;
	var_y   = 'data.geo.r0';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.geo.r0,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.edit_r0) ;

case 'import_r0'
	hout = zuiedit_import_mode('data.geo.r0','consigne') ;
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.import_r0) ;

case 'edit_z0'
	nom     = 'data.geo.z0' ;
	y       = evalin('base','data.geo.z0') ;
	texte_y = 'z0' ;
	var_y   = 'data.geo.z0';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.geo.z0,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.edit_z0) ;

case 'import_z0'
	hout = zuiedit_import_mode('data.geo.z0','consigne') ;
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.import_z0) ;

case 'edit_a'
	nom     = 'data.geo.a' ;
	y       = evalin('base','data.geo.a') ;
	texte_y = 'a' ;
	var_y   = 'data.geo.a';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.geo.a,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.edit_a) ;

case 'import_a'
	hout = zuiedit_import_mode('data.geo.a','consigne') ;
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.import_a) ;

case 'edit_e1'
	nom     = 'data.geo.e1' ;
	y       = evalin('base','data.geo.e1') ;
	texte_y = 'e1' ;
	var_y   = 'data.geo.e1';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.geo.e1,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.edit_e1) ;

case 'import_e1'
	hout = zuiedit_import_mode('data.geo.e1','consigne') ;
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.import_e1) ;

case 'edit_trh1'
	nom     = 'data.geo.trh1' ;
	y       = evalin('base','data.geo.trh1') ;
	texte_y = 'trh1' ;
	var_y   = 'data.geo.trh1';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.geo.trh1,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.edit_trh1) ;

case 'import_trh1'
	hout = zuiedit_import_mode('data.geo.trh1','consigne') ;
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.import_trh1) ;

case 'edit_trb1'
	nom     = 'data.geo.trb1' ;
	y       = evalin('base','data.geo.trb1') ;
	texte_y = 'trb1' ;
	var_y   = 'data.geo.trb1';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.geo.trb1,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.edit_trb1) ;

case 'import_trb1'
	hout = zuiedit_import_mode('data.geo.trb1','consigne') ;
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.import_trb1) ;

case 'edit_ind1'
	nom     = 'data.geo.ind1' ;
	y       = evalin('base','data.geo.ind1') ;
	texte_y = 'ind1' ;
	var_y   = 'data.geo.ind1';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.geo.ind1,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.edit_ind1) ;

case 'import_ind1'
	hout = zuiedit_import_mode('data.geo.ind1','consigne') ;
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.import_ind1) ;

case 'edit_b0'
	nom     = 'data.geo.b0' ;
	y       = evalin('base','data.geo.b0') ;
	texte_y = 'b0' ;
	var_y   = 'data.geo.b0';
	canal   = 1 ;
	hout = zuieditcons(nom,info.data.geo.b0,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	setappdata(hfig,'liste_tag',liste_tag) ;
	zuireset(h.edit_b0) ;

case 'import_b0'
	hout = zuiedit_import_mode('data.geo.b0','consigne') ;
	liste_tag = getappdata(hfig,'liste_tag') ;
	tag = get(hout,'tag') ;
	liste_tag{end+1} = tag ;
	zuireset(h.import_b0) ;

case 'import_rz'
	file     = evalin('base','param.edit.currentfile') ;
	liste_tag = getappdata(hfig,'liste_tag') ;
	fh        =  getappdata(0,'formulaire_handle') ;
	for i=1:length(liste_tag)
		hl = getfield(fh,liste_tag{i}) ;
 		if ishandle(hl)
 			herror = warndlg('separatrix parameter edition -> importation not allowed','Problem','modal') ;
 			zuireset(h.import_rz) ;
 			return
  		end
 	end	
	cr = zimport_sepa(hfig) ;
	zuireset(h.import_rz) ;
% RAZ
case {'init','raz'}
	zuiformvisible(hfig) ;
	zuiformreset(hfig) ;
	zuireset(h.raz) ;
	zuiuploadform(hfig) ;
	
% Annulation
case {'annulation','close'}
	zuireset(h.annulation) ;
	zuicloseone(hfig);	
	
% Validation
case 'validation'
	zuiformcache(hfig) ;
	zuireset(h.validation) ;
	zuidownloadform(hfig);
	var = evalin('base','data.geo.mode') ;
	val = zuidata(h.pop_mode) ;
	zassignin('base','data.geo.mode',ones(size(var))*val) ;
	zuisavenonok;
	
case 'pop_mode'
	%rien
otherwise
	warning('action not taken into account')
	   
end

