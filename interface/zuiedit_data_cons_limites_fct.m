% ZUIDEDIT_DATA_CONS_LIMITES_FCT  gestion des callbacks du formulaire de d'edition de consignes
%--------------------------------------------------------------
% fichier zuiedit_data_cons_limites_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de consignes
%
% syntaxe :
% 
% entrees :
%	action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 31/07/2001.
% 
% liste des modifications : 
%  * 27/09/2001 -> modification du libelle en x : texte_x = 'temps'
%
%--------------------------------------------------------------
function zuiedit_data_cons_limites_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_cons_limites') ;

% information pour l'assistant
zuicr(hfig,action) ;

% variables d'entrée de l'éditeur de consignes zuieditcons
x       = evalin('base','data.gene.temps') ;
texte_x = 'temps' ;
var_x   = 'void';
canal   = 1 ;

% selon ation
switch lower(action)

case 'edit_ip'
	nom     = 'data.cons.ip' ;
	y       = evalin('base','data.cons.ip') ;
	texte_y = 'Ip' ;
	var_y   = 'data.cons.ip';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_ip) ;

case 'import_ip'
	nom     = 'data.cons.ip' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_ip) ;
	
case 'edit_vloop'
	nom     = 'data.cons.vloop' ;
	y       = evalin('base','data.cons.vloop') ;
	texte_y = 'Vloop' ;
	var_y   = 'data.cons.vloop';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_vloop) ;

case 'import_vloop'
	nom     = 'data.cons.vloop' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_vloop) ;

case 'edit_flux'
	nom     = 'data.cons.flux' ;
	y       = evalin('base','data.cons.flux') ;
	texte_y = 'Flux' ;
	var_y   = 'data.cons.flux';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_flux) ;

case 'import_flux'
	nom     = 'data.cons.flux' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_flux) ;

case 'edit_ne1'
	nom     = 'data.cons.ne1' ;
	y       = evalin('base','data.cons.ne1') ;
	texte_y = 'Ne1' ;
	var_y   = 'data.cons.ne1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_ne1) ;

case 'import_ne1'
	nom     = 'data.cons.ne1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_ne1) ;

case 'edit_ge1'
	nom     = 'data.cons.ge1' ;
	y       = evalin('base','data.cons.ge1') ;
	texte_y = 'Ge1' ;
	var_y   = 'data.cons.ge1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_ge1) ;

case 'import_ge1'
	nom     = 'data.cons.ge1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_ge1) ;

case 'edit_te1'
	nom     = 'data.cons.te1' ;
	y       = evalin('base','data.cons.te1') ;
	texte_y = 'Te1' ;
	var_y   = 'data.cons.te1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_te1) ;

case 'import_te1'
	nom     = 'data.cons.te1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_te1) ;

case 'edit_qe1'
	nom     = 'data.cons.qe1' ;
	y       = evalin('base','data.cons.qe1') ;
	texte_y = 'Qe1' ;
	var_y   = 'data.cons.qe1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_qe1) ;

case 'import_qe1'
	nom     = 'data.cons.qe1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_qe1) ;

case 'edit_pe1'
	nom     = 'data.cons.pe1' ;
	y       = evalin('base','data.cons.pe1') ;
	texte_y = 'Pe1' ;
	var_y   = 'data.cons.pe1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_pe1) ;

case 'import_pe1'
	nom     = 'data.cons.pe1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_pe1) ;

case 'edit_ti1'
	nom     = 'data.cons.ti1' ;
	y       = evalin('base','data.cons.ti1') ;
	texte_y = 'Ti1' ;
	var_y   = 'data.cons.ti1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_ti1) ;

case 'import_ti1'
	nom     = 'data.cons.ti1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_ti1) ;

case 'edit_qi1'
	nom     = 'data.cons.qi1' ;
	y       = evalin('base','data.cons.qi1') ;
	texte_y = 'Qi1' ;
	var_y   = 'data.cons.qi1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_qi1) ;

case 'import_qi1'
	nom     = 'data.cons.qi1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_qi1) ;

case 'edit_pion1'
	nom     = 'data.cons.pion1' ;
	y       = evalin('base','data.cons.pion1') ;
	texte_y = 'Pion1' ;
	var_y   = 'data.cons.pion1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_pion1) ;

case 'import_pion1'
	nom     = 'data.cons.pion1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_pion1) ;

case 'edit_rot1'
	nom     = 'data.cons.rot1' ;
	y       = evalin('base','data.cons.rot1') ;
	texte_y = 'Rot1' ;
	var_y   = 'data.cons.rot1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_rot1) ;

case 'import_rot1'
	nom     = 'data.cons.rot1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_rot1) ;

case 'edit_frot1'
	nom     = 'data.cons.frot1' ;
	y       = evalin('base','data.cons.frot1') ;
	texte_y = 'Frot1' ;
	var_y   = 'data.cons.frot1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_frot1) ;

case 'import_frot1'
	nom     = 'data.cons.frot1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_frot1) ;

case 'edit_fe1'
	nom     = 'data.cons.fe1' ;
	y       = evalin('base','data.cons.fe1') ;
	texte_y = 'Fe1' ;
	var_y   = 'data.cons.fe1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_fe1) ;

case 'import_fe1'
	nom     = 'data.cons.fe1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_fe1) ;

case 'edit_fi1'
	nom     = 'data.cons.fi1' ;
	y       = evalin('base','data.cons.fi1') ;
	texte_y = 'Fi1' ;
	var_y   = 'data.cons.fi1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_fi1) ;

case 'import_fi1'
	nom     = 'data.cons.fi1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_fi1) ;

case 'edit_ffi1'
	nom     = 'data.cons.ffi1' ;
	y       = evalin('base','data.cons.ffi1') ;
	texte_y = 'FFi1' ;
	var_y   = 'data.cons.fi1';
	zuieditcons(nom,info.data.cons.ip,x,y,texte_x,texte_y,var_x,var_y,canal, ...
                    '',[],[],'','') ;
	zuireset(h.edit_ffi1) ;

case 'import_ffi1'
	nom     = 'data.cons.ffi1' ;
	zuiedit_import_mode(nom,'consigne') ;
	zuireset(h.import_ffi1) ;

case {'btn_quit','close'}
	zuicloseone(hfig);	
	
case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
otherwise
	warning('ation non prise en compte')
	   
end

