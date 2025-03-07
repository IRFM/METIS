% ZUIDEDIT_DATA_CONS_CHAUF_FCT  gestion des callbacks du formulaire d'edition de consignes
%--------------------------------------------------------------
% fichier zuiedit_data_cons_chauf_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de consignes
%
% syntaxe :
%	zuiedit_data_cons_chauf_fct(action)
%
% entrees :
%	action       =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 28/01/2004.
% 
% liste des modifications : 
% * 27/09/2001 -> modification du libelle en x : texte_x = 'temps'
% * 28/03/2002 -> correction pour Idn
% * 15/09/2003 -> edition angle poloidal et toroidal pour fce
% * 28/01/2004 -> suppression de la consigne dPdt pour IDN (plus utilisee)
% * 18/04/2006 -> passage en anglais
%--------------------------------------------------------------
function zuiedit_data_cons_chauf_hyb_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_cons_chauf_hyb') ;

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

nbfci = evalin('base','param.nombre.fci') ;
nbfce = evalin('base','param.nombre.fce') ;
nbhyb = evalin('base','param.nombre.hyb') ;
nbidn = evalin('base','param.nombre.idn') ;
a = [nbfci,nbfce,nbhyb,nbidn] ;
nb = max(sort(a)) ;
num ='';
for i=1:nb
[t,r] = strtok(lower(action),num2str(i)) ;
	if ~isempty(r)
		action = t ;
		num    = r ;
		break
	end
end
if ~isempty(num)
	numm1  = num2str(str2num(num)-1);
end

% selon ation
switch lower(t)

case 'edit_module_fci'
	nom     = 'data.cons.fci' ;
	canal   = 1 ;
	texte_y = 'fci' ;
	var_y   = nom ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
		texte_y = strcat(texte_y,'(:,',num,')') ;
		canal   = str2num(num) ;
	end
	y       = evalin('base',nom) ;
	code_retour = 'abs' ;
	hout = zuieditcons(nom,info.data.cons.fci,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(getfield(h,strcat('edit_module_fci',num2str(num)))) ;

case 'edit_module_fce'
	nom     = 'data.cons.fce' ;
	canal   = 1 ;
	texte_y = 'fce' ;
	var_y   = nom ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
		texte_y = strcat(texte_y,'(:,',num,')') ; 
		canal   = str2num(num) ;
	end
	y       = evalin('base',nom) ;
	code_retour = 'abs' ;
	hout = zuieditcons(nom,info.data.cons.fci,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   'pfce',liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(getfield(h,strcat('edit_module_fce',num2str(num)))) ;

case 'edit_module_hyb'
	nom     = 'data.cons.hyb' ;
	texte_y = 'hyb' ;
	canal   = 1 ;
	var_y       = nom ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
		texte_y = strcat(texte_y,'(:,',num,')') ; 
		canal   = str2num(num) ;
	end
	y           = evalin('base',nom) ;
	code_retour = 'abs' ;
	hout = zuieditcons(nom,info.data.cons.hyb,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(getfield(h,strcat('edit_module_hyb',num2str(num)))) ;

case 'edit_module_idn'
	nom     = 'data.cons.idn' ;
	texte_y = 'idn' ;
	canal   = 1 ;
	var_y   = nom ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
		texte_y = strcat(texte_y,'(:,',num,')') ; 
		canal   = str2num(num) ;
	end
	y       = evalin('base',nom) ;
	code_retour = 'idn' ;
	hout = zuieditcons(nom,info.data.cons.idn,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(getfield(h,strcat('edit_module_idn',num2str(num)))) ;

case 'edit_phase_fci'
	nom     = 'data.cons.fci' ;
	texte_y = 'fci' ;
	var_y   = nom ;
	canal   = 1 ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
		texte_y = strcat(texte_y,'(:,',num,')') ; 
		canal   = str2num(num) ;
	end
	y           = evalin('base',nom) ;
	code_retour = 'degres' ;
	hout = zuieditcons(nom,info.data.cons.fci,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(getfield(h,strcat('edit_phase_fci',num2str(num)))) ;

case 'edit_phase_fce'
	nom     = 'data.cons.fce' ;
	texte_y = 'fce' ;
	var_y   = nom ;
	canal   = 1 ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
		texte_y = strcat(texte_y,'(:,',num,')') ; 
		canal   = str2num(num) ;
	end
	y           = evalin('base',nom) ;
        
	%code_retour = 'degres' ;
	hout = zuieditcons(strcat('toro_',nom),info.data.cons.fce,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   'toro',liste_ref,var_ref,texte_prop,var_prop) ;	                   
	hout = zuieditcons(strcat('polo_',nom),info.data.cons.fce,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   'polo',liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(getfield(h,strcat('edit_phase_fce',num2str(num)))) ;

case 'edit_phase_hyb'
	nom     = 'data.cons.hyb' ;
	texte_y = 'hyb' ;
	var_y   = nom ;
	canal   = 1 ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
		texte_y = strcat(texte_y,'(:,',num,')') ; 
		canal   = str2num(num) ;
	end
	y           = evalin('base',nom) ;
	code_retour = 'angle' ;
	hout = zuieditcons(nom,info.data.cons.hyb,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(getfield(h,strcat('edit_phase_hyb',num2str(num)))) ;

case 'edit_phase_idn'
%	nom     = 'data.cons.idn' ;
%	texte_y = 'idn' ;
%	var_y   = nom 
%	canal   = 1 ;
%	if ~isempty(num)
%		nom     = strcat(nom,'(:,',num,')') ;
%		texte_y = strcat(texte_y,'(:,',num,')') ; 
%		canal   = str2num(num) ;
%	end
%	y           = evalin('base',nom) ;
%	code_retour = 'didndt' ;
%	hout = zuieditcons(nom,info.data.cons.idn,x,y,texte_x,texte_y,var_x,var_y,canal, ...
%	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(getfield(h,strcat('edit_phase_idn',num2str(num)))) ;

case 'import_fci'
	nom     = 'data.cons.fci' ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
	end
	hout = zuiedit_import_mode(nom,'consigne') ;
	zuireset(getfield(h,strcat('import_fci',num2str(num)))) ;

case 'import_fce'
	nom     = 'data.cons.fce' ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
	end
	hout = zuiedit_import_mode(nom,'consigne') ;
	zuireset(getfield(h,strcat('import_fce',num2str(num)))) ;

case 'import_hyb'
	nom     = 'data.cons.hyb' ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
	end
	hout = zuiedit_import_mode(nom,'consigne') ;
	zuireset(getfield(h,strcat('import_hyb',num2str(num)))) ;

case 'import_idn'
	nom     = 'data.cons.idn' ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
	end
	hout = zuiedit_import_mode(nom,'consigne') ;
	zuireset(getfield(h,strcat('import_idn',num2str(num)))) ;

case 'prec_fci'
	rep =questdlg('Copy values from previous channel ?', ...
	              'Confirmation', ...
	              'Yes','No','Yes') ;
	switch rep
	case 'Yes'
		nom     = strcat('data.cons.fci(:,',numm1,')') ;
		y       = evalin('base',nom) ;
		nom     = strcat('data.cons.fci(:,',num,')') ;
		zassignin('base',nom,y) ;
	case 'No'
	end
	zuireset(getfield(h,strcat('prec_fci',num2str(num)))) ;

case 'prec_fce'
	rep =questdlg('Copy values from previous channel ?', ...
	              'Confirmation', ...
	              'Yes','No','Yes') ;
	switch rep
	case 'Yes'
		nom     = strcat('data.cons.fce(:,',numm1,')') ;
		y       = evalin('base',nom) ;
		nom     = strcat('data.cons.fce(:,',num,')') ;
		zassignin('base',nom,y) ;
	case 'No'
	end
	zuireset(getfield(h,strcat('prec_fce',num2str(num)))) ;

case 'prec_hyb'
	rep =questdlg('Copy values from previous channel ?', ...
	              'Confirmation', ...
	              'Yes','No','Yes') ;
	switch rep
	case 'Yes'
		nom     = strcat('data.cons.hyb(:,',numm1,')') ;
		y       = evalin('base',nom) ;
		nom     = strcat('data.cons.hyb(:,',num,')') ;
		zassignin('base',nom,y) ;
	case 'No'
	end
	zuireset(getfield(h,strcat('prec_hyb',num2str(num)))) ;

case 'prec_idn'
	rep =questdlg('Copy values from previous channel ?', ...
	              'Confirmation', ...
	              'Yes','No','Yes') ;
	switch rep
	case 'Yes'
		nom     = strcat('data.cons.idn(:,',numm1,')') ;
		y       = evalin('base',nom) ;
		nom     = strcat('data.cons.idn(:,',num,')') ;
		zassignin('base',nom,y) ;
	case 'No'
	end
	zuireset(getfield(h,strcat('prec_idn',num2str(num)))) ;

case {'btn_quit','close'}
	zuicloseone(hfig) ;	
	
case 'init'
	zuiformvisible(hfig) ;
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;
	
otherwise
	warning('action non prise en compte')
end

   
end

