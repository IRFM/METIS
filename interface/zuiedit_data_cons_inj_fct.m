% ZUIDEDIT_DATA_CONS_INJ_FCT gestion callbacks du formulaire d'edition de consignes
%--------------------------------------------------------------
% fichier zuiedit_data_cons_inj_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de consignes
%
% syntaxe :
% 
% entrees :
%	action       =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 15/09/2003.
% 
% liste des modifications : 
%  * 13/09/2001 -> Pour les glacons, appel de l'editeur zuiedit_mode
%  * 27/09/2001 -> modification du libelle en x : texte_x = 'temps'
%  * 15/09/2003 -> ajout de ntnd
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
texte_x     = 'temps' ;
var_x       = 'void';
code_retour = '' ;
liste_ref   = {} ;
var_ref     = {} ;
texte_prop  = '' ;
var_prop    = '' ;

nbg = evalin('base','param.gene.nbg') ;
nb  = evalin('base','param.nombre.glacon') ;
num = '' ;
for i=1:max(nbg,nb)
	[t,r] = strtok(lower(action),num2str(i)) ;
	if ~isempty(r)
		action = t ;
		num    = r; 
		break
	end
end
if ~isempty(num)
	numm1  = num2str(str2num(num)-1);
end
% selon ation
switch lower(t)

% data.cons.c
% -----------
case 'edit_module_c'
	nom     = 'data.cons.c' ;
	texte_y = 'c' ;
	var_y   = nom ;
	canal   = 1 ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
		texte_y = strcat(texte_y,'(:,',num,')') ; 
		canal   = str2num(num) ;
	end
	y    = evalin('base',nom) ;
 	hout = zuieditcons(nom,info.data.cons.c,x,y,texte_x,texte_y,var_x,var_y,canal, ...
 	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(getfield(h,strcat('edit_module_c',num2str(num)))) ;
	
case 'prec_c'
	rep =questdlg('On recopie les valeurs du canal précédent ?', ...
	              'Confirmation', ...
	              'Oui','Non','Oui');
	switch rep
	case 'Oui'
		nom     = strcat('data.cons.c(:,',numm1,')') ;
		y       = evalin('base',nom) ;
		nom     = strcat('data.cons.c(:,',num,')') ;
		zassignin('base',nom,y) ;
	case 'Non'
	end
	zuireset(getfield(h,strcat('prec_c',num2str(num)))) ;

case 'import_c'
	nom     = 'data.cons.c' ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
	end
	hout = zuiedit_import_mode(nom,'consigne') ;
	zuireset(getfield(h,strcat('import_c',num2str(num)))) ;

% data.cons.pomp
% --------------
case 'edit_module_pomp'
	nom     = 'data.cons.pomp' ;
	texte_y = 'pomp' ;
	var_y   = nom ;
	canal   = 1;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
		texte_y = strcat(texte_y,'(:,',num,')') ;
		canal   = str2num(num) ;
	end
	y       = evalin('base',nom) ;
	hout = zuieditcons(nom,info.data.cons.pomp,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(getfield(h,strcat('edit_module_pomp',num2str(num)))) ;

case 'import_pomp'
	nom     = 'data.cons.pomp' ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
	end
	hout = zuiedit_import_mode(nom,'consigne') ;
	zuireset(getfield(h,strcat('import_pomp',num2str(num)))) ;

% data.cons.glacon
% ----------------
case 'edit_module_glacon'
	nom     = 'data.cons.glacon' ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
	end
	zuiedit_mode(nom,{0,1},{'attente','lancement'}) ;
	zuireset(getfield(h,strcat('edit_module_glacon',num2str(num)))) ;

case 'import_glacon'
	nom     = 'data.cons.glacon' ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
	end
	hout = zuiedit_import_mode(nom,'consigne') ;
	zuireset(getfield(h,strcat('import_glacon',num2str(num)))) ;
	
case 'prec_glacon'
	rep =questdlg('On recopie les valeurs du canal précédent ?', ...
	              'Confirmation', ...
	              'Oui','Non','Oui');
	switch rep
	case 'Oui'
		nom     = strcat('data.cons.glacon(:,',numm1,')') ;
		y       = evalin('base',nom) ;
		nom     = strcat('data.cons.glacon(:,',num,')') ;
		zassignin('base',nom,y) ;
	case 'Non'
	end
	zuireset(getfield(h,strcat('prec_glacon',num2str(num)))) ;

% data.cons.zeffm
% ---------------
case 'edit_module_zeffm'
	nom     = 'data.cons.zeffm' ;
	texte_y = 'zeffm' ;
	var_y   = nom ;
	canal   = 1 ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
		texte_y = strcat(texte_y,'(:,',num,')') ;
		canal   = str2num(num) ;
	end
	y       = evalin('base',nom) ;
	hout = zuieditcons(nom,info.data.cons.pomp,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   code_retour,liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(getfield(h,strcat('edit_module_zeffm',num2str(num)))) ;

case 'import_zeffm'
	nom     = 'data.cons.zeffm' ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
	end
	hout = zuiedit_import_mode(nom,'consigne') ;
	zuireset(getfield(h,strcat('import_zeffm',num2str(num)))) ;

% data.cons.nhnd
% ---------------
case 'edit_module_nhnd'
	nom     = 'data.cons.nhnd' ;
	texte_y = 'nhnd' ;
	var_y   = nom ;
	canal   = 1 ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
		texte_y = strcat(texte_y,'(:,',num,')') ;
		canal   = str2num(num) ;
	end
	y       = evalin('base',nom) ;
	hout = zuieditcons(nom,info.data.cons.pomp,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   'nhnd',liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(getfield(h,strcat('edit_module_nhnd',num2str(num)))) ;

case 'import_nhnd'
	nom     = 'data.cons.nhnd' ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
	end
	hout = zuiedit_import_mode(nom,'consigne') ;
	zuireset(getfield(h,strcat('import_nhnd',num2str(num)))) ;

% data.cons.ntnd
% ---------------
case 'edit_module_ntnd'
	nom     = 'data.cons.nhnd' ;
	texte_y = 'ntnd' ;
	var_y   = nom ;
	canal   = 1 ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
		texte_y = strcat(texte_y,'(:,',num,')') ;
		canal   = str2num(num) ;
	end
	y       = evalin('base',nom) ;
	hout = zuieditcons(nom,info.data.cons.pomp,x,y,texte_x,texte_y,var_x,var_y,canal, ...
	                   'ntnd',liste_ref,var_ref,texte_prop,var_prop) ;	                   
	zuireset(getfield(h,strcat('edit_module_ntnd',num2str(num)))) ;

case 'import_ntnd'
	nom     = 'data.cons.nhnd' ;
	if ~isempty(num)
		nom     = strcat(nom,'(:,',num,')') ;
	end
	hout = zuiedit_import_mode(nom,'consigne') ;
	zuireset(getfield(h,strcat('import_ntnd',num2str(num)))) ;

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
