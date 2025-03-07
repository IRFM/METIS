function zuiedit_source_pellet_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('edit_refval_pellet') ;

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
