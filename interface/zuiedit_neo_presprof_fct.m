function zuiedit_iondens_presprof_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_neo_presprof') ;

% information pour l'assistant
zuicr(hfig,action) ;

% variables d'entr� de l'�iteur de consignes zuieditcons
x           = evalin('base','param.gene.x') ;
kini        = evalin('base','param.gene.kmin') ;
texte_x     = 'x' ;
var_x       = 'void';
canal       = 1 ;
code_retour = '' ;
liste_ref   = {'Pe','Pion','Ne','Psi','Jmoy','\alphae','Zeff','empty'} ;
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
	warning('action not taking into account')

end

