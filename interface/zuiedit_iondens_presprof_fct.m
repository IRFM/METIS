function zuiedit_iondens_presprof_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_iondens_presprof') ;

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

case 'edit_ae'
	nom        = 'data.prof.ae' ;
	y          = evalin('base','data.prof.ae') ;
	texte_y    = 'ae' ;
	var_y      = 'data.prof.ae';
	texte_prop = 'zeffm' ;
	var_prop   = 'data.cons.zeffm' ;
	zuieditcons(nom,info.data.prof.ae,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ae) ;

case 'param_ae'
	nom        = 'data.prof.ae' ;
	text_modul = 'zeffm' ;
	var_modul  = 'data.cons.zeffm' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_ae) ;

case 'import_ae'
	nom     = 'data.prof.ae' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ae) ;

case 'dessin_ae'
	hout = zuiedit_profcmplx('data.prof.ae','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_ae)


case 'edit_zeff'
	nom        = 'data.prof.zeff' ;
	y          = evalin('base','data.prof.zeff') ;
	texte_y    = 'zeff' ;
	var_y      = 'data.prof.zeff';
	texte_prop = 'zeffm' ;
	var_prop   = 'data.cons.zeffm' ;
	zuieditcons(nom,info.data.prof.zeff,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_zeff) ;

case 'param_zeff'
	nom        = 'data.prof.zeff' ;
	text_modul = 'zeffm' ;
	var_modul  = 'data.cons.zeffm' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_zeff) ;

case 'import_zeff'
	nom     = 'data.prof.zeff' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_zeff) ;

case 'dessin_zeff'
	hout = zuiedit_profcmplx('data.prof.zeff','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_zeff)

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

case {'btn_quit','close'}
	zuicloseone(hfig);

case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);

otherwise
	warning('action not taking into account')

end

