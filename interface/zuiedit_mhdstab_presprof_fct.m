function zuiedit_mhdstab_presprof_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_mhdstab_presprof') ;

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

case 'edit_mhd_cd'
	nom        = 'data.prof.mhd_cd' ;
	y          = evalin('base','data.prof.mhd_cd') ;
	texte_y    = 'mhd_cd' ;
	var_y      = 'data.prof.mhd_cd';
	texte_prop = 'q0' ;
	var_prop   = 'data.cons.asser.q0' ;
	zuieditcons(nom,info.data.prof.mhd_cd,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_mhd_cd) ;

case 'param_mhd_cd'
	nom        = 'data.prof.mhd_cd' ;
	text_modul = 'q0' ;
	var_modul  = 'data.cons.asser.q0' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_mhd_cd) ;

case 'import_mhd_cd'
	nom     = 'data.prof.mhd_cd' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_mhd_cd') ;

case 'dessin_mhd_cd'
	hout = zuiedit_profcmplx('data.prof.mhd_cd','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_mhd_cd)

case {'btn_quit','close'}
	zuicloseone(hfig);

case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);

otherwise
	warning('action not taking into account')

end

