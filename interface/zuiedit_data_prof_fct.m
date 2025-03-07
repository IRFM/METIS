% ZUIDEDIT_DATA_PROF_FCT gestion callbacks du formulaire d'edition de profils
%--------------------------------------------------------------------
% fichier zuiedit_data_prof_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de d'edition de profils 
%
% syntaxe :
%	zuiedit_data_prof_fct(action)
%
% entrees :
%	action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 1.7, du 22/11/2001.
% 
% liste des modifications : 
%  * 22/11/2001 -> remplacement de gene.k par gene.kmin pour kini
%
%--------------------------------------------------------------
function zuiedit_data_prof_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

info    = zinfo;

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('data_prof') ;

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


case 'edit_psi'
	nom        = 'data.prof.psi' ;
	y          = evalin('base','data.prof.psi') ;
	texte_y    = 'psi' ;
	var_y      = 'data.prof.psi';
	texte_prop = 'li' ;
	var_prop   = 'data.cons.asser.li' ;
	hout_edit  = zuieditcons(nom,info.data.prof.psi,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                                 code_retour,liste_ref,var_ref,texte_prop,var_prop) 
	zuireset(h.edit_psi) ;

case 'param_psi'
	nom        = 'data.prof.psi' ;
	text_modul = 'li' ;
	var_modul  = 'data.cons.asser.li' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_psi) ;
	
case 'import_psi'
	nom     = 'data.prof.psi' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_psi) ;
	
case 'dessin_psi'
	hout = zuiedit_profcmplx('data.prof.psi','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_psi)


case 'edit_pe'
   tepe = evalin('base','param.edit.tepe','''''');
	if ~isempty(tepe)
		disp('Te profile window open -> You cannot assess Pe at the same time');
		zuireset(h.edit_pe) ;
		return
	end
	nom        = 'data.prof.pe' ;
	y          = evalin('base','data.prof.pe') ;
	texte_y    = 'pe' ;
	var_y      = 'data.prof.pe';
	texte_prop = 'te0' ;
	var_prop   = 'data.cons.asser.te0' ;
	zuieditcons(nom,info.data.prof.pe,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_pe) ;

case 'edit_te'
   tepe = evalin('base','param.edit.tepe','''''');
	if ~isempty(tepe)
		disp('Pe profile window open -> You cannot assess Te at the same time');
		zuireset(h.edit_te) ;
		return
	end
	nom        = 'data.prof.te' ;
	y          = evalin('base','data.prof.te') ;
	texte_y    = 'Te' ;
	var_y      = 'data.prof.te';
	texte_prop = 'te0' ;
	var_prop   = 'data.cons.asser.te0' ;
	zuieditcons(nom,info.data.prof.te,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_te) ;

case 'param_pe'
   tepe = evalin('base','param.edit.tepe','''''');
	if ~isempty(tepe)
		disp('Te profile window open -> You cannot assess Pe at the same time');
		zuireset(h.param_pe) ;
		return
	end
	nom        = 'data.prof.pe' ;
	text_modul = 'te0' ;
	var_modul  = 'data.cons.asser.te0' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_pe) ;

case 'param_te'
   tepe = evalin('base','param.edit.tepe','''''');
	if ~isempty(tepe)
		disp('Pe profile window open -> You cannot assess Te at the same time');
		zuireset(h.param_te) ;
		return
	end
	nom        = 'data.prof.te' ;
	text_modul = 'te0' ;
	var_modul  = 'data.cons.asser.te0' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_te) ;

case 'import_pe'
   tepe = evalin('base','param.edit.tepe','''''');
	if ~isempty(tepe)
		disp('Te profile window open -> You cannot assess Pe at the same time');
		zuireset(h.import_pe) ;
		return
	end
	nom     = 'data.prof.pe' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_pe) ;

case 'import_te'
   tepe = evalin('base','param.edit.tepe','''''');
	if ~isempty(tepe)
		disp('Pe profile window open -> You cannot assess Te at the same time');
		zuireset(h.import_te) ;
		return
	end
	nom     = 'data.prof.te' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_te) ;

case 'dessin_pe'
   tepe = evalin('base','param.edit.tepe','''''');
	if ~isempty(tepe)
		disp('Te profile window open -> You cannot assess Pe at the same time');
		zuireset(h.dessin_pe)
		return
	end
	hout = zuiedit_profcmplx('data.prof.pe','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_pe)


case 'dessin_te'
   tepe = evalin('base','param.edit.tepe','''''');
	if ~isempty(tepe)
		disp('Pe profile window open -> You cannot assess Te at the same time');
		zuireset(h.dessin_te)
		return
	end
	hout = zuiedit_profcmplx('data.prof.te','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_te)


case 'edit_pion'
   tipion = evalin('base','param.edit.tipion','''''');
	if ~isempty(tipion)
		disp('Ti profile window open -> You cannot assess Pion at the same time');
		zuireset(h.edit_pion) ;
		return
	end
	nom        = 'data.prof.pion' ;
	y          = evalin('base','data.prof.pion','''''') ;
	texte_y    = 'pion' ;
	var_y      = 'data.prof.pion';
	texte_prop = 'ti0' ;
	var_prop   = 'data.cons.asser.ti0' ;
	zuieditcons(nom,info.data.prof.pion,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_pion) ;

case 'edit_ti'
   tipion = evalin('base','param.edit.tipion','''''');
	if ~isempty(tipion)
		disp('Pion profile window open -> You cannot assess Ti at the same time');
		zuireset(h.edit_ti) ;
		return
	end
	nom        = 'data.prof.ti' ;
	y          = evalin('base','data.prof.ti') ;
	texte_y    = 'ti' ;
	var_y      = 'data.prof.ti';
	texte_prop = 'ti0' ;
	var_prop   = 'data.cons.asser.ti0' ;
	zuieditcons(nom,info.data.prof.ti,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ti) ;

case 'param_pion'
   tipion = evalin('base','param.edit.tipion','''''');
	if ~isempty(tipion)
		disp('Ti profile window open -> You cannot assess Pion at the same time');
		zuireset(h.param_pion) ;
		return
	end
	nom        = 'data.prof.pion' ;
	text_modul = 'ti0' ;
	var_modul  = 'data.cons.asser.ti0' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_pion) ;

case 'param_ti'
   tipion = evalin('base','param.edit.tipion','''''');
	if ~isempty(tipion)
		disp('Pion profile window open -> You cannot assess Ti at the same time');
		zuireset(h.param_ti) ;
		return
	end
	nom        = 'data.prof.ti' ;
	text_modul = 'ti0' ;
	var_modul  = 'data.cons.asser.ti0' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_ti) ;

case 'import_pion'
   tipion = evalin('base','param.edit.tipion','''''');
	if ~isempty(tipion)
		disp('Ti profile window open -> You cannot assess Pion at the same time');
		zuireset(h.import_pion) ;
		return
	end
	nom     = 'data.prof.pion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_pion) ;

case 'import_ti'
   tipion = evalin('base','param.edit.tipion','''''');
	if ~isempty(tipion)
		disp('Pion profile window open -> You cannot assess Ti at the same time');
		zuireset(h.import_ti) ;
		return
	end
	nom     = 'data.prof.ti' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ti) ;

case 'dessin_pion'
   tipion = evalin('base','param.edit.tipion','''''');
	if ~isempty(tipion)
		disp('Ti profile window open -> You cannot assess Pion at the same time');
		zuireset(h.dessin_pion)
		return
	end
	hout = zuiedit_profcmplx('data.prof.pion','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_pion)

case 'dessin_ti'
   tipion = evalin('base','param.edit.tipion','''''');
	if ~isempty(tipion)
		disp('Pion profile window open -> You cannot assess ti at the same time');
		zuireset(h.dessin_ti)
		return
	end
	hout = zuiedit_profcmplx('data.prof.ti','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_ti)


case 'edit_ne'
	nom        = 'data.prof.ne' ;
	y          = evalin('base','data.prof.ne') ;
	texte_y    = 'ne' ;
	var_y      = 'data.prof.ne';
	texte_prop = 'nemoy' ;
	var_prop   = 'data.cons.asser.nemoy' ;
	zuieditcons(nom,info.data.prof.ne,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ne) ;

case 'param_ne'
	nom        = 'data.prof.ne' ;
	text_modul = 'nemoy' ;
	var_modul  = 'data.cons.asser.nemoy' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_ne) ;

case 'import_ne'
	nom     = 'data.prof.ne' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ne) ;

case 'dessin_ne'
	hout = zuiedit_profcmplx('data.prof.ne','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_ne)


case 'edit_rot'
	nom        = 'data.prof.rot' ;
	y          = evalin('base','data.prof.rot') ;
	texte_y    = 'rot' ;
	var_y      = 'data.prof.rot';
	texte_prop = 'ti0' ;
	var_prop   = 'data.cons.asser.ti0' ;
	zuieditcons(nom,info.data.prof.rot,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_rot) ;

case 'param_rot'
	nom        = 'data.prof.rot' ;
	text_modul = 'ti0' ;
	var_modul  = 'data.cons.asser.ti0' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_rot) ;

case 'import_rot'
	nom     = 'data.prof.rot' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_rot) ;

case 'dessin_rot'
	hout = zuiedit_profcmplx('data.prof.rot','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_rot)


case 'edit_fluce'
	nom        = 'data.prof.fluce' ;
	y          = evalin('base','data.prof.fluce') ;
	texte_y    = 'fluce' ;
	var_y      = 'data.prof.fluce';
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info.data.prof.fluce,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_fluce) ;

case 'param_fluce'
	nom        = 'data.prof.fluce' ;
	text_modul = '' ;
	var_modul  = 'data.gene.temps' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_fluce) ;

case 'import_fluce'
	nom     = 'data.prof.fluce' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_fluce) ;

case 'dessin_fluce'
	hout = zuiedit_profcmplx('data.prof.fluce','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_fluce)


case 'edit_flucion'
	nom        = 'data.prof.flucion' ;
	y          = evalin('base','data.prof.flucion') ;
	texte_y    = 'flucion' ;
	var_y      = 'data.prof.flucion';
	texte_prop = '' ;
	var_prop   = '' ;
	zuieditcons(nom,info.data.prof.flucion,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_flucion) ;

case 'param_flucion'
	nom        = 'data.prof.flucion' ;
	text_modul = '' ;
	var_modul  = 'data.gene.temps' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_flucion) ;

case 'import_flucion'
	nom     = 'data.prof.flucion' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_flucion) ;

case 'dessin_flucion'
	hout = zuiedit_profcmplx('data.prof.flucion','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_flucion)


case 'edit_jmoy'
	nom        = 'data.prof.jmoy' ;
	y          = evalin('base','data.prof.jmoy') ;
	texte_y    = 'jmoy' ;
	var_y      = 'data.prof.jmoy' ;
	texte_prop = 'ip' ;
	var_prop   = 'data.cons.ip' ;
	zuieditcons(nom,info.data.prof.jmoy,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_jmoy) ;

case 'param_jmoy'
	nom        = 'data.prof.jmoy' ;
	text_modul = 'ip' ;
	var_modul  = 'data.cons.ip' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_jmoy) ;

case 'import_jmoy'
	nom     = 'data.prof.jmoy' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_jmoy) ;

case 'dessin_jmoy'
	hout = zuiedit_profcmplx('data.prof.jmoy','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_jmoy)


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


case 'edit_xdur'
	nom        = 'data.prof.xdur' ;
	y          = evalin('base','data.prof.xdur') ;
	texte_y    = 'xdur' ;
	var_y      = 'data.prof.xdur';
	texte_prop = 'hyb' ;
	var_prop   = 'data.cons.hyb' ;
	zuieditcons(nom,info.data.prof.xdur,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_xdur) ;

case 'param_xdur'
	nom        = 'data.prof.xdur' ;
	text_modul = 'hyb' ;
	var_modul  = 'data.cons.hyb' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_xdur) ;

case 'import_xdur'
	nom     = 'data.prof.xdur' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_xdur) ;

case 'dessin_xdur'
	hout = zuiedit_profcmplx('data.prof.xdur','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_xdur)

case 'edit_vtor_exp'
	nom        = 'data.prof.vtor_exp' ;
	y          = evalin('base','data.prof.vtor_exp') ;
	texte_y    = 'vtor_exp' ;
	var_y      = 'data.prof.vtor_exp';
	texte_prop = 'hyb' ;
	var_prop   = 'data.cons.hyb' ;
	zuieditcons(nom,info.data.prof.vtor_exp,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_vtor_exp) ;

case 'param_vtor_exp'
	nom        = 'data.prof.vtor_exp' ;
	text_modul = 'hyb' ;
	var_modul  = 'data.cons.hyb' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_vtor_exp) ;

case 'import_vtor_exp'
	nom     = 'data.prof.vtor_exp' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_vtor_exp) ;

case 'dessin_vtor_exp'
	hout = zuiedit_profcmplx('data.prof.vtor_exp','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_vtor_exp)


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


case 'edit_ej'
	nom        = 'data.prof.ej' ;
	y          = evalin('base','data.prof.ej') ;
	texte_y    = 'ej' ;
	var_y      = 'data.prof.ej';
	texte_prop = 'vloop' ;
	var_prop   = 'data.cons.vloop' ;
	zuieditcons(nom,info.data.prof.ej,x,y(kini,:),texte_x,texte_y,var_x,var_y,canal, ...
                    code_retour,liste_ref,var_ref,texte_prop,var_prop) ;
	zuireset(h.edit_ej) ;

case 'param_ej'
	nom        = 'data.prof.ej' ;
	text_modul = 'vloop' ;
	var_modul  = 'data.cons.vloop' ;
	zuiparam(nom,text_modul,var_modul) ;
	zuireset(h.param_ej) ;

case 'import_ej'
	nom     = 'data.prof.ej' ;
	zuiedit_import_mode(nom,'profil') ;
	zuireset(h.import_ej) ;

case 'dessin_ej'
	hout = zuiedit_profcmplx('data.prof.ej','data.gene.temps','param.gene.x') ;
	zuireset(h.dessin_ej)


case {'btn_quit','close'}
	zuicloseone(hfig);

case 'init'
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);

otherwise
	warning('action not taking into account')

end

