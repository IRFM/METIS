% ZUIDEDIT_PARAM_MODEXT_FCT gestion callbacks du formulaire des modules externes
%--------------------------------------------------------------
% fichier zuiedit_param_modext_ctrl.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	des modules externes , parametres  g��aux
% 	sous le mode edition du formulaire principal
%
% syntaxe :
%	zuiedit_param_modext_fct(action)
%
% entrees :
%  action : tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 3.0, du 10/01/2005.
% 
% liste des modifications : 
%
%   * 12/09/2001  -> changement de la liste des options pour les glacons  (J-F Artaud)
%   * 05/11/2001  -> ajout de l'idn dans la validation
%   * 14/03/2002  -> ajout de la stabilite MHD
%   * 19/03/2002  -> correction bug mode mhd
%   * 11/12/2002  -> interface anglais
%   * 03/09/2003 -> ajout gestion module ripple
%   * 10/11/2003 -> ajout des modules "machine" et "post"
%   * 10/01/2005 -> ajout du module cyclo
%   * 13/04/2006 -> ajout de la mise a jour du mode de calcul du module cyclo - petit patch (FI)
%
%--------------------------------------------------------------
function zuiedit_param_modext_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('ed_param_modext');

% information pour l'assistant
zuicr(hfig,action) ;
try
   hoc = getfield(h,action) ;
catch
   hoc =[];
end
  
fonction = evalin('base','param.fonction') ;

% liste_coef
%liste_coef = {'zjetto     ', ...
%              'zbgbs      ', ...
%              'zcoefsimple', ...
%	      'zitg       ',...
%	      'zetg       ',...
%	      'zetg_hoang ',...
%	      'zetg_stable',...
%	      'zweiland   ',...
%	      'zself      ',...
%	      'zrlwcoef   ',...
%	      'other      '} ;
%
%liste_coef =  zlistcoef;
liste_module =zlist_module;
liste_coef =liste_module.coef;

% selon ation
switch lower(action)

% fonctions equilibre
% equil
case {'chng_equi'}
	% disp('changer equi')
	liste_fct = liste_module.equi; %{'zequi_helena','other'} ;
	hout = zuichgfct('equi',liste_fct,h.equi,h.fct_equi) ;
	zuireset(hoc) ;

case {'para_equi'}
	% disp('parametres equi') ;
	zuifuninterface('equi') ;
	zuireset(hoc) ;
		
case {'mode_equi'}
	% disp('changer mode equi')
 	value = zuidata(h.mode_equi) ;
	if value==3
		zuienable(h.edit_equi) ;
	else
		zuidisable(h.edit_equi) ;
	end
		
case {'edit_equi'}
	% disp('edition mode equi')
	zuiedit_mode('data.mode.equi',{0,1,2,3}, ...
	            {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;
		
% neo
case {'chng_neo'}
	liste_fct = liste_module.neo;    %{'zneofile','zneo','other'} ;
	hout = zuichgfct('neo',liste_fct,h.neo,h.fct_neo) ;
	zuireset(hoc) ;

case {'para_neo'}
	% disp('parametres neo') ;
	zuifuninterface('neo') ;
	zuireset(hoc) ;
		
case {'mode_neo'}
	% disp('changer mode neo')
 	value = zuidata(h.mode_neo) ;
	if value==3
		zuienable(h.edit_neo) ;
	else
		zuidisable(h.edit_neo) ;
	end
		
case {'edit_neo'}
	% disp('edition mode neo')
	zuiedit_mode('data.mode.neo',{1,2}, ...
	            {'read from data','calculated'}) ;
	zuireset(hoc) ;
		
% rip
case {'chng_rip'}
	liste_fct = liste_module.rip;    %{'zripple_therm','other'} ;
	hout = zuichgfct('rip',liste_fct,h.rip,h.fct_rip) ;
	zuireset(hoc) ;

case {'para_rip'}
	zuifuninterface('rip') ;
	zuireset(hoc) ;
		
case {'mode_rip'}
 	value = zuidata(h.mode_rip) ;
	if value==3
		zuienable(h.edit_rip) ;
	else
		zuidisable(h.edit_rip) ;
	end
		
case {'edit_rip'}
	zuiedit_mode('data.mode.rip',{0,1,2}, ...
	            {'set to zero','read from data','calculated'}) ;
	zuireset(hoc) ;
		
% fonctions de MHD
% dds
case {'chng_mhd_dds'}
	liste_fct = liste_module.mhd.dds; %{'zddscrash','other'} ;
	hout = zuichgfct('mhd.dds',liste_fct,h.mhd_dds,h.fct_mhd_dds) ;
	zuireset(hoc) ;

case {'para_mhd_dds'}
	% disp('parametres mhd_dds')
	zuifuninterface('mhd.dds') ;
	zuireset(hoc) ;
		
case {'mode_mhd_dds'}
	% disp('changer mode mhd_dds')
 	value = zuidata(h.mode_mhd_dds) ;
	if value==2
		zuienable(h.edit_mhd_dds) ;
	else
		zuidisable(h.edit_mhd_dds) ;
	end
		
case {'edit_mhd_dds'}
	% disp('edition mode mhd_dds')
	zuiedit_mode('data.mode.mhd.dds',{0,1},{'set to zero','calculated'}) ;
	zuireset(hoc) ;
		
% elm
case {'chng_mhd_elm'}
	% disp('changer elm')
	liste_fct = liste_module.mhd.elm; %{'','other'} ;
	hout = zuichgfct('mhd.elm',liste_fct,h.mhd_elm,h.fct_mhd_elm) ;
	zuireset(hoc) ;

case {'para_mhd_elm'}
	% disp('parametres mhd_elm')
	zuifuninterface('mhd.elm') ;
	zuireset(hoc) ;
		
case {'mode_mhd_elm'}
	% disp('changer mode mhd_elm')
 	value = zuidata(h.mode_mhd_elm) ;
	if value==2
		zuienable(h.edit_mhd_elm) ;
	else
		zuidisable(h.edit_mhd_elm) ;
	end
		
case {'edit_mhd_elm'}
	% disp('edition mode mhd_elm')
	zuiedit_mode('data.mode.mhd.elm',{0,1},{'set to zero','calculated'}) ;
	zuireset(hoc) ;

% limite
case {'chng_mhd_limite'}
	liste_fct = liste_module.mhd.limite; %{'zlimmhd','other'} ;
	hout = zuichgfct('mhd.limite',liste_fct,h.mhd_limite,h.fct_mhd_limite) ;
	zuireset(hoc) ;

case {'para_mhd_limite'}
	% disp('parametres mhd_limite')
	zuifuninterface('mhd.limite') ;
	zuireset(hoc) ;
		
case {'mode_mhd_limite'}
	% disp('changer mode mhd_limite')
 	value = zuidata(h.mode_mhd_limite) ;
	if value==2
		zuienable(h.edit_mhd_limite) ;
	else
		zuidisable(h.edit_mhd_limite) ;
	end
		
case {'edit_mhd_limite'}
	% disp('edition mode limite')
	zuiedit_mode('data.mode.mhd.limite',{0,1},{'set to zero','calculated'}) ;
	zuireset(hoc);

% stabilite
case {'chng_mhd_stab'}
	liste_fct = liste_module.mhd.stab; %{'zlimmhd','other'} ;
	hout = zuichgfct('mhd.stab',liste_fct,h.mhd_stab,h.fct_mhd_stab) ;
	zuireset(hoc) ;

case {'para_mhd_stab'}
	% disp('parametres mhd_stab')
	zuifuninterface('mhd.stab') ;
	zuireset(hoc) ;
		
case {'mode_mhd_stab'}
	% disp('changer mode mhd_stab')
 	value = zuidata(h.mode_mhd_stab) ;
	if value==3
		zuienable(h.edit_mhd_stab) ;
	else
		zuidisable(h.edit_mhd_stab) ;
	end
		
case {'edit_mhd_stab'}
	% disp('edition mode stab')
	zuiedit_mode('data.mode.mhd.stab',{0,1,2,3,4}, ...
	            {'set to zero','read from data','calculated','copy previous value',...
					'calcule en post-processing'}) ;
	zuireset(hoc);
		
% fonctions de sources
% fci
case {'chng_fci'}
	% disp('changer fci')
	liste_fct = liste_module.fci; %{'zfcifile','zfcisimple','other'} ;
	hout = zuichgfct('fci',liste_fct,h.fci,h.fct_fci) ;
	zuireset(hoc) ;

case {'para_fci'}
	% disp('parametres fci')
	zuifuninterface('fci') ;
	zuireset(hoc) ;
		
case {'mode_fci'}
	% disp('changer mode fci')
 	value = zuidata(h.mode_fci) ;
	if value==3
		zuienable(h.edit_fci) ;
	else
		zuidisable(h.edit_fci) ;
	end
		
case {'edit_fci'}
	% disp('edition mode fci')
	zuiedit_mode('data.mode.fci',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;
		
% fce
case {'chng_fce'}
	liste_fct = liste_module.fce; %{'zfcesimple','zremafile','other'} ;
	hout = zuichgfct('fce',liste_fct,h.fce,h.fct_fce) ;
	zuireset(hoc) ;

case {'para_fce'}
	% disp('parametres fce')
	zuifuninterface('fce') ;
	zuireset(hoc) ;
		
case {'mode_fce'}
	% disp('changer mode fce')
 	value = zuidata(h.mode_fce) ;
	if value==3
		zuienable(h.edit_fce) ;
	else
		zuidisable(h.edit_fce) ;
	end
		
case {'edit_fce'}
	% disp('edition mode fce')
	zuiedit_mode('data.mode.fce',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;
		
% hyb
case {'chng_hyb'}
	liste_fct = liste_module.hyb; %{'zhybsimple','zdelphe','zdke','other'} ;
	hout = zuichgfct('hyb',liste_fct,h.hyb,h.fct_hyb) ;
	zuireset(hoc) ;

case {'para_hyb'}
	% disp('parametres hyb')
	zuifuninterface('hyb') ;
	zuireset(hoc) ;
		
case {'mode_hyb'}
	% disp('changer mode hyb')
 	value = zuidata(h.mode_hyb) ;
	if value==3
		zuienable(h.edit_hyb) ;
	else
		zuidisable(h.edit_hyb) ;
	end
		
case {'edit_hyb'}
	% disp('edition mode hyb')
	zuiedit_mode('data.mode.hyb',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;
		
% idn
case {'chng_idn'}
	liste_fct = liste_module.idn; %{'zidnsimple','zsinbad2temps','zsinbadfile','zidnfast','other'} ;
	hout = zuichgfct('idn',liste_fct,h.idn,h.fct_idn) ;
	zuireset(hoc) ;

case {'para_idn'}
	% disp('parametres idn')
	zuifuninterface('idn') ;
	zuireset(hoc) ;
		
case {'mode_idn'}
	% disp('changer mode idn')
 	value = zuidata(h.mode_idn) ;
	if value==3
		zuienable(h.edit_idn) ;
	else
		zuidisable(h.edit_idn) ;
	end
		
case {'edit_idn'}
	% disp('edition mode idn')
	zuiedit_mode('data.mode.idn',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;
		
% n0
case {'chng_n0'}
	liste_fct = liste_module.n0; %{'zneutres','zjonas','other'} ;
	hout = zuichgfct('n0',liste_fct,h.n0,h.fct_n0) ;
	zuireset(hoc) ;

case {'para_n0'}
	% disp('parametres n0')
	zuifuninterface('n0') ;
	zuireset(hoc) ;
		
case {'mode_n0'}
	% disp('changer mode n0')
 	value = zuidata(h.mode_n0) ;
	if value==3
		zuienable(h.edit_n0) ;
	else
		zuidisable(h.edit_n0) ;
	end
		
case {'edit_n0'}
	% disp('edition mode n0')
	zuiedit_mode('data.mode.n0',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;
		
% bord
case {'chng_bord'}
	liste_fct = liste_module.bord; %{'zrecycle','zmurgaz','other'} ;
	hout = zuichgfct('bord',liste_fct,h.bord,h.fct_bord) ;
	zuireset(hoc) ;

case {'para_bord'}
	% disp('changer bord')
	zuifuninterface('bord') ;
	zuireset(hoc) ;
		
case {'mode_bord'}
	% disp('changer mode bord')
 	value = zuidata(h.mode_bord) ;
	if value==3
		zuienable(h.edit_bord) ;
	else
		zuidisable(h.edit_bord) ;
	end
		
case {'edit_bord'}
	% disp('edition mode bord')
	zuiedit_mode('data.mode.bord',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;

% glacon
case {'chng_glacon'}
	liste_fct = liste_module.glacon;%{'zglaquelc','other'} ;
	hout = zuichgfct('pellet',liste_fct,h.glacon,h.fct_glacon) ;
	zuireset(hoc) ;

case {'para_glacon'}
	% disp('changer glacon')
	zuifuninterface('pellet') ;
	zuireset(hoc) ;
		
case {'mode_glacon'}
	% disp('changer mode glacon')
 	value = zuidata(h.mode_glacon) ;
	if value > 1
		zuienable(h.edit_glacon) ;
	else
		zuidisable(h.edit_glacon) ;
	end
		
case {'edit_glacon'}
	% disp('edition mode glacon')
	zuiedit_mode('data.mode.glacon',{0,1},...
                    {'off','on'}) ;
	zuireset(hoc) ;

% fusion
case {'chng_fus'}
	% disp('changer fus')
	liste_fct = liste_module.fus; %{'zfusion','zfusion_spot','other'} ;
	hout = zuichgfct('fus',liste_fct,h.fus,h.fct_fus) ;
	zuireset(hoc) ;

case {'para_fus'}
	% disp('changer fus')
	zuifuninterface('fus') ;
	zuireset(hoc) ;
		
case {'mode_fus'}
	% disp('changer mode fus')
 	value = zuidata(h.mode_fus) ;
	if value==3
		zuienable(h.edit_fus) ;
	else
		zuidisable(h.edit_fus) ;
	end
		
case {'edit_fus'}
	% disp('edition mode fus')
	zuiedit_mode('data.mode.fus',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;

% cyclo
case {'chng_cyclo'}
	% disp('changer cyclo')
	liste_fct =liste_module.cyclo; %{'zalbajar','zcytran','autre'} ;
	hout = zuichgfct('cyclo',liste_fct,h.cyclo,h.fct_cyclo) ;
	zuireset(hoc) ;

case {'para_cyclo'}
	% disp('changer cyclo')
	zuifuninterface('cyclo') ;
	zuireset(hoc) ;
		
case {'mode_cyclo'}
	% disp('changer mode cyclo')
 	value = zuidata(h.mode_cyclo) ;
	if value==3
		zuienable(h.edit_cyclo) ;
	else
		zuidisable(h.edit_cyclo) ;
	end
		
case {'edit_cyclo'}
	% disp('edition mode cyclo')
	zuiedit_mode('data.mode.cyclo',{0,1,2,3},...
                    {'set to zero','read from data','calculated','k->k+1'}) ;
	zuireset(hoc) ;

% coefficients de transport
% coefa
case {'chng_coefa'}
	hout = zuichgfct('coefa',liste_coef,h.coefa,h.fct_coefa) ;
	zuireset(hoc) ;

case {'para_coefa'}
	% disp('changer coefa')
	zuifuninterface('coefa') ;
	zuireset(hoc) ;
		
case {'mode_coefa'}
	% disp('changer mode coefa')
 	value = zuidata(h.mode_coefa) ;
	if value==3
		zuienable(h.edit_coefa) ;
	else
		zuidisable(h.edit_coefa) ;
	end
		
case {'edit_coefa'}
	% disp('edition mode coefa')
	zuireset(hoc) ;
		
% coefb
case {'chng_coefb'}
	% disp('changer coefb')
	hout = zuichgfct('coefb',liste_coef,h.coefb,h.fct_coefb) ;
	zuireset(hoc) ;

case {'para_coefb'}
	% disp('changer coefb')
	zuifuninterface('coefb') ;
	zuireset(hoc) ;
		
case {'mode_coefb'}
	% disp('changer mode coefb')
 	value = zuidata(h.mode_coefb) ;
	if value==3
		zuienable(h.edit_coefb) ;
	else
		zuidisable(h.edit_coefb) ;
	end
		
case {'edit_coefb'}
	% disp('edition mode coefb')
	zuireset(hoc) ;
		
% coefc
case {'chng_coefc'}
	hout = zuichgfct('coefc',liste_coef,h.coefc,h.fct_coefc) ;
	zuireset(hoc) ;

case {'para_coefc'}
	% disp('changer coefc')
	zuifuninterface('coefc') ;
	zuireset(hoc) ;
		
case {'mode_coefc'}
	% disp('changer mode coefc')
 	value = zuidata(h.mode_coefc) ;
	if value==3
		zuienable(h.edit_coefc) ;
	else
		zuidisable(h.edit_coefc) ;
	end
		
case {'edit_coefc'}
	% disp('edition mode coefc')
	zuireset(hoc) ;
		
% coefd
case {'chng_coefd'}
	hout = zuichgfct('coefd',liste_coef,h.coefd,h.fct_coefd) ;
	zuireset(hoc) ;

case {'para_coefd'}
	% disp('changer coefd')
	zuifuninterface('coefd') ;
	zuireset(hoc) ;
		
case {'mode_coefd'}
	% disp('changer mode coefd')
 	value = zuidata(h.mode_coefd) ;
	if value==3
		zuienable(h.edit_coefd) ;
	else
		zuidisable(h.edit_coefd) ;
	end
		
case {'edit_coefd'}
	% disp('edition mode coefd')
	zuireset(hoc) ;
		
% coefe
case {'chng_coefe'}
	hout = zuichgfct('coefe',liste_coef,h.coefe,h.fct_coefe) ;
	zuireset(hoc) ;

case {'para_coefe'}
	% disp('changer coefe')
	zuifuninterface('coefe') ;
	zuireset(hoc) ;
		
case {'mode_coefe'}
	% disp('changer mode coefe')
 	value = zuidata(h.mode_coefe) ;
	if value==3
		zuienable(h.edit_coefe) ;
	else
		zuidisable(h.edit_coefe) ;
	end
		
case {'edit_coefe'}
	% disp('edition mode coefe')
	zuireset(hoc) ;
		
% coeff
case {'chng_coeff'}
	hout = zuichgfct('coeff',liste_coef,h.coeff,h.fct_coeff) ;
	zuireset(hoc) ;

case {'para_coeff'}
	% disp('changer coeff')
	zuifuninterface('coeff');
	zuireset(hoc);
		
case {'mode_coeff'}
	% disp('changer mode coeff')
 	value = zuidata(h.mode_coeff) ;
	if value==3
		zuienable(h.edit_coeff) ;
	else
		zuidisable(h.edit_coeff) ;
	end
		
case {'edit_coeff'}
	% disp('edition mode coeff')
	zuireset(hoc) ;
		
% divers
% impur
case {'chng_impur'}
	% disp('changer impur')
	liste_fct =liste_module.impur; %{'zinebcompo','zitc','zinebcompovb','other'} ;
	hout = zuichgfct('impur',liste_fct,h.impur,h.fct_impur) ;
	zuireset(hoc) ;

case {'para_impur'}
	% disp('changer impur')
	zuifuninterface('impur') ;
	zuireset(hoc) ;
		
case {'mode_impur'}
	% disp('changer mode impur')
 	value = zuidata(h.mode_impur) ;
	if value==3
		zuienable(h.edit_impur) ;
	else
		zuidisable(h.edit_impur) ;
	end
		
case {'edit_impur'}
	% disp('edition mode impur')
	zuiedit_mode('data.mode.impur',{0,1,2,3},...
                    {'set to zero','read from data','calculated','copy previous value'}) ;
	zuireset(hoc) ;

% plot
case {'chng_plot'}
	% disp('changer plot')
	liste_fct = liste_module.plot; %{'zplot','other'} ;
	hout = zuichgfct('plot',liste_fct,h.plot,h.fct_plot) ;
	zuireset(hoc) ;

case {'para_plot'}
	% disp('changer plot')
	zuifuninterface('plot');
	zuireset(hoc);
		
case {'mode_plot'}
	% disp('changer mode plot')
 	value = zuidata(h.mode_plot) ;
	if value==3
		zuienable(h.edit_plot) ;
	else
		zuidisable(h.edit_plot) ;
	end
		
case {'edit_plot'}
	% disp('edition mode plot')
	zuiedit_mode('data.mode.plot',{0,1},...
                    {'off','on'}) ;
	zuireset(hoc);

% asser
case {'chng_asser'}
	% disp('changer asser')
	liste_fct = liste_module.asser; %{'zasserip','other'} ;
	hout = zuichgfct('asser',liste_fct,h.asser,h.fct_asser) ;
	zuireset(hoc) ;

case {'para_asser'}
	% disp('changer plot')
	zuifuninterface('asser');
	zuireset(hoc);
		
case {'mode_asser'}
	% disp('changer mode plot')
 	value = zuidata(h.mode_asser) ;
	if value==3
		zuienable(h.edit_asser) ;
	else
		zuidisable(h.edit_asser) ;
	end
		
case {'edit_asser'}
	% disp('edition mode plot')
	zuiedit_mode('data.mode.asser',{1,2},...
                    {'read from data','calculated'}) ;
	zuireset(hoc);


% post
case {'chng_post'}
	liste_fct =  liste_module.post; %zlistpost;
	hout = zuichgfct('post',liste_fct,h.post,h.fct_post) ;
	zuireset(hoc) ;

case {'para_post'}
	zuifuninterface('post') ;
	zuireset(hoc) ;
		
case {'mode_post'}
	value = zuidata(h.mode_post) ;
	if value==2
		zuienable(h.edit_post) ;
	else
		zuidisable(h.edit_post) ;
	end
		
case {'edit_post'}
	zuiedit_mode('data.mode.post',{0,2,3},{'no calculation','calculated','copy previous value'}) ;
	zuireset(hoc) ;

% machine
case {'chng_machine'}
	liste_fct =   liste_module.machine; %zlistdevice;
	hout = zuichgfct('machine',liste_fct,h.machine,h.fct_machine) ;
	zuireset(hoc) ;

case {'para_machine'}
	zuifuninterface('machine') ;
	zuireset(hoc) ;
		


% RAZ
case {'init','raz'}
	if ishandle(hfig)
		zuiformvisible(hfig);
	else
		zuiedit_param_modext;
		[hfig,h] = zuiformhandle('ed_param_modext');
	end
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;
	zuidisable(h.edit_equi) ;
	zuidisable(h.edit_neo) ;
	zuidisable(h.edit_rip) ;
	zuidisable(h.edit_mhd_dds) ;
	zuidisable(h.edit_mhd_elm) ;
	zuidisable(h.edit_mhd_limite) ;
	zuidisable(h.edit_mhd_stab) ;
	zuidisable(h.edit_fci) ;
	zuidisable(h.edit_fce) ;
	zuidisable(h.edit_hyb) ;
	zuidisable(h.edit_n0) ;
	zuidisable(h.edit_bord) ;
	zuidisable(h.edit_glacon) ;
	zuidisable(h.edit_fus) ;
	zuidisable(h.edit_impur) ;
	zuidisable(h.edit_plot) ;
	zuidisable(h.edit_asser) ;
	zuidisable(h.edit_post) ;
	zuidisable(h.edit_machine) ;
	zuireset(h.raz) ;
	
% Annulation
case {'annulation','close'}
		zuicloseone(hfig);	

% Validation
case 'validation'
	zuiformcache(hfig) ;
	zuidownloadform(hfig);
	var = evalin('base','data.mode.equi') ;
	val = zuidata(h.mode_equi) ;
	if val~= 3
		zassignin('base','data.mode.equi',ones(size(var))*val) ;
	end			

	var = evalin('base','data.mode.neo') ;
	val = zuidata(h.mode_neo) ;
	if val~= 3
		zassignin('base','data.mode.neo',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.rip') ;
	val = zuidata(h.mode_rip) ;
	if val~= 3
		zassignin('base','data.mode.rip',ones(size(var))*val) ;
	end			

	var = evalin('base','data.mode.mhd.dds') ;
	val = zuidata(h.mode_mhd_dds) ;
	if val~= 3
		zassignin('base','data.mode.mhd.dds',ones(size(var))*val) ;
	end			

	var = evalin('base','data.mode.mhd.elm') ;
	val = zuidata(h.mode_mhd_elm) ;
	if val~= 3
		zassignin('base','data.mode.mhd.elm',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.mhd.limite') ;
	val = zuidata(h.mode_mhd_limite) ;
	if val~= 3
		zassignin('base','data.mode.mhd.limite',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.mhd.stab') ;
	val = zuidata(h.mode_mhd_stab) ;
	if val~= 3
		zassignin('base','data.mode.mhd.stab',ones(size(var))*val) ;
	end			

		var = evalin('base','data.mode.fci') ;
	val = zuidata(h.mode_fci) ;
	if val~= 3
		zassignin('base','data.mode.fci',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.fce') ;
	val = zuidata(h.mode_fce) ;
	if val~= 3
		zassignin('base','data.mode.fce',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.hyb') ;
	val = zuidata(h.mode_hyb) ;
	if val~= 3
		zassignin('base','data.mode.hyb',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.idn') ;
	val = zuidata(h.mode_idn) ;
	if val~= 3
		zassignin('base','data.mode.idn',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.n0') ;
	val = zuidata(h.mode_n0) ;
	if val~= 3
		zassignin('base','data.mode.n0',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.bord') ;
	val = zuidata(h.mode_bord) ;
	if val~= 3
		zassignin('base','data.mode.bord',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.glacon') ;
	val = zuidata(h.mode_glacon) ;
	if val <= 1
		zassignin('base','data.mode.glacon',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.fus') ;
	val = zuidata(h.mode_fus) ;
	if val~= 3
		zassignin('base','data.mode.fus',ones(size(var))*val) ;
	end			

	var = evalin('base','data.mode.cyclo') ;
	val = zuidata(h.mode_cyclo) ;
	if val~= 3
		zassignin('base','data.mode.cyclo',ones(size(var))*val) ;
	end			

	var = evalin('base','data.mode.impur') ;
	val = zuidata(h.mode_impur) ;
	if val~= 3
		zassignin('base','data.mode.impur',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.plot') ;
	val = zuidata(h.mode_plot) ;
	if val~= 3
		zassignin('base','data.mode.plot',ones(size(var))*val) ;
	end			
		
	var = evalin('base','data.mode.asser') ;
	val = zuidata(h.mode_asser) ;
	if val~= 3
		zassignin('base','data.mode.asser',ones(size(var))*val) ;
	end
			
	var = evalin('base','data.mode.post') ;
	val = zuidata(h.mode_post) ;
	if val~= 2
		zassignin('base','data.mode.post',ones(size(var))*(val.*2)) ;
	end			
		
	zuisavenonok;
	zuireset(h.validation);
	
otherwise
	warning('action not taken into account')
	   
end

%% fonction donnant la liste des fonctions post
%function liste = zlistpost
%
%try
%   nom            = evalin('base','param.from.machine');
%   fonctionpost   = strcat('zpost_',lower(nom),'_new');
%   if ~isempty(which(fonctionpost))
%      liste          = {fonctionpost,'other'};
%   else
%      liste = {'other'}
%   end
%catch
%   liste = {'other'}
%end

%% fonction donnant la liste des fonctions machine
%function liste = zlistdevice
%
% try
%   nom            = evalin('base','param.from.machine');
%   fonctionpost   = strcat('zdevice_',lower(nom));
%   if ~isempty(which(fonctionpost))
%      liste          = {fonctionpost,'other'};
%   else
%      liste = {'other'}
%   end
%catch
%   liste = {'other'}
%end
