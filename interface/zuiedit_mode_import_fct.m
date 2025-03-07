% ZUIDEDIT_MODE_IMPORT_FCT gestion des callbacks du formulaire importation valeurs d'edition de mode
%--------------------------------------------------------------
% fichier zuiedit_mode_import_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks  de zuiedit_mode_import
%
% syntaxe :
%	zuiedit_mode_import_fct(action)
%
% entrees :
%	action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 18/09/2003.
% 
% liste des modifications : 
%   * 14/03/2002 -> ajout de la stabilite MHD
%   * 18/09/2003 -> ajout mode.rotc
%
%--------------------------------------------------------------
function zuiedit_mode_import_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% disp('callback de zuiedit_mode_import_fct: ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('mode_import') ;

% information pour l'assistant
zuicr(hfig,action) ;

% handle du formulaire appelant
hfig_pere = getappdata(hfig,'pere') ;


% selon ation
switch lower(action)

case 'push_impur'
	y = evalin('base','data.mode.impur') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_impur);
	zuicloseone(hfig) ;

case 'push_psi'
	y = evalin('base','data.mode.psi') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_psi);
	zuicloseone(hfig) ;

case 'push_nel'
	y = evalin('base','data.mode.nel') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_nel);
	zuicloseone(hfig) ;

case 'push_pe'
	y = evalin('base','data.mode.pe') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_pe);
	zuicloseone(hfig) ;

case 'push_pion'
	y = evalin('base','data.mode.pion') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_pion);
	zuicloseone(hfig) ;

case 'push_equi'
	y = evalin('base','data.mode.equi') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_equi);
	zuicloseone(hfig) ;

case 'push_neo'
	y = evalin('base','data.mode.neo') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_neo);
	zuicloseone(hfig) ;

case 'push_fluce'
	y = evalin('base','data.mode.fluce') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_fluce);
	zuicloseone(hfig) ;

case 'push_flucion'
	y = evalin('base','data.mode.flucion') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_flucion);
	zuicloseone(hfig) ;

case 'push_rot'
	y = evalin('base','data.mode.rot') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_rot);
	zuicloseone(hfig) ;

case 'push_cons_psi'
	y = evalin('base','data.mode.cons.psi') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_cons_psi);
	zuicloseone(hfig) ;

case 'push_cons_ne'
	y = evalin('base','data.mode.cons.ne') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_cons_ne);
	zuicloseone(hfig) ;

case 'push_cons_pe'
	y = evalin('base','data.mode.cons.pe') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_cons_pe);
	zuicloseone(hfig) ;

case 'push_cons_pion'
	y = evalin('base','data.mode.cons.pion') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_cons_pion);
	zuicloseone(hfig) ;

case 'push_cons_fluce'
	y = evalin('base','data.mode.cons.fluce') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_cons_fluce);
	zuicloseone(hfig) ;

case 'push_cons_flucion'
	y = evalin('base','data.mode.cons.flucion') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_cons_flucion);
	zuicloseone(hfig) ;

case 'push_cons_rot'
	y = evalin('base','data.mode.cons.rot') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_cons_rot);
	zuicloseone(hfig) ;

case 'push_cons_zeffm'
	y = evalin('base','data.mode.cons.zeffm') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_cons_zeffm);
	zuicloseone(hfig) ;

case 'push_bord_ne'
	y = evalin('base','data.mode.consbord.ne') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_bord_ne);
	zuicloseone(hfig) ;

case 'push_bord_te'
	y = evalin('base','data.mode.consbord.te') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_bord_te);
	zuicloseone(hfig) ;

case 'push_bord_ti'
	y = evalin('base','data.mode.consbord.ti') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_bord_ti);
	zuicloseone(hfig) ;

case 'push_bord_fluxge'
	y = evalin('base','data.mode.consbord.fluxge') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_bord_fluxge);
	zuicloseone(hfig) ;

case 'push_bord_fluxqe'
	y = evalin('base','data.mode.consbord.fluxqe') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_bord_fluxqe);
	zuicloseone(hfig) ;

case 'push_bord_fluxqi'
	y = evalin('base','data.mode.consbord.fluxqi') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_bord_fluxqi);
	zuicloseone(hfig) ;

case 'push_mhd_dds'
	y = evalin('base','data.mode.mhd.dds') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_mhd_dds);
	zuicloseone(hfig) ;

case 'push_mhd_elm'
	y = evalin('base','data.mode.mhd.elm') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_mhd_elm);
	zuicloseone(hfig) ;

case 'push_mhd_limite'
	y = evalin('base','data.mode.mhd.limite') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_mhd_limite);
	zuicloseone(hfig) ;

case 'push_mhd_stab'
	y = evalin('base','data.mode.mhd.stab') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_mhd_stab);
	zuicloseone(hfig) ;

case 'push_fci'
	y = evalin('base','data.mode.fci') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_fci);
	zuicloseone(hfig) ;

case 'push_fce'
	y = evalin('base','data.mode.fce') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_fce);
	zuicloseone(hfig) ;

case 'push_hyb'
	y = evalin('base','data.mode.hyb') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_hyb);
	zuicloseone(hfig) ;

case 'push_idn'
	y = evalin('base','data.mode.idn') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_idn);
	zuicloseone(hfig) ;

case 'push_n0'
	y = evalin('base','data.mode.n0') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_n0);
	zuicloseone(hfig) ;

case 'push_bord'
	y = evalin('base','data.mode.bord') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_bord);
	zuicloseone(hfig) ;

case 'push_glacon'
	y = evalin('base','data.mode.glacon') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_glacon);
	zuicloseone(hfig) ;

case 'push_fus'
	y = evalin('base','data.mode.fus') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_fus);
	zuicloseone(hfig) ;

case 'push_ohm'
	y = evalin('base','data.mode.ohm') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ohm);
	zuicloseone(hfig) ;

case 'push_qneo'
	y = evalin('base','data.mode.qneo') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_qneo);
	zuicloseone(hfig) ;

case 'push_qei'
	y = evalin('base','data.mode.qei') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_qei);
	zuicloseone(hfig) ;

case 'push_prad'
	y = evalin('base','data.mode.prad') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_prad);
	zuicloseone(hfig) ;

case 'push_brem'
	y = evalin('base','data.mode.brem') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_brem);
	zuicloseone(hfig) ;

case 'push_cyclo'
	y = evalin('base','data.mode.cyclo') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_cyclo);
	zuicloseone(hfig) ;

case 'push_ee'
	y = evalin('base','data.mode.ee') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ee);
	zuicloseone(hfig) ;

case 'push_ei'
	y = evalin('base','data.mode.ei') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ei);
	zuicloseone(hfig) ;

case 'push_en'
	y = evalin('base','data.mode.en') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_en);
	zuicloseone(hfig) ;

case 'push_ej'
	y = evalin('base','data.mode.ej') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ej);
	zuicloseone(hfig) ;

case 'push_ve'
	y = evalin('base','data.mode.ve') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ve);
	zuicloseone(hfig) ;

case 'push_ep'
	y = evalin('base','data.mode.ep') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ep);
	zuicloseone(hfig) ;

case 'push_ie'
	y = evalin('base','data.mode.ie') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ie);
	zuicloseone(hfig) ;

case 'push_ii'
	y = evalin('base','data.mode.ii') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ii);
	zuicloseone(hfig) ;

case 'push_in'
	y = evalin('base','data.mode.in') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_in);
	zuicloseone(hfig) ;

case 'push_ij'
	y = evalin('base','data.mode.ij') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ij);
	zuicloseone(hfig) ;

case 'push_vi'
	y = evalin('base','data.mode.vi') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_vi);
	zuicloseone(hfig) ;

case 'push_ip'
	y = evalin('base','data.mode.ip') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ip);
	zuicloseone(hfig) ;

case 'push_ne'
	y = evalin('base','data.mode.ne') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ne);
	zuicloseone(hfig) ;

case 'push_ni'
	y = evalin('base','data.mode.ni') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ni);
	zuicloseone(hfig) ;

case 'push_nn'
	y = evalin('base','data.mode.nn') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_nn);
	zuicloseone(hfig) ;

case 'push_nj'
	y = evalin('base','data.mode.nj') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_nj);
	zuicloseone(hfig) ;

case 'push_vn'
	y = evalin('base','data.mode.vn') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_vn);
	zuicloseone(hfig) ;

case 'push_fefe'
	y = evalin('base','data.mode.fefe') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_fefe);
	zuicloseone(hfig) ;

case 'push_fev'
	y = evalin('base','data.mode.fev') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_fev);
	zuicloseone(hfig) ;

case 'push_fifi'
	y = evalin('base','data.mode.fifi') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_fifi);
	zuicloseone(hfig) ;

case 'push_fiv'
	y = evalin('base','data.mode.fiv') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_fiv);
	zuicloseone(hfig) ;

case 'push_rotc'
	y = evalin('base','data.mode.rotc') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_rotc);
	zuicloseone(hfig) ;

case 'push_rotv'
	y = evalin('base','data.mode.rotv') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_rotv);
	zuicloseone(hfig) ;
case 'push_eta'
	y = evalin('base','data.mode.eta') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_eta);
	zuicloseone(hfig) ;

case 'push_jboot'
	y = evalin('base','data.mode.jboot') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_jboot);
	zuicloseone(hfig) ;

case 'push_plot'
	y = evalin('base','data.mode.plot') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_plot);
	zuicloseone(hfig) ;

case 'push_zeff'
	y = evalin('base','data.mode.zeff') ;
	setappdata(hfig_pere,'y_import',y) ;
	zuireset(h.push_zeff);
	zuicloseone(hfig) ;

case 'push_ae'
	y = evalin('base','data.mode.ae') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_ae);
	zuicloseone(hfig) ;

case 'push_asser'
	y = evalin('base','data.mode.asser') ;
	setappdata(hfig_pere,'y_import',y) 
	zuireset(h.push_asser);
	zuicloseone(hfig) ;

case {'close'}
	zuicloseone(hfig);	
	
case {'init'}
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
	zuireset(h.raz);
	
case {'btn_annul','close'}
	y = [] ;
	setappdata(hfig_pere,'y_import',y) ;
	
	zuicloseone(hfig) ;
	
otherwise
	warning('ation non prise en compte')
	   
end

