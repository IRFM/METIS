% ZUIDEDIT_CONFIG_MODE_FCT   gestion des callbacksdu formulaire de configuration de mode
%--------------------------------------------------------------
% fichier zuiedit_config_mode_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	de configuration de mode
%
% syntaxe :
%	zuiedit_config_mode_fct(action)
%
% entrees :
%  action       =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 3.0, du 10/01/2005.
% 
% liste des modifications : 
% * 27/10/2001 -> changement de certaines legendes
% * 14/03/2002 -> ajout de la stabilite MHD
% * 19/03/2002 -> correction bug mode mhd
% * 10/12/2002 -> interface anglais
% * 03/09/2003 -> ajout gestion du mode ripple
% * 18/09/2003 -> bug mode.rotc
% * 23/03/2004 -> bug data.mode.cons.rot
% * 10/01/2005 -> modification des mode de cyclo
%
%--------------------------------------------------------------
function zuiedit_config_mode_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% disp('callback de zuiedit_config_mode_fct: ')
% disp(action)

% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('config_mode') ;

% information pour l'assistant
zuicr(hfig,action) ;

% selon ation
switch lower(action)

case 'radio_impur'
	zuiedit_mode('data.mode.impur',{1,2},{'read from data','calculated'}) ;
	zuireset(h.radio_impur);

case 'radio_psi'
	zuiedit_mode('data.mode.psi',{1,2},{'read from data','calculated'}) ;
	zuireset(h.radio_psi);

case 'radio_nel'
	zuiedit_mode('data.mode.nel',{1,2,3,4,5},{'read from data','calculated','total density controled','toltal density controled + shape ngr/nbar','toltal density controled + shape Wiesen'}) ;
	zuireset(h.radio_nel);

case 'radio_pe'
	zuiedit_mode('data.mode.pe',{1,2},{'read from data','calculated'}) ;
	zuireset(h.radio_pe);

case 'radio_pion'
	zuiedit_mode('data.mode.pion',{1,2},{'read from data','calculated'}) ;
	zuireset(h.radio_pion);

case 'radio_equi'
	zuiedit_mode('data.mode.equi',{1,2,3},{'read from data','calculated','k -> k+1'}) ;
	zuireset(h.radio_equi);

case 'radio_neo'
	zuiedit_mode('data.mode.neo',{1,2},{'read from data','calculated'}) ;
	zuireset(h.radio_neo);

case 'radio_rip'
	zuiedit_mode('data.mode.rip',{1,2},{'read from data','calculated'}) ;
	zuireset(h.radio_rip);

case 'radio_fluce'
	zuiedit_mode('data.mode.fluce',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_fluce);

case 'radio_flucion'
	zuiedit_mode('data.mode.flucion',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_flucion);

case 'radio_rot'
	zuiedit_mode('data.mode.rot',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_rot);

case 'radio_cons_psi'
	zuiedit_mode('data.mode.cons.psi',{0,1,2,3},{'Ip','Vloop','Psi(1)','Psi(1) free boundary'}) ;
	zuireset(h.radio_cons_psi);

case 'radio_cons_ne'
	zuiedit_mode('data.mode.cons.ne',{0,1},{'Ne(1)','ge(1)'}) ;
	zuireset(h.radio_cons_ne);

case 'radio_cons_pe'
	zuiedit_mode('data.mode.cons.pe',{0,1,2},{'Te(1)','qe(1)','Pe(1)'}) ;
	zuireset(h.radio_cons_pe);

case 'radio_cons_pion'
	zuiedit_mode('data.mode.cons.pion',{0,1,2},{'Ti(1)','qi(1)','Pion(1)'}) ;
	zuireset(h.radio_cons_pion);

case 'radio_cons_fluce'
	zuiedit_mode('data.mode.cons.fluce',{0,1},{'set to zero','read from data'}) ;
	zuireset(h.radio_cons_fluce);

case 'radio_cons_flucion'
	zuiedit_mode('data.mode.cons.flucion',{0,1},{'set to zero','read from data'}) ;
	zuireset(h.radio_cons_flucion);

case 'radio_cons_rot'
	zuiedit_mode('data.mode.cons.rot',{0,1},{'Rot(1)','edge flux of Rot'}) ;
	zuireset(h.radio_cons_rot);

case 'radio_cons_zeffm'
	zuiedit_mode('data.mode.cons.zeffm',{0,1},{'read from data','0d scaling'}) ;
	zuireset(h.radio_cons_zeffm);

case 'radio_bord_ne'
	zuiedit_mode('data.mode.consbord.ne',{0,1},{'read from data','edge function'}) ;
	zuireset(h.radio_bord_ne);

case 'radio_bord_te'
	zuiedit_mode('data.mode.consbord.te',{0,1},{'read from data','edge function'}) ;
	zuireset(h.radio_bord_te);

case 'radio_bord_ti'
	zuiedit_mode('data.mode.consbord.ti',{0,1},{'read from data','edge function'}) ;
	zuireset(h.radio_bord_ti);

case 'radio_bord_fluxge'
	zuiedit_mode('data.mode.consbord.fluxge',{0,1},{'read from data','edge function'}) ;
	zuireset(h.radio_bord_fluxge);

case 'radio_bord_fluxqe'
	zuiedit_mode('data.mode.consbord.fluxqe',{0,1},{'read from data','edge function'}) ;
	zuireset(h.radio_bord_fluxqe);

case 'radio_bord_fluxqi'
	zuiedit_mode('data.mode.consbord.fluxqi',{0,1},{'read from data','edge function'}) ;
	zuireset(h.radio_bord_fluxqi);

case 'radio_mhd_dds'
	zuiedit_mode('data.mode.mhd.dds',{0,1},{'set to zero','calculated mode'}) ;
	zuireset(h.radio_mhd_dds);

case 'radio_mhd_elm'
	zuiedit_mode('data.mode.mhd.elm',{0,1},{'set to zero','calculated mode'}) ;
	zuireset(h.radio_mhd_elm);

case 'radio_mhd_limite'
	zuiedit_mode('data.mode.mhd.limite',{0,1},{'set to zero','calculated mode'}) ;
	zuireset(h.radio_mhd_limite);

case 'radio_mhd_stab'
	zuiedit_mode('data.mode.mhd.stab',{0,1,2,3,4},{'set to zero','read from data',...
	             'calculated','k->k+1','post-processing mode'}) ;
	zuireset(h.radio_mhd_stab);

case 'radio_fci'
	zuiedit_mode('data.mode.fci',{0,1,2,3}, ...
                    {'set to zero','read from data','calculated','k->k+1'}) ;
	zuireset(h.radio_fci);

case 'radio_fce'
	zuiedit_mode('data.mode.fce',{0,1,2,3}, ...
                    {'set to zero','read from data','calculated ','k->k+1'}) ;
	zuireset(h.radio_fce);

case 'radio_hyb'
	zuiedit_mode('data.mode.hyb',{0,1,2,3}, ...
                    {'set to zero','read from data','calculated','k->k+1'}) ;
	zuireset(h.radio_hyb);

case 'radio_idn'
	zuiedit_mode('data.mode.idn',{0,1,2,3}, ...
                    {'set to zero','read from data','calculated','k->k+1'}) ;
	zuireset(h.radio_idn);

case 'radio_n0'
	zuiedit_mode('data.mode.n0',{0,1,2,3}, ...
                    {'set to zero','read from data','calculated','k->k+1'}) ;
	zuireset(h.radio_n0);

case 'radio_bord'
	zuiedit_mode('data.mode.bord',{0,1,2,3}, ...
                    {'set to zero','read from data','calculated','k->k+1'}) ;
	zuireset(h.radio_bord);

case 'radio_glacon'
	zuiedit_mode('data.mode.glacon',{0,1},{'set to zero','read from data'}) ;
	zuireset(h.radio_glacon);

case 'radio_fus'
	zuiedit_mode('data.mode.fus',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_fus);

case 'radio_ohm'
	zuiedit_mode('data.mode.ohm',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_ohm);

case 'radio_qneo'
	zuiedit_mode('data.mode.qneo',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_qneo);

case 'radio_qei'
	zuiedit_mode('data.mode.qei',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_qei);

case 'radio_prad'
	zuiedit_mode('data.mode.prad',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_prad);

case 'radio_brem'
	zuiedit_mode('data.mode.brem',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_brem);

case 'radio_cyclo'
	zuiedit_mode('data.mode.cyclo',{0,1,2,3},{'set to zero','read from data','calculated','k -> k+1'}) ;
	zuireset(h.radio_cyclo);

case 'radio_ee'
	zuiedit_mode('data.mode.ee',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_ee);

case 'radio_ei'
	zuiedit_mode('data.mode.ei',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_ei);

case 'radio_en'
	zuiedit_mode('data.mode.en',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_en);

case 'radio_ej'
	zuiedit_mode('data.mode.ej',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_ej);

case 'radio_ve'
	zuiedit_mode('data.mode.ve',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_ve);

case 'radio_ep'
	zuiedit_mode('data.mode.ep',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_ep);

case 'radio_ie'
	zuiedit_mode('data.mode.ie',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_ie);

case 'radio_ii'
	zuiedit_mode('data.mode.ii',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_ii);

case 'radio_in'
	zuiedit_mode('data.mode.in',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_in);

case 'radio_ij'
	zuiedit_mode('data.mode.ij',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_ij);

case 'radio_vi'
	zuiedit_mode('data.mode.vi',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_vi);

case 'radio_ip'
	zuiedit_mode('data.mode.ip',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_ip);

case 'radio_ne'
	zuiedit_mode('data.mode.ne',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_ne);

case 'radio_ni'
	zuiedit_mode('data.mode.ni',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_ni);

case 'radio_nn'
	zuiedit_mode('data.mode.nn',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_nn);

case 'radio_nj'
	zuiedit_mode('data.mode.nj',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_nj);

case 'radio_vn'
	zuiedit_mode('data.mode.vn',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_vn);

case 'radio_fefe'
	zuiedit_mode('data.mode.fefe',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_fefe);

case 'radio_fev'
	zuiedit_mode('data.mode.fev',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_fev);

case 'radio_fifi'
	zuiedit_mode('data.mode.fifi',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_fifi);

case 'radio_fiv'
	zuiedit_mode('data.mode.fiv',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_fiv);

case 'radio_rotc'
	zuiedit_mode('data.mode.rotc',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_rotc);

case 'radio_rotv'
	zuiedit_mode('data.mode.rotv',{0,1,2},{'set to zero','read from data','calculated'}) ;
	zuireset(h.radio_rotv);

case 'radio_eta'
	zuiedit_mode('data.mode.eta',{1,2},{'read from data','calculated'}) ;
	zuireset(h.radio_eta);

case 'radio_jboot'
	zuiedit_mode('data.mode.jboot',{1,2},{'read from data','calculated'}) ;
	zuireset(h.radio_jboot);

case 'radio_plot'
	zuiedit_mode('data.mode.plot',{0,1},{'off','on'}) ;
	zuireset(h.radio_plot);

case 'radio_zeff'
	zuiedit_mode('data.mode.zeff',{1,2},{'read from data','calculated'}) ;
	zuireset(h.radio_zeff);

case 'radio_ae'
	zuiedit_mode('data.mode.ae',{1,2},{'read from data','calculated'}) ;
	zuireset(h.radio_ae);

case 'radio_asser'
	zuiedit_mode('data.mode.asser',{1,2},{'read from data','calculated'}) ;
	zuireset(h.radio_asser);
	
case {'init','raz'}
	zuiformvisible(hfig);
	zuiformreset(hfig);
	zuiuploadform(hfig);
	
	zuireset(h.raz);
	
case {'btn_quit','close'}
	zuireset(h.btn_quit);
	zuiformcache(hfig) ;
        zuicloseone(hfig);
	
otherwise
	warning('action not taking into account')
	   
end

