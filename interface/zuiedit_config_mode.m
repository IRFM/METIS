% ZUIEDIT_MODE   formulaire configuration de mode de calcul
%--------------------------------------------------------------
% fichier zuiedit_mode.m ->  
%		zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	creation du formulaire de configuration de mode de calcul
%
% syntaxe  :
%	hout=zuiedit_config_mode ;
%
% entree :
%
% sorties :
%	hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 18/09/2003.
% 
% liste des modifications : 
% * 14/03/2002 -> ajout de la stabilite MHD
% * 03/09/2003 -> ajout du ripple
% * 18/09/2003 -> bug sur mode.rotc
%
%--------------------------------------------------------------
function hform=zuiedit_config_mode

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle('config_mode');
if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

% on recupere les info
info = zinfo ;

% formulaire d'edition de modes
% avec la sous structure from
form={};

% Titre
sepa = {'separation_comm','frame','',3,''};
form{1} = {sepa};

%colj = {'jump','jump',' ',[],''};
colj = {'sources','text','',[],'',''} ;
col1 = {'nom_mode','text@full','computation modes',[],''};
form{length(form)+1} = {colj,col1};

sepa = {'separation_comm','frame','',3,''};
form{1} = {sepa};

% S�aration
col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''} ;
col5 = {'void','frame','void',3,''} ;
col6 = {'void','frame','void',3,''} ;
col7 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% mode pde
col1 = {'void','jump','void',[],''} ;
col2 = {'void','jump','void',[],''} ;
col3 = {'void','jump','void',[],''} ;
col4 = {'void','jump','void',[],''} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
col7 = {'pde','text@full','P.D.E',1,'',''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

%col1  = {'radio_impur'  ,'radio','impur'  ,0,info.data.mode.impur,[],''} ;
col3  = {'radio_psi'    ,'radio','psi'    ,0,info.data.mode.psi,[],''} ;
col4  = {'radio_nel'    ,'radio','nel'    ,0,info.data.mode.nel,[],''} ;
col5  = {'radio_pe'     ,'radio','pe'     ,0,info.data.mode.pe,[],''} ;
form{length(form)+1} = {col3,col4,col5,colj,colj,colj,colj};
%col1  = {'radio_equi'   ,'radio','equi'   ,0,info.data.mode.equi,[],''} ;
%col2  = {'radio_neo'    ,'radio','neo'    ,0,info.data.mode.neo,[],''} ;
%col3 = {'radio_rip'    ,'radio','ripple'    ,0,info.data.mode.rot,[],''} ;
%col4 = {'radio_ae'          ,'radio','ae'     ,0,info.data.mode.ae,[],''} ;
col3  = {'radio_pion'   ,'radio','pion'   ,0,info.data.mode.pion,[],''} ;
col4  = {'radio_rot'    ,'radio','rot'    ,0,info.data.mode.rot,[],''} ;
col5  = {'radio_fluce'  ,'radio','fluce'  ,0,info.data.mode.fluce,[],''} ;
col6  = {'radio_flucion','radio','flucion',0,info.data.mode.flucion,[],''} ;
form{length(form)+1} = {col3,col4,col5,col6,colj,colj,colj};

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'void','frame','void',3,''} ;
col7 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% conditions aux limites
col1 = {'void','jump','void',[],''} ;
col2 = {'void','jump','void',[],''} ;
col3 = {'void','jump','void',[],''} ;
col4 = {'void','jump','void',[],''} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
col7 = {'sources','text@full','boundary condition model',1,'',''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

col1 = {'radio_cons_psi'    ,'radio','psi'    ,0,info.data.mode.cons.psi,[],''} ;
col2 = {'radio_cons_ne'     ,'radio','ne'     ,0,info.data.mode.cons.ne,[],''} ;
col3 = {'radio_cons_pe'     ,'radio','pe'     ,0,info.data.mode.cons.pe,[],''} ;
col4 = {'radio_cons_pion'   ,'radio','pion'   ,0,info.data.mode.cons.pion,[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,colj,colj,colj};
col5 = {'radio_cons_fluce'  ,'radio','fluce'  ,0,info.data.mode.cons.fluce,[],''} ;
col6 = {'radio_cons_flucion','radio','flucion',0,info.data.mode.cons.flucion,[],''} ;
col7 = {'radio_cons_rot'    ,'radio','rot'    ,0,info.data.mode.cons.rot,[],''} ;
%col8 = {'radio_cons_zeffm'  ,'radio','zeffm'  ,0,info.data.mode.cons.zeffm,[],''} ;
form{length(form)+1} = {col5,col6,col7,colj,colj,colj,colj};

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'void','frame','void',3,''} ;
col7 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% bord
col1 = {'void','jump','void',[],''} ;
col2 = {'void','jump','void',[],''} ;
col3 = {'void','jump','void',[],''} ;
col4 = {'void','jump','void',[],''} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
col7 = {'sources','text@full','edge calculation',1,'',''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

col1 = {'radio_bord_ne'    ,'radio','ne'    ,0,info.data.mode.consbord.ne,[],''} ;
col2 = {'radio_bord_te'    ,'radio','te'    ,0,info.data.mode.consbord.te,[],''} ;
col3 = {'radio_bord_ti'    ,'radio','ti'    ,0,info.data.mode.consbord.ti,[],''} ;
col4 = {'radio_bord_fluxge','radio','fluxge',0,info.data.mode.consbord.fluxge,[],''} ;
col5 = {'radio_bord_fluxqe','radio','fluxqe',0,info.data.mode.consbord.fluxqe,[],''} ;
col6 = {'radio_bord_fluxqi','radio','fluxqi',0,info.data.mode.consbord.fluxqi,[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,colj};

%col1 = {'void','frame','void',3,''} ;
%col2 = {'void','frame','void',3,''} ;
%col3 = {'void','frame','void',3,''} ;
%col4 = {'void','frame','void',3,''};
%col5 = {'void','frame','void',3,''} ;
%col6 = {'void','frame','void',3,''} ;
%col7 = {'separation_comm','frame@full','',3,''};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% mhd
%col1 = {'void','jump','void',[],''} ;
%col2 = {'void','jump','void',[],''} ;
%col3 = {'void','jump','void',[],''} ;
%col4 = {'void','jump','void',[],''} ;
%col5 = {'void','jump','void',[],''} ;
%col6 = {'void','jump','void',[],''} ;
%col7 = {'mhd','text@full','mhd',1,'',''};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

%col1 = {'radio_mhd_dds'   ,'radio','sawteeth'      ,0,info.data.mode.mhd.dds,[],''} ;
%col2 = {'radio_mhd_elm'   ,'radio','elm'      ,0,info.data.mode.mhd.elm,[],''} ;
%col3 = {'radio_mhd_limite','radio','limit'   ,0,info.data.mode.mhd.limite,[],''} ;
%col4 = {'radio_mhd_stab'  ,'radio','stability',0,info.data.mode.mhd.stab,[],''} ;
%form{length(form)+1} = {col1,col2,col3,col4,colj,colj,colj};

%col1 = {'void','frame','void',3,''} ;
%col2 = {'void','frame','void',3,''} ;
%col3 = {'void','frame','void',3,''} ;
%col4 = {'void','frame','void',3,''};
%col5 = {'void','frame','void',3,''} ;
%col6 = {'void','frame','void',3,''} ;
%col7 = {'separation_comm','frame@full','',3,''};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% sources
%col1 = {'void','jump','void',[],''} ;
%col2 = {'void','jump','void',[],''} ;
%col3 = {'void','jump','void',[],''} ;
%col4 = {'void','jump','void',[],''} ;
%col5 = {'void','jump','void',[],''} ;
%col6 = {'void','jump','void',[],''} ;
%col7 = {'sources','text@full','sources',1,'',''};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

%col1  = {'radio_fci'   ,'radio','icrh'   ,0,info.data.mode.fci,[],''} ;
%col2  = {'radio_fce'   ,'radio','ecrh'   ,0,info.data.mode.fce,[],''} ;
%col3  = {'radio_hyb'   ,'radio','lh'   ,0,info.data.mode.hyb,[],''} ;
%col4  = {'radio_idn'   ,'radio','nbi'   ,0,info.data.mode.idn,[],''} ;
%col5  = {'radio_n0'    ,'radio','n0'    ,0,info.data.mode.n0,[],''} ;
%form{length(form)+1} = {col1,col2,col3,col4,col5,colj,colj};
%col6  = {'radio_bord'  ,'radio','edge'  ,0,info.data.mode.bord,[],''} ;
%col7  = {'radio_glacon','radio','pellet',0,info.data.mode.glacon,[],''} ;
%col8  = {'radio_fus'   ,'radio','fus'   ,0,info.data.mode.fus,[],''} ;
%col9  = {'radio_ohm'   ,'radio','ohm'   ,0,info.data.mode.ohm,[],''} ;
%col10 = {'radio_qneo'  ,'radio','qneo'  ,0,info.data.mode.qneo,[],''} ;
%form{length(form)+1} = {col6,col7,col8,col9,col10,colj,colj};
%col11 = {'radio_qei'   ,'radio','qei'   ,0,info.data.mode.qei,[],''} ;
%col12 = {'radio_prad'  ,'radio','prad'  ,0,info.data.mode.prad,[],''} ;
%col13 = {'radio_brem'  ,'radio','brem'  ,0,info.data.mode.brem,[],''} ;
%col14 = {'radio_cyclo' ,'radio','cyclo' ,0,info.data.mode.cyclo,[],''} ;
%form{length(form)+1} = {col11,col12,col13,col14,colj,colj,colj};

%col1 = {'void','frame','void',3,''} ;
%col2 = {'void','frame','void',3,''} ;
%col3 = {'void','frame','void',3,''} ;
%col4 = {'void','frame','void',3,''};
%col5 = {'void','frame','void',3,''} ;
%col6 = {'void','frame','void',3,''} ;
%col7 = {'separation_comm','frame@full','',3,''};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% coefficient de transport
%col1 = {'void','jump','void',[],''} ;
%col2 = {'void','jump','void',[],''} ;
%col3 = {'void','jump','void',[],''} ;
%col4 = {'void','jump','void',[],''} ;
%col5 = {'void','jump','void',[],''} ;
%col6 = {'void','jump','void',[],''} ;
%col7 = {'sources','text@full','Transport Coefficients',1,'',''};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

%col1  = {'radio_ee','radio','ee',0,info.data.mode.ee,[],''} ;
%col2  = {'radio_ei','radio','ei',0,info.data.mode.ei,[],''} ;
%col3  = {'radio_en','radio','en',0,info.data.mode.en,[],''} ;
%col4  = {'radio_ej','radio','ej',0,info.data.mode.ej,[],''} ;
%col5  = {'radio_ve','radio','ve',0,info.data.mode.ve,[],''} ;
%col6  = {'radio_ep','radio','ep',0,info.data.mode.ep,[],''} ;
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,colj};
%col7  = {'radio_ie','radio','ie',0,info.data.mode.ie,[],''} ;
%col8  = {'radio_ii','radio','ii',0,info.data.mode.ii,[],''} ;
%col9  = {'radio_in','radio','in',0,info.data.mode.in,[],''} ;
%col10 = {'radio_ij','radio','ij',0,info.data.mode.ij,[],''} ;
%col11 = {'radio_vi','radio','vi',0,info.data.mode.vi,[],''} ;
%col12 = {'radio_ip','radio','ip',0,info.data.mode.ip,[],''} ;
%form{length(form)+1} = {col7,col8,col9,col10,col11,col12,colj};
%col13 = {'radio_ne','radio','ne',0,info.data.mode.ne,[],''} ;
%col14 = {'radio_ni','radio','ni',0,info.data.mode.ni,[],''} ;
%col15 = {'radio_nn','radio','nn',0,info.data.mode.nn,[],''} ;
%col16 = {'radio_nj','radio','nj',0,info.data.mode.nj,[],''} ;
%col17 = {'radio_vn','radio','vn',0,info.data.mode.vn,[],''} ;
%form{length(form)+1} = {col13,col14,col15,col16,col17,colj,colj};
%col18 = {'radio_fefe' ,'radio','fefe' ,0,info.data.mode.fefe,[],''} ;
%col19 = {'radio_fev'  ,'radio','fev'  ,0,info.data.mode.fev,[],''} ;
%col20 = {'radio_fifi' ,'radio','fifi' ,0,info.data.mode.fifi,[],''} ;
%col21 = {'radio_fiv'  ,'radio','fiv'  ,0,info.data.mode.fiv,[],''} ;
%col22 = {'radio_rotc' ,'radio','rot' ,0,info.data.mode.rotc,[],''} ;
%%col23 = {'radio_rotv' ,'radio','rotv' ,0,info.data.mode.rotv,[],''} ;
%form{length(form)+1} = {col18,col19,col20,col21,col22,col23,colj};
%col23 = {'radio_eta'  ,'radio','eta'  ,0,info.data.mode.eta,[],''} ;
%col24 = {'radio_jboot','radio','jboot',0,info.data.mode.jboot,[],''} ;
%form{length(form)+1} = {col23,col24,colj,colj,colj,colj,colj};

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'void','frame','void',3,''} ;
col7 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% autres
%col1 = {'void','jump','void',[],''} ;
%col2 = {'void','jump','void',[],''} ;
%col3 = {'void','jump','void',[],''} ;
%col4 = {'void','jump','void',[],''} ;
%col5 = {'void','jump','void',[],''} ;
%col6 = {'void','jump','void',[],''} ;
%col7 = {'sources','text@full','others',1,'',''};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

%col1  = {'radio_plot'        ,'radio','plot'   ,0,info.data.mode.plot,[],''} ;
%col2  = {'radio_zeff'        ,'radio','zeff'   ,0,info.data.mode.zeff,[],''} ;
%col3  = {'radio_cons_zeffm'  ,'radio','zeffm'  ,0,info.data.mode.cons.zeffm,[],''} ;
%col3  = {'radio_ae'          ,'radio','ae'     ,0,info.data.mode.ae,[],''} ;
%col4  = {'radio_asser'       ,'radio','asser'  ,0,info.data.mode.asser,[],''} ;
%form{length(form)+1} = {col1,col2,col3,col4,col5,colj,colj};

%col1 = {'void','frame','void',3,''} ;
%col2 = {'void','frame','void',3,''} ;
%col3 = {'void','frame','void',3,''} ;
%col4 = {'void','frame','void',3,''};
%col5 = {'void','frame','void',3,''} ;
%col6 = {'void','frame','void',3,''} ;
%col7 = {'separation_comm','frame@full','',3,''};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% bouton Quit
% -----------
comm{1}={'btn_quit','radio@center','Close',0,''};

% Formulaire
hout=zuicreeform(' ','config_mode','zuiedit_config_mode_fct','',form,comm) ;

