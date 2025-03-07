% ZUIEDIT_MODE_IMPORT importe valeurs du mode de calcul choisi vers le formulaire d'edition de mode
%--------------------------------------------------------------
% fichier zuiedit_mode_import.m ->  
%	zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	choix d'un mode de calcul pour importer ses valeurs
%	vers le formulaire d'edition de mode
%
% syntaxe  :
%	hout = zuiedit_mode_import(hfig_pere) ;
%
% entree :
%	hfig_pere = handle du formulaire 'appelant'
%
% sortie
%	hout      = handle du formulaire cree
%
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 18/09/2003.
% 
% liste des modifications : 
%
%   * 14/03/2002 -> ajout de la stabilite MHD
%   * 10/12/2002 -> interface en anglais
%   * 18/09/2003 -> ajout mode.rotc
%
%--------------------------------------------------------------
function hout = zuiedit_mode_import(hfig_pere)

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle('mode_import') ;
if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

% on recupere les info
info = zinfo ;

% formulaire 
form={};

% Titre
sepa = {'separation_comm','frame','',3,''};
form{1} = {sepa};

colj = {'jump','jump','',[],''};
col1 = {'nom_mode','text@full','calcul mode',[],''};
form{length(form)+1} = {col1};

sepa = {'separation_comm','frame','',3,''};
form{1} = {sepa};

% Séparation
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

col1  = {'push_impur'  ,'push','impur'  ,0,info.data.mode.impur,[],''} ;
col2  = {'push_psi'    ,'push','psi'    ,0,info.data.mode.psi,[],''} ;
col3  = {'push_nel'    ,'push','nel'    ,0,info.data.mode.nel,[],''} ;
col4  = {'push_pe'     ,'push','pe'     ,0,info.data.mode.pe,[],''} ;
col5  = {'push_pion'   ,'push','pion'   ,0,info.data.mode.pion,[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,colj,colj};
col6  = {'push_equi'   ,'push','equi'   ,0,info.data.mode.equi,[],''} ;
col7  = {'push_neo'    ,'push','neo'    ,0,info.data.mode.neo,[],''} ;
col8  = {'push_fluce'  ,'push','fluce'  ,0,info.data.mode.fluce,[],''} ;
col9  = {'push_flucion','push','flucion',0,info.data.mode.flucion,[],''} ;
col10 = {'push_rot'    ,'push','rot'    ,0,info.data.mode.rot,[],''} ;
form{length(form)+1} = {col6,col7,col8,col9,col10,colj,colj};

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''} ;
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
col7 = {'sources','text@full','boundary conditions',1,'',''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

col1 = {'push_cons_psi'    ,'push','psi'    ,0,info.data.mode.cons.psi,[],''} ;
col2 = {'push_cons_ne'     ,'push','ne'     ,0,info.data.mode.cons.ne,[],''} ;
col3 = {'push_cons_pe'     ,'push','pe'     ,0,info.data.mode.cons.pe,[],''} ;
col4 = {'push_cons_pion'   ,'push','pion'   ,0,info.data.mode.cons.pion,[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,colj,colj,colj};
col5 = {'push_cons_fluce'  ,'push','fluce'  ,0,info.data.mode.cons.fluce,[],''} ;
col6 = {'push_cons_flucion','push','flucion',0,info.data.mode.cons.flucion,[],''} ;
col7 = {'push_cons_rot'    ,'push','rot'    ,0,info.data.mode.cons.rot,[],''} ;
%col8 = {'push_cons_zeffm'  ,'push','zeffm'  ,0,info.data.mode.cons.zeffm,[],''} ;
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
col7 = {'sources','text@full','bulk-edge link',1,'',''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

col1 = {'push_bord_ne'    ,'push','ne'    ,0,info.data.mode.consbord.ne,[],''} ;
col2 = {'push_bord_te'    ,'push','te'    ,0,info.data.mode.consbord.te,[],''} ;
col3 = {'push_bord_ti'    ,'push','ti'    ,0,info.data.mode.consbord.ti,[],''} ;
col4 = {'push_bord_fluxge','push','fluxge',0,info.data.mode.consbord.fluxge,[],''} ;
col5 = {'push_bord_fluxqe','push','fluxqe',0,info.data.mode.consbord.fluxqe,[],''} ;
col6 = {'push_bord_fluxqi','push','fluxqi',0,info.data.mode.consbord.fluxqi,[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,colj};

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'void','frame','void',3,''} ;
col7 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% mhd
col1 = {'void','jump','void',[],''} ;
col2 = {'void','jump','void',[],''} ;
col3 = {'void','jump','void',[],''} ;
col4 = {'void','jump','void',[],''} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
col7 = {'mhd','text@full','mhd',1,'',''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

col1 = {'push_mhd_dds'   ,'push','sawteeth'   ,0,info.data.mode.mhd.dds,[],''} ;
col2 = {'push_mhd_elm'   ,'push','elm'   ,0,info.data.mode.mhd.elm,[],''} ;
col3 = {'push_mhd_limite','push','boundary',0,info.data.mode.mhd.limite,[],''} ;
col4 = {'push_mhd_stab'  ,'push','stability',0,info.data.mode.mhd.stab,[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,colj,colj,colj};

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'void','frame','void',3,''} ;
col7 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% sources
col1 = {'void','jump','void',[],''} ;
col2 = {'void','jump','void',[],''} ;
col3 = {'void','jump','void',[],''} ;
col4 = {'void','jump','void',[],''} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
col7 = {'sources','text@full','sources',1,'',''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

col1  = {'push_fci'   ,'push','fci'   ,0,info.data.mode.fci,[],''} ;
col2  = {'push_fce'   ,'push','fce'   ,0,info.data.mode.fce,[],''} ;
col3  = {'push_hyb'   ,'push','lh'   ,0,info.data.mode.hyb,[],''} ;
col4  = {'push_idn'   ,'push','nbi'   ,0,info.data.mode.idn,[],''} ;
col5  = {'push_n0'    ,'push','n0'    ,0,info.data.mode.n0,[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,colj,colj};
col6  = {'push_bord'  ,'push','edge'  ,0,info.data.mode.bord,[],''} ;
col7  = {'push_glacon','push','pellet',0,info.data.mode.glacon,[],''} ;
col8  = {'push_fus'   ,'push','fus'   ,0,info.data.mode.fus,[],''} ;
col9  = {'push_ohm'   ,'push','ohm'   ,0,info.data.mode.ohm,[],''} ;
col10 = {'push_qneo'  ,'push','qneo'  ,0,info.data.mode.qneo,[],''} ;
form{length(form)+1} = {col6,col7,col8,col9,col10,colj,colj};
col11 = {'push_qei'   ,'push','qei'   ,0,info.data.mode.qei,[],''} ;
col12 = {'push_prad'  ,'push','prad'  ,0,info.data.mode.prad,[],''} ;
col13 = {'push_brem'  ,'push','brem'  ,0,info.data.mode.brem,[],''} ;
col14 = {'push_cyclo' ,'push','cyclo' ,0,info.data.mode.cyclo,[],''} ;
form{length(form)+1} = {col11,col12,col13,col14,colj,colj,colj};

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'void','frame','void',3,''} ;
col7 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% coefficient de transport
col1 = {'void','jump','void',[],''} ;
col2 = {'void','jump','void',[],''} ;
col3 = {'void','jump','void',[],''} ;
col4 = {'void','jump','void',[],''} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
col7 = {'sources','text@full','transport coefficients',1,'',''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

col1  = {'push_ee','push','ee',0,info.data.mode.ee,[],''} ;
col2  = {'push_ei','push','ei',0,info.data.mode.ei,[],''} ;
col3  = {'push_en','push','en',0,info.data.mode.en,[],''} ;
col4  = {'push_ej','push','ej',0,info.data.mode.ej,[],''} ;
col5  = {'push_ve','push','ve',0,info.data.mode.ve,[],''} ;
col6  = {'push_ep','push','ep',0,info.data.mode.ep,[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,colj};
col7  = {'push_ie','push','ie',0,info.data.mode.ie,[],''} ;
col8  = {'push_ii','push','ii',0,info.data.mode.ii,[],''} ;
col9  = {'push_in','push','in',0,info.data.mode.in,[],''} ;
col10 = {'push_ij','push','ij',0,info.data.mode.ij,[],''} ;
col11 = {'push_vi','push','vi',0,info.data.mode.vi,[],''} ;
col12 = {'push_ip','push','ip',0,info.data.mode.ip,[],''} ;
form{length(form)+1} = {col7,col8,col9,col10,col11,col12,colj};
col13 = {'push_ne','push','ne',0,info.data.mode.ne,[],''} ;
col14 = {'push_ni','push','ni',0,info.data.mode.ni,[],''} ;
col15 = {'push_nn','push','nn',0,info.data.mode.nn,[],''} ;
col16 = {'push_nj','push','nj',0,info.data.mode.nj,[],''} ;
col17 = {'push_vn','push','vn',0,info.data.mode.vn,[],''} ;
form{length(form)+1} = {col13,col14,col15,col16,col17,colj,colj};
col18 = {'push_fefe' ,'push','fefe' ,0,info.data.mode.fefe,[],''} ;
col19 = {'push_fev'  ,'push','fev'  ,0,info.data.mode.fev,[],''} ;
col20 = {'push_fifi' ,'push','fifi' ,0,info.data.mode.fifi,[],''} ;
col21 = {'push_fiv'  ,'push','fiv'  ,0,info.data.mode.fiv,[],''} ;
col22 = {'push_rotc' ,'push','rot' ,0,info.data.mode.rotc,[],''} ;
col23 = {'push_rotv' ,'push','rotv' ,0,info.data.mode.rotv,[],''} ;
form{length(form)+1} = {col18,col19,col20,col21,col22,col23,colj};
col23 = {'push_eta'  ,'push','eta'  ,0,info.data.mode.eta,[],''} ;
col24 = {'push_jboot','push','jboot',0,info.data.mode.jboot,[],''} ;
form{length(form)+1} = {col23,col24,colj,colj,colj,colj,colj};

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'void','frame','void',3,''} ;
col7 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% autres
col1 = {'void','jump','void',[],''} ;
col2 = {'void','jump','void',[],''} ;
col3 = {'void','jump','void',[],''} ;
col4 = {'void','jump','void',[],''} ;
col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;
col7 = {'sources','text@full','autres',1,'',''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

col1  = {'push_plot'        ,'push','plot'   ,0,info.data.mode.plot,[],''} ;
col2  = {'push_zeff'        ,'push','zeff'   ,0,info.data.mode.zeff,[],''} ;
col3  = {'push_ae'          ,'push','ae'     ,0,info.data.mode.ae,[],''} ;
col4  = {'push_asser'       ,'push','asser'  ,0,info.data.mode.asser,[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,colj,colj};

col1 = {'void','frame','void',3,''} ;
col2 = {'void','frame','void',3,''} ;
col3 = {'void','frame','void',3,''} ;
col4 = {'void','frame','void',3,''};
col5 = {'void','frame','void',3,''} ;
col6 = {'void','frame','void',3,''} ;
col7 = {'separation_comm','frame@full','',3,''};
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

% bouton Annulation
btn{1} = {'btn_annul','push','cancelling',0,'to leave the formular'} ;

% Formulaire
hout=zuicreeform(' ','mode_import','zuiedit_mode_import_fct','',form,btn) ;

[hform,hui] = zuiformhandle('mode_import') ;
setappdata(hform,'pere',hfig_pere) ;
