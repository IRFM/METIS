function hout=zuiedit_transport_coeff_mode

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle('ed_transport_coeff_mode');
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
col1 = {'nom_mode','text@full','Transport coefficient computation mode',[],''};
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

% coefficient de transport
%col1 = {'void','jump','void',[],''} ;
%col2 = {'void','jump','void',[],''} ;
%col3 = {'void','jump','void',[],''} ;
%col4 = {'void','jump','void',[],''} ;
%col5 = {'void','jump','void',[],''} ;
%col6 = {'void','jump','void',[],''} ;
%col7 = {'sources','text@full','Transport Coefficients',1,'',''};
%form{length(form)+1} = {col1,col2,col3,col4,col5,col6,col7} ;

col1  = {'radio_ee','radio','ee',0,info.data.mode.ee,[],''} ;
col2  = {'radio_ei','radio','ei',0,info.data.mode.ei,[],''} ;
col3  = {'radio_en','radio','en',0,info.data.mode.en,[],''} ;
col4  = {'radio_ej','radio','ej',0,info.data.mode.ej,[],''} ;
col5  = {'radio_ve','radio','ve',0,info.data.mode.ve,[],''} ;
col6  = {'radio_ep','radio','ep',0,info.data.mode.ep,[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5,col6,colj};
col7  = {'radio_ie','radio','ie',0,info.data.mode.ie,[],''} ;
col8  = {'radio_ii','radio','ii',0,info.data.mode.ii,[],''} ;
col9  = {'radio_in','radio','in',0,info.data.mode.in,[],''} ;
col10 = {'radio_ij','radio','ij',0,info.data.mode.ij,[],''} ;
col11 = {'radio_vi','radio','vi',0,info.data.mode.vi,[],''} ;
col12 = {'radio_ip','radio','ip',0,info.data.mode.ip,[],''} ;
form{length(form)+1} = {col7,col8,col9,col10,col11,col12,colj};
col13 = {'radio_ne','radio','ne',0,info.data.mode.ne,[],''} ;
col14 = {'radio_ni','radio','ni',0,info.data.mode.ni,[],''} ;
col15 = {'radio_nn','radio','nn',0,info.data.mode.nn,[],''} ;
col16 = {'radio_nj','radio','nj',0,info.data.mode.nj,[],''} ;
col17 = {'radio_vn','radio','vn',0,info.data.mode.vn,[],''} ;
form{length(form)+1} = {col13,col14,col15,col16,col17,colj,colj};
col18 = {'radio_fefe' ,'radio','fefe' ,0,info.data.mode.fefe,[],''} ;
col19 = {'radio_fev'  ,'radio','fev'  ,0,info.data.mode.fev,[],''} ;
col20 = {'radio_fifi' ,'radio','fifi' ,0,info.data.mode.fifi,[],''} ;
col21 = {'radio_fiv'  ,'radio','fiv'  ,0,info.data.mode.fiv,[],''} ;
col22 = {'radio_rotc' ,'radio','rot' ,0,info.data.mode.rotc,[],''} ;
col23 = {'radio_rotv' ,'radio','rotv' ,0,info.data.mode.rotv,[],''} ;
form{length(form)+1} = {col18,col19,col20,col21,col22,col23,colj};
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
hout=zuicreeform(' ','ed_transport_coeff_mode','zuiedit_transport_coeff_mode_fct','',form,comm) ;

