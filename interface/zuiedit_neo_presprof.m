function hout=zuiedit_iondens_neo

% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle('ed_neo_presprof') ;
if ishandle(hout)
        zuiformvisible(hout) ;
	return
end

% infos pour les tooltips
info = zinfo ;

% formulaire avec la sous structure from
form={};

% Titre
% -----
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text@full','profiles',[],''};
form{length(form)+1} = {col1};

% S�aration
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% data.source.ohm
% ---------------
col1 = {'ohm'       ,'text' ,'ohm',10,info.data.source.ohm,''} ;
col2 = {'edit_ohm'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_ohm','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ohm','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4 };

% data.source.qei
% ---------------
col1 = {'qei'       ,'text' ,'qei',10,info.data.source.qei,''} ;
col2 = {'edit_qei'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_qei','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_qei','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.jboot
% -----------------
col1 = {'jboot'       ,'text' ,'jboot',10,info.data.source.jboot,''} ;
col2 = {'edit_jboot'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_jboot','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_jboot','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% separation
% ----------
sepa = {'separation_comm','frame','',1,''};
form{length(form)+1} = {sepa};
colj = {'jump','jump',' ',[],''};
colj = {'jump','jump',' ',[],''};
form{length(form)+1} = {colj,colj};

% Bouton Quit
% ----------
comm{1}={'btn_quit','radio@center','Close',0,''};

% Formulaire
% ----------
hout=zuicreeform('Edition','ed_neo_presprof','zuiedit_neo_presprof_fct','',form,comm) ;


