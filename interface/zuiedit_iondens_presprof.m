function hout=zuiedit_iondens_presprof

% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle('ed_iondens_presprof') ;
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

% data.prof.ae
% ------------
col1 = {'ae'        ,'text' ,'ae'    ,10,info.data.prof.ae,''} ;
col2 = {'edit_ae'  ,'radio','edit',0,'edit a profile',[],''} ;
col3 = {'param_ae' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_ae','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_ae','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.zeff
% --------------
col1 = {'zeff'       ,'text' ,'zeff'    ,10,info.data.prof.zeff,''} ;
col2 = {'edit_zeff'  ,'radio','edit',0,'edit a profile',[],''} ;
col3 = {'param_zeff' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_zeff','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_zeff','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.source.prad
% ----------------
col1 = {'prad'       ,'text' ,'prad',10,info.data.source.prad,''} ;
col2 = {'edit_prad'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'sources'    ,'text' ,''    ,10  ,''                ,''} ;
col4 = {'import_prad','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_prad','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.source.brem
% ----------------
col1 = {'brem'       ,'text' ,'brem',10,info.data.source.brem,''} ;
col2 = {'edit_brem'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'sources'    ,'text' ,''    ,10,'',''} ;
col4 = {'import_brem','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_brem','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.source.cyclo
% -----------------
col1 = {'cyclo'       ,'text' ,'cyclo',10,info.data.source.cyclo,''} ;
col2 = {'edit_cyclo'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'sources','text','',10,'',''} ;
col4 = {'import_cyclo','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_cyclo','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

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
hout=zuicreeform('Edition','ed_iondens_presprof','zuiedit_iondens_presprof_fct','',form,comm) ;


