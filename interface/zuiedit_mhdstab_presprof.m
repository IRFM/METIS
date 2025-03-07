function hout=zuiedit_mhdstab_presprof

% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle('ed_mhdstab_presprof') ;
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

% data.prof.mhd_cd
% ----------------
col1 = {'mhd_cd'       ,'text' ,'mhd_cd'    ,10,info.data.prof.mhd_cd,''} ;
col2 = {'edit_mhd_cd'  ,'radio','edit',0,'edit a profile',[],''} ;
col3 = {'param_mhd_cd' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_mhd_cd','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_mhd_cd','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% Bouton Quit
% ----------
comm{1}={'btn_quit','radio@center','Close',0,''};

% Formulaire
% ----------
hout=zuicreeform('Edition','ed_mhdstab_presprof','zuiedit_mhdstab_presprof_fct','',form,comm) ;


