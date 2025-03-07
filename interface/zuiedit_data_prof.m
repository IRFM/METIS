% ZUIEDIT_DATA_PROF formulaire d'�ition des profils
%---------------------------------------------------
% fichier zuiedit_data_prof.m ->  
%		zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	creation du formulaire d'editions des profils
%
% syntaxe  :
%	hout=zuiedit_data_prof;
%
% entree :
%
% sorties :
%	hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 2.0, du 10/12/2002.
% 
% liste des modifications : 
% 
%  * 30/08/2001 -> correction hout et hform
%  * 10/12/2002 -> interface en anglais
%
%--------------------------------------------------------------
function hout=zuiedit_data_prof

% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle('data_prof') ;
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

% data.prof.psi
% -------------
col1 = {'psi'       ,'text' ,'psi'    ,10,info.data.prof.psi,''} ;
col2 = {'edit_psi'  ,'radio','edit',0  ,'edit a profile',[],''} ;
col3 = {'param_psi' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_psi','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_psi','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.pe
% ------------
col1 = {'pe'       ,'text' ,'pe'    ,10,info.data.prof.pe,''} ;
col2 = {'edit_pe'  ,'radio','edit',0,'edit a profile',[],''} ;
col3 = {'param_pe' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_pe','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_pe','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.te
% ------------
col1 = {'te'       ,'text' ,'Te'    ,10,info.data.prof.te,''} ;
col2 = {'edit_te'  ,'radio','edit',0,'edit a profile',[],''} ;
col3 = {'param_te' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_te','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_te','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.pion
% --------------
col1 = {'pion'       ,'text' ,'pion'    ,10,info.data.prof.pion,''} ;
col2 = {'edit_pion'  ,'radio','edit',0,'edit a profile',[],''} ;
col3 = {'param_pion' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_pion','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_pion','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.ti
% --------------
col1 = {'ti'       ,'text' ,'Ti'    ,10,info.data.prof.ti,''} ;
col2 = {'edit_ti'  ,'radio','edit',0,'edit a profile',[],''} ;
col3 = {'param_ti' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_ti','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_ti','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.ne
% ------------
col1 = {'ne'       ,'text' ,'ne'    ,10,info.data.prof.ne,''} ;
col2 = {'edit_ne'  ,'radio','edit',0,'edit a profile',[],''} ;
col3 = {'param_ne' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_ne','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_ne','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.rot
% -------------
col1 = {'rot'       ,'text' ,'rot'    ,10,info.data.prof.rot,''} ;
col2 = {'edit_rot'  ,'radio','edit',0,'edit a profile',[],''} ;
col3 = {'param_rot' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_rot','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_rot','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.vtor_exp
% -------------
col1 = {'vtor_exp'       ,'text' ,'vtor_exp'    ,10,info.data.prof.vtor_exp,''} ;
col2 = {'edit_vtor_exp'  ,'radio','edit',0,'edit a profile',[],''} ;
col3 = {'param_vtor_exp' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_vtor_exp','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_vtor_exp','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.fluce
% ---------------
col1 = {'fluce'       ,'text' ,'fluce'    ,10,info.data.prof.fluce,''} ;
col2 = {'edit_fluce'  ,'radio','edit',0,'edit a profile',[],'','enable','off'} ;
col3 = {'param_fluce' ,'radio','fitting parameter',0,'fitted profile',[],'','enable','off'} ;
col4 = {'import_fluce','radio','load',0,'load the profile',[],'','enable','off'} ;
col5 = {'dessin_fluce','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.flucion
% -----------------
col1 = {'flucion'       ,'text' ,'flucion'    ,10,info.data.prof.flucion,''} ;
col2 = {'edit_flucion'  ,'radio','edit',0,'edit a profile',[],'','enable','off'} ;
col3 = {'param_flucion' ,'radio','fitting parameter',0,'fitted profile',[],'','enable','off'} ;
col4 = {'import_flucion','radio','load',0,'load the profile',[],'','enable','off'} ;
col5 = {'dessin_flucion','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.jmoy
% --------------
col1 = {'jmoy'       ,'text' ,'jmoy'    ,10,info.data.prof.jmoy,''} ;
col2 = {'edit_jmoy'  ,'radio','edit',0,'edit a profile',[],''} ;
col3 = {'param_jmoy' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_jmoy','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_jmoy','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.ae
% ------------
%col1 = {'ae'        ,'text' ,'ae'    ,10,info.data.prof.ae,''} ;
%col2 = {'edit_ae'  ,'radio','edit',0,'edit a profile',[],''} ;
%col3 = {'param_ae' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
%col4 = {'import_ae','radio','load',0,'load the profile',[],''} ;
%col5 = {'dessin_ae','radio','modify',0,'modify the profile',[],''} ;
%form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.zeff
% --------------
%col1 = {'zeff'       ,'text' ,'zeff'    ,10,info.data.prof.zeff,''} ;
%col2 = {'edit_zeff'  ,'radio','edit',0,'edit a profile',[],''} ;
%col3 = {'param_zeff' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
%col4 = {'import_zeff','radio','load',0,'load the profile',[],''} ;
%col5 = {'dessin_zeff','radio','modify',0,'modify the profile',[],''} ;
%form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.xdur
% --------------
col1 = {'xdur'       ,'text' ,'xdur'    ,10,info.data.prof.xdur,''} ;
col2 = {'edit_xdur'  ,'radio','edit',0,'edit a profile',[],''} ;
col3 = {'param_xdur' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
col4 = {'import_xdur','radio','load',0,'load the profile',[],''} ;
col5 = {'dessin_xdur','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.mhd_cd
% ----------------
%col1 = {'mhd_cd'       ,'text' ,'mhd_cd'    ,10,info.data.prof.mhd_cd,''} ;
%col2 = {'edit_mhd_cd'  ,'radio','edit',0,'edit a profile',[],''} ;
%col3 = {'param_mhd_cd' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
%col4 = {'import_mhd_cd','radio','load',0,'load the profile',[],''} ;
%col5 = {'dessin_mhd_cd','radio','modify',0,'modify the profile',[],''} ;
%form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.prof.ej
% ------------
%col1 = {'ej'       ,'text' ,'ej'    ,10,info.data.prof.ej,''} ;
%col2 = {'edit_ej'  ,'radio','edit',0,'edit a profile',[],''} ;
%col3 = {'param_ej' ,'radio','fitting parameter',0,'fitted profile',[],''} ;
%col4 = {'import_ej','radio','load',0,'load the profile',[],''} ;
%col5 = {'dessin_ej','radio','modify',0,'modify the profile',[],''} ;
%form{length(form)+1} = {col1,col2,col3,col4,col5} ;

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
hout=zuicreeform('Edition','data_prof','zuiedit_data_prof_fct','',form,comm) ;


