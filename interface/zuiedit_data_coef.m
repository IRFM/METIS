% ZUIEDIT_DATA_COEF formulaire d'�ition des profils des coefficients de transport
%---------------------------------------------------------------------------------
% fichier zuiedit_data_coef.m ->  
%	zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%
%
% syntaxe  :
%	hout=zuiedit_data_coef;
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
% * 11/12/2002 : interface en anglais
%--------------------------------------------------------------
function hform=zuiedit_data_coef

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle('data_coef') ;
if ishandle(hform)
        zuiformvisible(hform) ;
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

col1 = {'libel','text@full','transport coefficient profile',[],''};
form{length(form)+1} = {col1};

% S�aration
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% data.coef.eta
% -------------
col1 = {'eta'       ,'text' ,'eta'    ,10,info.data.coef.eta,''} ;
col2 = {'edit_eta'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_eta','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_eta','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.ee
% ------------
col1 = {'ee'       ,'text' ,'ee'    ,10,info.data.coef.ee,''} ;
col2 = {'edit_ee'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_ee','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ee','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.ei
% ------------
col1 = {'ei'       ,'text' ,'ei'    ,10,info.data.coef.ei,''} ;
col2 = {'edit_ei'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_ei','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ei','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.en
% ------------
col1 = {'en'       ,'text' ,'en'    ,10,info.data.coef.en,''} ;
col2 = {'edit_en'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_en','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_en','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.ej
% ------------
col1 = {'ej'       ,'text' ,'ej'    ,10,info.data.coef.ej,''} ;
col2 = {'edit_ej'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_ej','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ej','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.ve
% ------------
col1 = {'ve'       ,'text' ,'ve'    ,10,info.data.coef.ve,''} ;
col2 = {'edit_ve'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_ve','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ve','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.ep
% ------------
col1 = {'ep'       ,'text' ,'ep'    ,10,info.data.coef.ep,''} ;
col2 = {'edit_ep'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_ep','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ep','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.ie
% ------------
col1 = {'ie'        ,'text' ,'ie'    ,10,info.data.coef.ie,''} ;
col2 = {'edit_ie'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_ie','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ie','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.ii
% ------------
col1 = {'ii'        ,'text' ,'ii'    ,10,info.data.coef.ii,''} ;
col2 = {'edit_ii'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_ii','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ii','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.in
% ------------
col1 = {'in'       ,'text' ,'in'    ,10,info.data.coef.in,''} ;
col2 = {'edit_in'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_in','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_in','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.ij
% ------------
col1 = {'ij'       ,'text' ,'ij'    ,10,info.data.coef.ij,''} ;
col2 = {'edit_ij'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_ij','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ij','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.vi
% ------------
col1 = {'vi'       ,'text' ,'vi'    ,10,info.data.coef.vi,''} ;
col2 = {'edit_vi'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_vi','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_vi','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.ip
% ------------
col1 = {'ip'       ,'text' ,'ip'    ,10,info.data.coef.ip,''} ;
col2 = {'edit_ip'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_ip','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ip','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.ne
% ------------
col1 = {'ne'       ,'text' ,'ne'    ,10,info.data.coef.ne,''} ;
col2 = {'edit_ne'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_ne','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ne','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.ni
% ------------
col1 = {'ni'       ,'text' ,'ni'    ,10,info.data.coef.ni,''} ;
col2 = {'edit_ni'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_ni','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ni','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.nn
% ------------
col1 = {'nn'       ,'text' ,'nn'    ,10,info.data.coef.nn,''} ;
col2 = {'edit_nn'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_nn','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_nn','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.nj
% ------------
col1 = {'nj'       ,'text' ,'nj'    ,10,info.data.coef.nj,''} ;
col2 = {'edit_nj'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_nj','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_nj','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.vn
% ------------
col1 = {'vn'       ,'text' ,'vn'    ,10,info.data.coef.vn,''} ;
col2 = {'edit_vn'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_vn','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_vn','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.fev
% -------------
col1 = {'fev'       ,'text' ,'fev'    ,10,info.data.coef.fev,''} ;
col2 = {'edit_fev'  ,'radio','edit',0,'edit the profile',[],'','enable','off'} ;
col3 = {'import_ne','radio','load',0,'load the profile',[],'','enable','off'} ;
col4 = {'dessin_ne','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.fefe
% --------------
col1 = {'fefe'       ,'text' ,'fefe'    ,10,info.data.coef.fefe,''} ;
col2 = {'edit_fefe'  ,'radio','edit',0,'edit the profile',[],'','enable','off'} ;
col3 = {'import_fefe','radio','load',0,'load the profile',[],'','enable','off'} ;
col4 = {'dessin_fefe','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.fiv
% -------------
col1 = {'fiv'       ,'text' ,'fiv'    ,10,info.data.coef.fiv,''} ;
col2 = {'edit_fiv'  ,'radio','edit',0,'edit the profile',[],'','enable','off'} ;
col3 = {'import_fiv','radio','load',0,'load the profile',[],'','enable','off'} ;
col4 = {'dessin_fiv','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.fifi
% --------------
col1 = {'fifi'       ,'text' ,'fifi'    ,10,info.data.coef.fifi} ;
col2 = {'edit_fifi'  ,'radio','edit',0,'edit the profile',[],'','enable','off'} ;
col3 = {'import_fifi','radio','load',0,'load the profile',[],'','enable','off'} ;
col4 = {'dessin_fifi','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.rotv
% --------------
col1 = {'rotv'       ,'text' ,'rotv'    ,10,info.data.coef.rotv,''} ;
col2 = {'edit_rotv'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_rotv','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_rotv','radio','modify',0,'modify the stored profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.coef.rot
% -------------
col1 = {'rot'       ,'text' ,'rot'    ,10,info.data.coef.rot,''} ;
col2 = {'edit_rot'  ,'radio','edit',0,'edit the profile',[],''} ;
col3 = {'import_rot','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_rot','radio','modify',0,'modify the stored profile',[],''} ;
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
hout=zuicreeform('Edition','data_coef','zuiedit_data_coef_fct','',form,comm) ;


