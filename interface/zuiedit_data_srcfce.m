%  ZUIEDIT_DATA_SRCFCE formulaire d'�ition des profils source-fce
%%------------------------------------------------------------------
% fichier zuiedit_data_srcfce.m ->  
%		    zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	creation du formulaire d'edition des profils source-FCE
%
% syntaxe  :
%	hout=zuiedit_data_srcfce;
%
% entree :
%
% sorties :
%	hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 03/09/2003.
% 
% liste des modifications : 
% * 11/12/2002 : interface en anglais
% * 02/09/2003 -> ajout des editions de source.??.q %
%--------------------------------------------------------------
function hout=zuiedit_data_srcfce

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle('data_srcfce') ;
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

col1 = {'libel','text@full','ERCH source',[],''};
form{length(form)+1} = {col1};

% S�aration
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% data.source.fce.el
% ------------------
col1 = {'el'       ,'text' ,'el',10,info.data.source.fce.el,''} ;
col2 = {'edit_el'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_el','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_el','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.fce.ion
% -------------------
col1 = {'ion'       ,'text' ,'ion',10,info.data.source.fce.ion,''} ;
col2 = {'edit_ion'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_ion','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ion','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.fce.ne
% ------------------
col1 = {'ne'       ,'text' ,'ne',10,info.data.source.fce.ne,''} ;
col2 = {'edit_ne'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_ne','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_ne','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.fce.j
% -----------------
col1 = {'j'       ,'text' ,'j',10,info.data.source.fce.j,''} ;
col2 = {'edit_j'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_j','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_j','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.fce.w
% -----------------
col1 = {'w'       ,'text' ,'w',10,info.data.source.fce.w,''} ;
col2 = {'edit_w'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_w','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_w','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.fce.wb
% -----------------
col1 = {'wb'       ,'text' ,'wb',10,info.data.source.fce.wb,''} ;
col2 = {'edit_wb'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_wb','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_wb','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.fce.q
% -----------------
col1 = {'q'       ,'text' ,'q',10,info.data.source.fce.q,''} ;
col2 = {'edit_q'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_q','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_q','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.fce.fluce
% ---------------------
col1 = {'fluce'       ,'text' ,'fluce',10,info.data.source.fce.fluce,''} ;
col2 = {'edit_fluce'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_fluce','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_fluce','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.fce.flucion
% -----------------------
col1 = {'flucion'       ,'text' ,'flucion',10,info.data.source.fce.flucion,''} ;
col2 = {'edit_flucion'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_flucion','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_flucion','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.fce.psupra
% ----------------------
col1 = {'psupra'       ,'text' ,'psupra',10,info.data.source.fce.psupra,''} ;
col2 = {'edit_psupra'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_psupra','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_psupra','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.fce.paniso
% ----------------------
col1 = {'paniso'       ,'text' ,'paniso',10,info.data.source.fce.paniso,''} ;
col2 = {'edit_paniso'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_paniso','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_paniso','radio','modify',0,'modify the profile',[],''} ;
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
hout=zuicreeform('Edition','data_srcfce','zuiedit_data_srcfce_fct','',form,comm) ;


