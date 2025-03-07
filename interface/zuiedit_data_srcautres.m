% ZUIEDIT_DATA_SRCAUTRES formulaire d'�ition des profils source rayonnement et n�classique
%-------------------------------------------------------------------------------------------
% fichier zuiedit_data_srcautres.m ->  
%		   zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	creation du formulaire d'�ition des profils source rayonnement et n�classique
%
% syntaxe  :
%	hout=zuiedit_data_srcautres;
%
% entree :
%
% sorties :
%	hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 2.0, du 11/12/2002.
% 
% liste des modifications : 
% * 11/12/2002 : interface en anglais
%--------------------------------------------------------------
function hout=zuiedit_data_srcautres

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle('data_srcautres') ;
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

col1 = {'libel','text@full','radiative and neoclassical source',[],''};
form{length(form)+1} = {col1};

% S�aration
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% data.source.prad
% ----------------
col1 = {'prad'       ,'text' ,'prad',10,info.data.source.prad,''} ;
col2 = {'edit_prad'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_prad','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_prad','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.brem
% ----------------
col1 = {'brem'       ,'text' ,'brem',10,info.data.source.brem,''} ;
col2 = {'edit_brem'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_brem','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_brem','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

% data.source.cyclo
% -----------------
col1 = {'cyclo'       ,'text' ,'cyclo',10,info.data.source.cyclo,''} ;
col2 = {'edit_cyclo'  ,'radio','edit',0  ,'edit the profile',[],''} ;
col3 = {'import_cyclo','radio','load',0,'load the profile',[],''} ;
col4 = {'dessin_cyclo','radio','modify',0,'modify the profile',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4} ;

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
comm{1}={'btn_quit','radio@center','Close',0,'to close the window'};

% Formulaire
% ----------
hout=zuicreeform('Edition','data_srcautres','zuiedit_data_srcautres_fct','',form,comm) ;


