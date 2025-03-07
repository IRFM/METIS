% ZUIEDIT_DATA_ASSER formulaire d'�ition des consignes des asservissements
%-------------------------------------------------------------------------------
%
% fichier zuiedit_data_asser.m ->  
%	zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	creation du formulaire d'�ition des consignes des asservissements
%
% syntaxe  :
%   	hout=zuiedit_data_asser;
%
% entree :
%
% sorties :
%   hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 2.0, du 10/12/2002.
% 
% liste des modifications : 
% * 10/12/2002 interface en anglais
%--------------------------------------------------------------
function hout=zuiedit_data_asser

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle('data_asser') ;
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

col1 = {'libel','text@full','Feedback control',[],''};
form{length(form)+1} = {col1};

colj=  {'jump1','jump','jump',[],''};

% S�aration
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% data.cons.asser.nl0
% -------------------
col1 = {'nl0'       ,'text' ,'nl0'    ,10,info.data.cons.asser.nl0,''} ;
col2 = {'edit_nl0'  ,'radio','edit',0,'edit the value',[],''} ;
col3 = {'import_nl0','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.asser.ne0
% -------------------
col1 = {'ne0'       ,'text' ,'ne0'    ,10,info.data.cons.asser.ne0,''} ;
col2 = {'edit_ne0'  ,'radio','edit',0,'edit the value',[],''} ;
col3 = {'import_ne0','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.asser.ne1
% -------------------
col1 = {'ne1'       ,'text' ,'ne1'    ,10,info.data.cons.asser.ne1,''} ;
col2 = {'edit_ne1'  ,'radio','edit',0,'edit the value',[],''} ;
col3 = {'import_ne1','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.asser.nemoy
% ---------------------
col1 = {'nemoy'       ,'text' ,'nemoy'   ,10,info.data.cons.asser.nemoy,''} ;
col2 = {'edit_nemoy'  ,'radio','edit',0,'edit the value',[],''} ;
col3 = {'import_nemoy','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.asser.te0
% -------------------
col1 = {'te0'       ,'text' ,'te0'    ,10,info.data.cons.asser.te0,''} ;
col2 = {'edit_te0'  ,'radio','edit',0,'edit the value',[],''} ;
col3 = {'import_te0','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.asser.te1
% -------------------
col1 = {'te1'       ,'text' ,'te1'    ,10,info.data.cons.asser.te1,''} ;
col2 = {'edit_te1'  ,'radio','edit',0,'edit the value',[],''} ;
col3 = {'import_te1','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.asser.ti0
% -------------------
col1 = {'ti0'       ,'text' ,'ti0'    ,10,info.data.cons.asser.ti0,''} ;
col2 = {'edit_ti0'  ,'radio','edit',0,'edit the value',[],''} ;
col3 = {'import_ti0','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.asser.ti1
% -------------------
col1 = {'ti1'       ,'text' ,'ti1'    ,10,info.data.cons.asser.ti1,''} ;
col2 = {'edit_ti1'  ,'radio','edit',0,'edit the value',[],''} ;
col3 = {'import_ti1','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.asser.li
% ------------------
col1 = {'li'       ,'text' ,'li'    ,10,info.data.cons.asser.li,''} ;
col2 = {'edit_li'  ,'radio','edit',0,'edit the value',[],''} ;
col3 = {'import_li','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% data.cons.asser.vloop	=> N'EXISTE PAS 
% ---------------------
%col1 = {'vloop'       ,'text' ,'vloop'    ,10,info.data.cons.asser.vloop,''} ;
%col2 = {'edit_vloop'  ,'radio','edit',0,'edit the value',[],''} ;
%col3 = {'import_vloop','radio','load',0,'load the value',[],''} ;
%form{length(form)+1} = {col1,col2,col3} ;

% data.cons.asser.q0
% ------------------
col1 = {'q0'       ,'text' ,'q0'    ,10,info.data.cons.asser.q0,''} ;
col2 = {'edit_q0'  ,'radio','edit',0,'edit the value',[],''} ;
col3 = {'import_q0','radio','load',0,'load the value',[],''} ;
form{length(form)+1} = {col1,col2,col3} ;

% S�aration
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% data.cons.asser.c
% -----------------
nb = evalin('base','param.gene.nbg') ;
for i=1:nb
	if i==1
		col1 = {'c','text','c',5,info.data.cons.c,''} ;
	else
		col1 = {'c','text',' '  ,5,'',''} ;
	end
	st   = sprintf('channel %s',num2str(i)) ;
	col2 = {'num' ,'text',st,10,'number ?',''} ;
	st   = strcat('edit_module_c',num2str(i)) ;
	col3 = {st,'radio','edit',0,'edit the value',[],''} ;
	st   = strcat('prec_c',num2str(i)) ;
	if i>1
		col4 = {st,'radio','next channel',0,'copy last channel',[],''} ;
	else
		col4 = {st,'text',' ',0,'',[],''} ;
	end
	st   = strcat('import_c',num2str(i)) ;
	col5 = {st,'radio','load',0,'load the value',[],''} ;
	form{length(form)+1} = {col1,col2,col3,col4,col5} ;
end


% separation
% ----------
sepa = {'separation_comm','frame','',1,''};
form{length(form)+1} = {sepa};
colj = {'jump','jump',' ',[],''};
colj = {'jump','jump',' ',[],''};
form{length(form)+1} = {colj,colj};

% Bouton Quit
% ----------
comm{1}={'btn_quit','radio@center','Closer',0,'To close the window'};

% Formulaire
% ----------
hout=zuicreeform('Edition','data_asser','zuiedit_data_asser_fct','',form,comm) ;


