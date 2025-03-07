% ZUIEDIT_DATA_CONS_INJ formulaire d'édition des consignes d'injection/extraction matiere
%                       et réglage de la composition du plasma
%
%-------------------------------------------------------------------------
% fichier zuiedit_data_cons_chauf.m ->  
%		zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	 d'édition des consignes d'injection/extraction matiere
%	et réglage de la composition du plasma
%
% syntaxe  :
%	hout=zuiedit_data_cons_inj;
%
% entree :
%
% sorties :
%   hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 15/03/2003.
% 
% liste des modifications : 
%
% * 30/08/2001 -> correction hout et hform (J-F Artaud)
% * 18/12/2002 -> Interface en anglais
% * 15/09/2003 -> ajout de ntnd
%
%--------------------------------------------------------------
function hout=zuiedit_data_cons_inj

% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle('data_cons_inj') ;
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

col1 = {'libel','text@full','Plasma composition',[],''};
form{length(form)+1} = {col1};

col1 = {'libel','text@full','particules sources/Zeff',[],''};
form{length(form)+1} = {col1};

colj=  {'jump1','jump','jump',[],''};

% Séparation
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};


% data.cons.c
% -----------
%nb = evalin('base','param.gene.nbg') ;
%for i=1:nb
%	if i==1
%		col1 = {'c','text','c',10,info.data.cons.c,''} ;
%	else
%		col1 = {'c','text',' '  ,10,'',''} ;
%	end
%	st   = sprintf('channel %s',num2str(i)) ;
%	col2 = {'num' ,'text',st,10,'Channel ?',''} ;
%	st   = strcat('edit_module_c',num2str(i)) ;
%	col3 = {st,'radio','Edit',0,'Edit the initial value',[],''} ;
%	st   = strcat('prec_c',num2str(i)) ;
%	if i>1
%		col4 = {st,'radio','previous channel',0,'Put the previous value',[],''} ;
%	else
%		%col4 = {st,'radio','canal précédent',0,'Recopie le canal précédent',[],'','enable','off'} ;
%		col4 = {st,'text',' ',0,'',[],''} ;
%	end
%	st   = strcat('import_c',num2str(i)) ;
%	col5 = {st,'radio','load',0,'load initila value',[],''} ;
%	form{length(form)+1} = {col1,col2,col3,col4,col5} ;
%end

% data.cons.pomp
% --------------
%col1 = {'pomp','text','pomp',10,info.data.cons.pomp,''} ;
%col2 = {'txt' ,'text',' ',0,'',[],''} ;
%col3 = {'edit_module_pomp','radio','Edit',0,'Edit the initial value',[],''} ;
%col4 = {'prec_pomp','text',' ',0,'',[],''} ;
%col5 = {'import_pomp','radio','load',0,'Put the previous value',[],''} ;
%form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.cons.glacon
% ----------------
%nb = evalin('base','param.nombre.glacon') ;
%for i=1:nb
%	if i==1
%		col1 = {'pellet','text','pellet',10,info.data.cons.glacon,''} ;
%	else
%		col1 = {'pellet','text',' '  ,10,'',''} ;
%        end
%	st   = sprintf('channel %s',num2str(i)) ;
%	col2 = {'num' ,'text',st,10,'number ?',''} ;
%	st   = strcat('edit_module_glacon',num2str(i)) ;
%	col3 = {st,'radio','Edit',0,'Edit the initial value',[],''} ;
%	st   = strcat('prec_glacon',num2str(i)) ;
%	if i>1
%		col4 = {st,'radio','previous channel',0,'Put the previous value',[],''} ;
%	else
%		col4 = {st,'text',' ',0,'',[],''} ;
%	end
%	st   = strcat('import_glacon',num2str(i)) ;
%	col5 = {st,'radio','load',0,'load the initial value',[],''} ;
%	form{length(form)+1} = {col1,col2,col3,col4,col5} ;
%end

% data.cons.zeffm
% ---------------
col1 = {'zeffm','text','zeffm',10,info.data.cons.zeffm,''} ;
col2 = {'txt','text',' ',0,'',[],''} ;
col3 = {'edit_module_zeffm','radio','Edit',0,'Edit the initila value',[],''} ;
col4 = {'prec_zeffm','text',' ',0,'',[],''} ;
col5 = {'import_zeffm','radio','load',0,'load the initial value',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.cons.nhnd
% --------------
col1 = {'nhnd','text','nhnd',10,info.data.cons.nhnd,''} ;
col2 = {'txt' ,'text',' ',0,'',[],''} ;
col3 = {'edit_module_nhnd','radio','Edit',0,'Edit the initial value',[],''} ;
col4 = {'prec_nhnd','text',' ',0,'',[],''} ;
col5 = {'import_nhnd','radio','load',0,'load the initial value',[],''} ;
form{length(form)+1} = {col1,col2,col3,col4,col5} ;

% data.cons.ntnd
% --------------
col1 = {'ntnd','text','ntnd',10,info.data.cons.nhnd,''} ;
col2 = {'txt' ,'text',' ',0,'',[],''} ;
col3 = {'edit_module_ntnd','radio','Edit',0,'Edit the initial value',[],''} ;
col4 = {'prec_ntnd','text',' ',0,'',[],''} ;
col5 = {'import_ntnd','radio','load',0,'load the initial value',[],''} ;
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
hout=zuicreeform('Edition','data_cons_inj','zuiedit_data_cons_inj_fct','',form,comm) ;


