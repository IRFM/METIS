% ZUIEDIT_DATA_CONS_CHAUF formulaire d'�ition des consignes de chauffage
%-------------------------------------------------------------------------
% fichier zuiedit_data_cons_chauf.m ->  
%		zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	creation du formulaire d'�ition des consignes de chauffage
%
% syntaxe  :
%	hout=zuiedit_data_cons_chauf;
%
% entree :
%
% sorties :
%	hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 28/01/2004.
% 
% liste des modifications : 
% * 10/12/2002 : interface en anglais
% * 15/09/2003 -> edition angle poloidal et toroidal pour fce
% * 28/01/2004 -> suppression de la consigne dPdt pour IDN (plus utilisee)
%--------------------------------------------------------------
function hout=zuiedit_data_cons_chauf

% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle('data_cons_chauf') ;
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

col1 = {'libel','text@full','heat source references',[],''};
form{length(form)+1} = {col1};

colj=  {'jump1','jump','jump',[],''};

% S�aration
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};


% data.cons.fci
% -------------
nb = evalin('base','param.nombre.fci') ;
for i=1:nb
	if i==1
		col1 = {'fci','text','icrh',10,info.data.cons.fci,''} ;
	else
		col1 = {'fci','text',' '  ,10,'',''} ;
	end
	st   = sprintf('channel %s',num2str(i)) ;
	col2 = {'num' ,'text',st,10,'number ?',''} ;
	st   = strcat('edit_module_fci',num2str(i)) ;
	col3 = {st,'radio','edit power',0,'Edit the power module',[],''} ;
	st   = strcat('edit_phase_fci',num2str(i)) ;
	col4 = {st ,'radio','edit phase ',0,'Edit the phase module',[],''} ;
	st   = strcat('prec_fci',num2str(i)) ;
	if i>1
		col5 = {st,'radio','previous channel',0,'copy the previous channel',[],''} ;
	else
		%col5 = {st,'radio','canal pr��ent',0,'Recopie le canal pr��ent',[],'','enable','off'} ;
		col5 = {st,'text',' ',0,'',[],''} ;
	end
	st   = strcat('import_fci',num2str(i)) ;
	col6 = {st,'radio','load',0,'load value',[],''} ;
	form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;
end

% data.cons.fce
% -------------
nb = evalin('base','param.nombre.fce') ;
for i=1:nb
	if i==1
		col1 = {'fce','text','ecrh',10,info.data.cons.fce,''} ;
	else
		col1 = {'fce','text',' '  ,10,'',''} ;
	end
	st   = sprintf('channel %s',num2str(i)) ;
	col2 = {'num' ,'text',st,10,'number ?',''} ;
	st   = strcat('edit_module_fce',num2str(i)) ;
	col3 = {st,'radio','edit power',0,'Edit the power module',[],''} ;
	st   = strcat('edit_phase_fce',num2str(i)) ;
	col4 = {st,'radio','edit angles ',0,'Edit the toroidal and poloidal launched angles',[],''} ;
	st   = strcat('prec_fce',num2str(i)) ;
	if i>1
		col5 = {st,'radio','previous channel',0,'copy the previous channel',[],''} ;
	else
		col5 = {st,'text',' ',0,'',[],''} ;
	end
	st   = strcat('import_fce',num2str(i)) ;
	col6 = {st,'radio','load',0,'load value',[],''} ;
	form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;
end

% data.cons.hyb
% -------------
nb = evalin('base','param.nombre.hyb') ;
for i=1:nb
	if i==1
		col1 = {'hyb','text','lh',10,info.data.cons.hyb,''} ;
	else
		col1 = {'hyb','text',' '  ,10,'',''} ;
	end
	st   = sprintf('channel %s',num2str(i)) ;
	col2 = {'num' ,'text',st,10,'number ?',''} ;
	st   = strcat('edit_module_hyb',num2str(i)) ;
	col3 = {st,'radio','edit power',0,'Edit the power module',[],''} ;
	st   = strcat('edit_phase_hyb',num2str(i)) ;
	col4 = {st,'radio','edit N // ',0,'Edit the parallel wave number',[],''} ;
	st   = strcat('prec_hyb',num2str(i)) ;
	if i>1
		col5 = {st,'radio','previous channel',0,'copy the previous channel',[],''} ;
	else
		col5 = {st,'text',' ',0,'',[],''} ;
	end
	st   = strcat('import_hyb',num2str(i)) ;
	col6 = {st,'radio','load',0,'load value',[],''} ;
	form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;
end

% data.cons.idn
% -------------
nb = evalin('base','param.nombre.idn') ;
for i=1:nb
	if i==1
		col1 = {'idn','text','nbi',10,info.data.cons.idn,''} ;
	else
		col1 = {'idn','text',' '  ,10,'',''} ;
	end
	st   = sprintf('channel %s',num2str(i)) ;
	col2 = {'num' ,'text',st,10,'number ?',''} ;
	st   = strcat('edit_module_idn',num2str(i)) ;
	col3 = {st,'radio','power',0,'Edit the power module',[],''} ;
	st   = strcat('edit_phase_idn',num2str(i)) ;
	col4 = {st,'radio','- ',0,'Not used',[],''} ;
	st   = strcat('prec_idn',num2str(i)) ;
	if i>1
		col5 = {st,'radio','previous channel',0,'copy the previous channel',[],''} ;
	else
		col5 = {st,'text',' ',0,'',[],''} ;
	end
	st   = strcat('import_idn',num2str(i)) ;
	col6 = {st,'radio','load',0,'load value',[],''} ;
	form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;
end

% data.cons.ext
% -------------
%nb = evalin('base','param.nombre.ext') ;
%for i=1:nb
%	if i==1
%		col1 = {'idn','text','ext',10,info.data.cons.ext,''} ;
%	else
%		col1 = {'idn','text',' '  ,10,info.data.cons.ext,''} ;
%	end
%	st   = strcat('canal ',num2str(i)) ;
%	col2 = {'num' ,'text',st,10,'num�o ??',''} ;
%	st   = strcat('edit_module_ext',num2str(i)) ;
%	col3 = {'edit_module_ext','radio','editer module',0,'Editer le module de la consigne',[],''} ;
%	st   = strcat('edit_phase_ext',num2str(i)) ;
%	col4 = {'edit_phase_ext' ,'radio','editer phase ',0,'Editer la phase de la consigne',[],''} ;
%	st   = strcat('prec_ext',num2str(i)) ;
%	if i>1
%		col5 = {'prec_ext','radio','canal pr��ent',0,'Recopie le canal pr��ent',[],''} ;
%	else
%		col5 = {'prec_ext','radio','canal pr��ent',0,'Recopie le canal pr��ent',[],'','enable','off'} ;
%	end
%	st   = strcat('import_ext',num2str(i)) ;
%	col6 = {'import_ext','radio','importer',0,'Importer la consigne',[],''} ;
%	form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;
%end

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
hout=zuicreeform('Edition','data_cons_chauf','zuiedit_data_cons_chauf_fct','',form,comm) ;


