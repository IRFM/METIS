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
function hout=zuiedit_data_cons_chauf_idn

% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle('data_cons_chauf_idn') ;
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
hout=zuicreeform('Edition','data_cons_chauf_idn','zuiedit_data_cons_chauf_idn_fct','',form,comm) ;


