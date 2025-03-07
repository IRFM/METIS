% ZUIEDIT_PARAM_FROM   formulaire 'informations' sous 'Edition'
%--------------------------------------------------------------
% fichier zuiedit_param_from.m ->  
%		zuicreeform : creation du formulaire
%		zuiedit_param_from_fct	: fonction de control des callbacks
%		zuiedit_param_from_crtl	: fonction de control des donnees
% 
% fonction Matlab 5 :
%	fonction de creation de GUI pour le formulaire
%	d'informations (FROM) , commandes de parametres  
%	sous le mode edition du formulaire principal
%
% syntaxe  :
%   	zuiedit_param_from ;
%
% entrees :
%
% sorties :
%
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 28/08/2001.
% 
% liste des modifications : 
%  * 28/08/2001   -> ajout de l'affichage du cr�teur (J-F Artaud)
%  * 27/09/2001   -> modification de chaque appel 'evalin'
%                    on initialise la variable a vide si elle n'exite pas dans la base :
%                    evalin('base','param.from....','[]')
%
%--------------------------------------------------------------

function zuiedit_param_from

[hfig,h] = zuiformhandle('from');
if ishandle(hfig)
	zuiformvisible(hfig) ;
	return
end

info=zinfo ;

%--------------------------------------------------------------------
% le formulaire
nb1 = 10;
nb2 = 15 ;
form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel_loadfile','text@full','Global simulation parameters ',[],''};
form{length(form)+1} = {col1};

col1 = {'libel_loadfile','text@full','informations (from)',[],''};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% Informations : param.from  
% -------------------------
colj = {'void','jump','void',[],''} ;

val  = evalin('base','param.from.machine','[]');
col1 = {'libel','text','machine',nb1,info.param.from.machine,''};
col2 = {'libel','text',val,nb2,info.param.from.machine,''};
form{length(form)+1} = {col1,col2};

numshot = evalin('base','param.from.shot.num','[]') ;
numchoc   = fix(numshot) ;
val  = evalin('base','param.from.shot.num','[]');
col1 = {'libel','text','numero de choc',nb1,info.param.from.shot.num,''};
col2 = {'libel','text',numchoc,nb2,info.param.from.shot.num,''};
form{length(form)+1} = {col1,col2};

% decodage de la date
col1 = {'libel','text','date du choc',nb1,info.param.from.shot.date,''};
try
	an    = evalin('base','param.from.shot.date(1)') ;
	mois  = evalin('base','param.from.shot.date(2)') ;
	jour  = evalin('base','param.from.shot.date(3)') ;
	heure = evalin('base','param.from.shot.date(4)') ;
	min   = evalin('base','param.from.shot.date(5)') ;
	sec   = fix(evalin('base','param.from.shot.date(6)')) ;

	val  = sprintf('%d/%d/%d �%d:%d:%d',jour,mois,an,heure,min,sec) ;
	col2 = {'libel','text',val,nb2,info.param.from.shot.date,''};
catch
	date = evalin('base','param.from.shot.date','[]') ;
	col2 = {'libel','text',date,nb2,info.param.from.shot.date,''};
end
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.from.shot.info','[]') ;
if ~isempty(val)
	if isstruct(val)
   		col1 = {'libel','text','shot information',nb1,info.param.from.shot.info,''};
		col2 = {'radio_info','radio','View informations',0,' Open information viewer'} ;
	else
   		col1 = {'libel','text','info. sur choc',nb1,info.param.from.shot.info,''};
   		col2 = {'libel','text',val,nb2,info.param.from.shot.info,''};
	end
   	form{length(form)+1} = {col1,col2};
end

% decodage de la date
col1 = {'libel','text','file creation data',nb1,info.param.from.creation.date,''};
try
	jour  = evalin('base','param.from.creation.date(3)') ;
	mois  = evalin('base','param.from.creation.date(2)') ;
	an    = evalin('base','param.from.creation.date(1)') ;
	heure = evalin('base','param.from.creation.date(4)') ;
	min   = evalin('base','param.from.creation.date(5)') ;
	sec   = fix(evalin('base','param.from.creation.date(6)')) ;

	val  = sprintf('%d/%d/%d �%d:%d:%d',jour,mois,an,heure,min,sec) ;
	col2 = {'libel','text',val,nb2,info.param.from.creation.date,''};
catch
	date = evalin('base','param.from.creation.date','[]') ;
	col2 = {'libel','text',date,nb2,info.param.from.creation.date,''};
end
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.from.creation.user','[]') ;
if ~strcmp(val,'')
   ind = find(char(val)==10) ;
   if ~isempty(ind)
      val = val(1:ind-1) ;
   end
   col1 = {'libel','text','nom user',nb1,info.param.from.creation.user,''};
   col2 = {'user','edit',val,nb2,info.param.from.creation.user,[],'param.from.creation.user'};
   form{length(form)+1} = {col1,col2};
end

val  = evalin('base','param.from.creation.info','[]') ;
if ~isempty(val)
	if isstruct(val)
   		col1 = {'libel','text','info. creation',nb1,info.param.from.creation.info,''};
		col2 = {'radio_creainfo','radio','View informations',0,' Open information viewer'} ;
   		form{length(form)+1} = {col1,col2};
	else
		if ~strcmp(val,'')
   			col1 = {'libel','text','Data source information',nb1,info.param.from.creation.info,''};
   			col2 = {'libel','text',val,nb2,info.param.from.creation.info,''};
   			form{length(form)+1} = {col1,col2};
		end
	end
end

val  = evalin('base','param.from.creation.com','[]') ;
%if ~strcmp(val,'')
	col1 = {'libel','text','comments',nb1,info.param.from.creation.com,''};
	col2 = {'libel','edit',val,nb2,info.param.from.creation.com,[],'param.from.creation.com'};
	form{length(form)+1} = {col1,col2};
%end

val  = evalin('base','param.from.source.desc','[]') ;
if ~strcmp(val,'')
   	col1 = {'libel','text','description',nb1,info.param.from.source.desc,''};
	desc = '-' ;
   	for i=1:length(val)
		desc = strcat(desc,char(val(i)),'-') ;
	end
   	col2 = {'libel','text',desc,nb2,info.param.from.source.desc,''};
   	form{length(form)+1} = {col1,col2};
end

% ajout J-F Artaud du 28/08/2001
% info sur le createur
try
   val  = evalin('base','param.from.createur') ;
catch
   val ='';
end

if ~isempty(val)
	col1 = {'libel_createur','text','Data constructor',nb1,info.param.from.createur,''};
	col2 = {'radio_createur','radio',val,0,'Display constructor parameters'} ;
	form{length(form)+1} = {col1,col2};
end

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% param.from.sample.signal
% ------------------------
nb1 = 10 ;
nb2 = 10;
col1 = {'signal','text@full','temporal signal',1,info.param.from.sample.signal,''};
form{length(form)+1} = {col1} ;

val = evalin('base','param.from.sample.signal.ondelette','[]') ;
col1 = {'libel','text','wavelet',nb1,'',''};
col2 = {'libel','text',val,nb2,'',''};

tps = evalin('base','param.from.sample.signal.defaut.temps','[]') ;
if isnan(tps)
	tps = 'NaN' ;
end
col3 = {'libel','text','time',nb1,'',''};
col4 = {'libel','text',tps,nb2,'',''};

val = evalin('base','param.from.sample.signal.defaut.espace','[]') ;
col5 = {'libel','text','defaut.espace',nb1,info.param.from.sample.signal,''};
col6 = {'libel','text',val,nb2,'',''};

form{length(form)+1} = {col1,col2,col3,col4,col5,col6};

val = evalin('base','param.from.sample.signal.defaut.inf','[]') ;
%if isempty(val)
%	val = '[ ]' ;
%end
col1 = {'libel','text','defaut.inf',nb1,info.param.from.sample.signal,''};
col2 = {'libel','text',val,nb2,'',''};

val = evalin('base','param.from.sample.signal.plus','[]') ;
col3 = {'libel','text','plus',nb1,info.param.from.sample.signal,''};
col4 = {'libel','text',val,nb2,'',''};

col5 = {'void','jump','void',[],''} ;
col6 = {'void','jump','void',[],''} ;

form{length(form)+1} = {col1,col2,col3,col4,col5,col6};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% param.from.sample.groupe
% ------------------------
col1 = {'groupe','text@full','profile signal',1,info.param.from.sample.groupe,''};
form{length(form)+1} = {col1} ;

val = evalin('base','param.from.sample.groupe.ondelette','[]') ;
col1 = {'libel','text','wavelet',nb1,'',''};
col2 = {'libel','text',val,nb2,'',''};

tps = evalin('base','param.from.sample.groupe.defaut.temps','[]') ;
if isnan(tps)
	tps = 'NaN' ;
end
col3 = {'libel','text','time',nb1,'',''};
col4 = {'libel','text',tps,nb2,'',''};

val = evalin('base','param.from.sample.groupe.defaut.espace','[]') ;
col5 = {'libel','text','defaut.espace',nb1,info.param.from.sample.groupe,''};
col6 = {'libel','text',val,nb2,'',''};

form{length(form)+1} = {col1,col2,col3,col4,col5,col6} ;

val = evalin('base','param.from.sample.groupe.defaut.inf','[]') ;
%if isempty(val)
%	val = '[ ]' ;
%end
col1 = {'libel','text','defaut.inf',nb1,info.param.from.sample.groupe,''};
col2 = {'libel','text',val,nb2,'',''};

val = evalin('base','param.from.sample.groupe.plus','[]') ;
col3 = {'libel','text','plus',nb1,info.param.from.sample.groupe,''};
col4 = {'libel','text',val,nb2,'',''};

val = evalin('base','param.from.sample.groupe.energie','[]') ;
col5 = {'libel','text','threshold energy',nb1,info.param.from.sample.groupe,''};
col6 = {'libel','text',val,nb2,'',''};

form{length(form)+1} = {col1,col2,col3,col4,col5,col6};

hout=zuicreeform('Edition','from','zuiedit_param_from_fct','',form);
