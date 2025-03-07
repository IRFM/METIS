% ZUIEDIT_DATA_CONS formulaire principal d'�ition des consignes
%--------------------------------------------------------------------
% fichier zuiedit_data_cons.m ->  
%	zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	creation du formulaire principal d'�ition des consignes
%
% syntaxe  :
%	hout=zuiedit_data_cons;
%
% entree :
%
% sorties :
%	hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 2.0, du 10/12/2001.
% 
% liste des modifications : 
% * 14/03/2002 -> ajout de la stabilite MHD
% * 10/12/2002 -> interface anglais
%  
%--------------------------------------------------------------
function hout=zuiedit_data_cons

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle('data_cons') ;
if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

% formulaire avec la sous structure from
form={};

% Titre
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text','Plasma references',[],''};
form{length(form)+1} = {col1};

% S�aration
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% Radio boutons
col1 = {'limites' ,'radio','boundary conditions',0,'boundary condition values',[],''} ;
form{length(form)+1} = {col1};
sepa ={'separation_comm','frame','',3,''} ;

col1 = {'chauffages','radio','Heat sources',0,'Heat source values',[],''} ;
form{length(form)+1} = {col1};
sepa ={'separation_comm','frame','',3,''} ;

col1 = {'inject' ,'radio','matter (inj/ext) - plasma composition',0, ...
	'adjustement of matter injection/extraction and plasma composition',[],''} ;
form{length(form)+1} = {col1};
sepa ={'separation_comm','frame','',3,''} ;

col1 = {'mhd' ,'radio','MHD (data.cons.stab)',0, ...
	'MHD parameters (toroidal mode number)',[],''} ;
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame','',1,''};
form{length(form)+1} = {sepa};
colj = {'jump','jump',' ',[],''};
colj = {'jump','jump',' ',[],''};
form{length(form)+1} = {colj,colj};

% Bouton Quit
% ----------
comm{1}={'btn_quit','radio@center','Close',0,'close the window'};

% Formulaire
hout=zuicreeform('Edition','data_cons','zuiedit_data_cons_fct','',form,comm) ;

