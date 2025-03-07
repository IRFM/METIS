% ZUIEDIT_DATA_ASSER formulaire d'édition des consignes des asservissements
%
%-------------------------------------------------------------------------------
% fonction Matlab 5 :
%
% fichier zuiedit_data_asser.m ->  
%		zuicreeform : creation du formulaire
% 
%
% syntaxe  :
%  
%   	hout=zuiedit_data_asser;
%
% entree :
%
% sorties :
%   hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% 
% version 1.7, du 13/12/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
function hout=zuiedit_data_asser2

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle('data_asser') ;
if ishandle(hform)
        zuiformvisible(hform) ;
	hout = hform;
	return
end

% infos pour les tooltips
info = zinfo ;

% recuperation de la liste des consignes d'aaservissements
infoasser = info.data.cons.asser;
consasser = evalin('base','data.cons.asser');
if isempty(consasser) |isempty(infoasser)
   errordlg('no data', ...
            'Error within interface creation');
   return
end
nomasser  = fieldnames(infoasser);

% formulaire avec la sous structure from
form={};

% Titre
% -----
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text@full','Feedback control references',[],''};
form{length(form)+1} = {col1};

colj=  {'jump1','jump','jump',[],''};

% Séparation
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};


% debut de la boucle sur les champs
% cas des consigne simple
for kl = 1:length(nomasser)
     % extraction des donnees
     nomc  = nomasser{kl};
     infoc = getfield(infoasser,nomc);
     consc = getfield(consasser,nomc);

     % selon la dimension
     if size(consc,2) == 1
         % creation du formulaire associe
     	 col1 = {nomc,'text' ,nomc,10,infoc,''} ;
	 col2 = {strcat('edit_',nomc)  ,'radio','edit',0,'Edit the reference',[],''} ;
	 col3 = {strcat('import_',nomc),'radio','load',0,'load the reference',[],''} ;
	 form{length(form)+1} = {col1,col2,col3} ;
     end
end

% debut de la boucle sur les champs
% cas des consigne multiple
for kl = 1:length(nomasser)
     % extraction des donnees
     nomc  = nomasser{kl};
     infoc = getfield(infoasser,nomc);
     consc = getfield(consasser,nomc);

     % selon la dimension
     if size(consc,2) >1
	 % Séparation
	 % ----------
	 sepa ={'separation_comm','frame','',3,''};
	 form{length(form)+1} = {sepa};
         % creation du formulaire associe
         nb = size(consc,2);
	 for i=1:nb
	 	 if i==1
	 		 col1 = {nomc,'text',nomc,5,infoc,''} ;
	 	 else
	 		 col1 = {nomc,'text',' '  ,5,'',''} ;
	 	 end
	 	 st   = sprintf('channel %s',num2str(i)) ;
	 	 col2 = {'num' ,'text',st,10,'channel ??',''} ;
	 	 st   = strcat('edit_',nomc,'_',num2str(i)) ;
	 	 col3 = {st,'radio','edit',0,'edit the reference',[],''} ;
	 	 st   = strcat('prec_',nomc,'_',num2str(i)) ;
	 	 if i>1
	 		 col4 = {st,'radio','last channel',0,'Rewrite the last channel',[],''} ;
	 	 else
	 		 col4 = {st,'text',' ',0,'',[],''} ;
	 	 end
	 	 st   = strcat('import_',nomc,'_',num2str(i)) ;
	 	 col5 = {st,'radio','load',0,'load the reference',[],''} ;
	 	 form{length(form)+1} = {col1,col2,col3,col4,col5} ;
	 end
     end
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
hout=zuicreeform('Edition','data_asser','zuiedit_data_asser_fct2','',form,comm) ;


