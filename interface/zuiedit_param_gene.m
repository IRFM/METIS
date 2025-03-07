% ZUIEDIT_PARAM_GENE formulaire 'parametres g��aux' sous 'Edition'
%--------------------------------------------------------------
% fichier zuiedit_param_gene.m ->
%		zuicreeform	: creation du formulaire
%
% fonction Matlab 5 :
%	fonction de creation de GUI pour le formulaire
%	des parametres  g��aux
% 	sous le mode edition du formulaire principal
%
% syntaxe  :
%  zuiedit_parma_gene ;
%
% entrees :
%
% sorties :
% 			
% fonction ecrite par C. Passeron, poste 61 19
% version 2.0, du 11/12/2002.
% 
% liste des modifications : 
% * 03/10/2001 -> ajout de la configuration deu profiler (J-F Artaud)
% * 11/12/2002 -> interface en anglais 
% 
%--------------------------------------------------------------

function zuiedit_param_gene

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('ed_param_gene');

if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

info=zinfo ;

form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text','Global simulation parameters',[],''};
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

col1 = {'radio_conv','radio' ,'convergence',0,'solver convergence parameter'};
form{length(form)+1} = {col1};

%col1 = {'radio_multi','radio' ,'multiplicator',0,'neoclassical coefficient multiplicator'};
%form{length(form)+1} = {col1};

col1 = {'radio_sauv','radio' ,'save',0,'ouput file parameters'};
form{length(form)+1} = {col1};

col1 = {'radio_info','radio' ,'information',0,'inside CRONOS parameter used during the run'};
form{length(form)+1} = {col1};

col1 = {'radio_config','radio' ,'equations',0,'parameters for diffusion equations '};
form{length(form)+1} = {col1};

%col1 = {'radio_exec','radio' ,'execution',0,'Initial and final time of the run'};
%form{length(form)+1} = {col1};

col1 = {'radio_profile','radio' ,'profiler',0,'"profiler" MatLab parameter'};
form{length(form)+1} = {col1};

col1 = {'radio_from','radio' ,'origin',0,'origin of the simulation '};
form{length(form)+1} = {col1};

%--------------------------------------------------------------------
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

hout=zuicreeform('Edition','ed_param_gene','zuiedit_param_gene_fct','',form,comm);
