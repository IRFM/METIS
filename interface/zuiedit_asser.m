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

function zuiedit_asser

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('ed_asser');

if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

info=zinfo ;

form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text','Feedback control',[],''};
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

col1 = {'radio_assercalc','radio' ,'calculation mode',0,''};
form{length(form)+1} = {col1};

col1 = {'radio_asserref','radio' ,'reference values',0,''};
form{length(form)+1} = {col1};

col1 = {'radio_simulink'    ,'radio' ,'Simulink'                ,0 ,'Open Simulink @ Cronos dialog',[],''};
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

hout=zuicreeform('Edition','ed_asser','zuiedit_asser_fct','',form,comm);
