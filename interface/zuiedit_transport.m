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

function zuiedit_transport

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('ed_tansport');

if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

info=zinfo ;

form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text','Transport equations',[],''};
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

col1 = {'radio_transmode'   ,'radio' ,'computation mode',0 ,''};
form{length(form)+1} = {col1};

col1 = {'limites' ,'radio','boundary conditions values',0,'boundary condition values',[],''} ;
form{length(form)+1} = {col1};

%col1 = {'radio_coef','radio' ,'Transport coefficients' ,0 ,'Edit transport coefficients'};
%form{length(form)+1} = {col1};

col1 = {'transport_coeff_calc','radio' ,'coefficient model',0,''};
form{length(form)+1} = {col1};

col1 = {'transport_coeff_mode','radio' ,'coefficient computation mode' ,0,''};
form{length(form)+1} = {col1};

col1 = {'transport_coeff_prof','radio','coefficient profile' ,0,''};
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

hout=zuicreeform('Edition','ed_transport','zuiedit_transport_fct','',form,comm);
