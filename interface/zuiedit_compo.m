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
%--------------------------------------------------------------

function zuiedit_compo

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('ed_plasma_composition');

if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

info=zinfo ;

form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel','text','Plasma composition',[],''};
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

col1 = {'radio_compo','radio' ,'ions species',0,'ions species parameters'};
form{length(form)+1} = {col1};

col1 = {'inject' ,'radio','ions densities',0, ...
        'adjustement of matter injection/extraction and plasma composition',[],''} ;
form{length(form)+1} = {col1};

col1  = {'radio_zeff'        ,'radio','zeff computing mode'   ,0,info.data.mode.zeff,[                     ],''} ;
form{length(form)+1} = {col1};

col1  = {'radio_ae'          ,'radio','ae computing mode'     ,0,info.data.mode.ae,[],''} ;
form{length(form)+1} = {col1};

col1  = {'radio_prad'          ,'radio','prad computing mode'     ,0,info.data.mode.prad,[],''} ;
form{length(form)+1} = {col1};

col1  = {'radio_brem'          ,'radio','brem computing mode'     ,0,info.data.mode.brem,[],''} ;
form{length(form)+1} = {col1};

col1 = {'radio_scal'   ,'radio' ,'Composition Update'         ,0 ,'Update the plasma composition, ion pressure and density when zeff or gaz composition have been modified ',[],''};
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
%comm{1}={'btn_quit','radio@center','Close',0,'close the window'};
comm{1}={'btn_quit','radio@center','Close',0,''};

hout=zuicreeform('Edition','ed_plasma_composition','zuiedit_compo_fct','',form,comm);
