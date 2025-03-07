% ZUIEDIT  fonction de creation du formulaire "Edition"
%--------------------------------------------------------------------------
%
% fichier zuiedit.m ->  
%			zuicreeform	: creation du formulaire
%			zuiedit_fct	: gestion des callbacks des uicontrols
%					  du formulaire
%
%
% fonction Matlab 5 :
% 
% Cette fonction cree la feuille principale/edition  (interface graphique) sous 
% la forme d'un formulaire. Le formulaire est constitue de lignes.
% Chaque ligne peut avoir un nombre de colonnes quelconques. Les lignes
% qui se suivent et qui ont un nombre de colonnes egal, voient leurs colonnes
% alignees verticalement.
%
% syntaxe  :
%	zuiedit ;
%
% entrees
%
% sorties :
%  				
% fonction ecrite par C. Passeron, poste 61 19
% version 3.0, du 18/12/2004.
% 
% liste des modifications : 
%  * 11/12/2002 -> ajout de la mise a jour de la composition plasma
%  * 18/12/2004 -> ajout du botton pour simulink
%  * 05/02/2006 -> remaniement de l'interface
%--------------------------------------------------------------

function zuiedit

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('edit');
if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

% 1eres lignes : Rappel du nom du fichier de travail
% --------------------------------------------------
form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel_loadfile','text@full',' Filename',[],''};
form{length(form)+1} = {col1};

filename = evalin('base','param.edit.currentfile','[]') ;
if isempty(filename)
	herror = errordlg('you must load a file','CAUTION');
	zwaitfor(herror) ;
	zuicloseone(hform) ;
	% on revient a la fenetre zuidirect
	[hfig,h] = zuiformhandle('direct');
	if ishandle(hfig)
		zuiformvisible(hfig) ;	
	end
	return
end

col1 = {'text_loadfile','text',filename,[],'Filename',[],[]};
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% 1er frame : Commandes de Parametre
% ----------------------------------
col1 = {'libel','text','PARAMETERS',[],''};
form{length(form)+1} = {col1};

% ---------------------------------------------
col1 = {'radio_time'    ,'radio' ,'Time'    ,0,'time parameters'};
form{length(form)+1} = {col1};

col1 = {'radio_compo'   ,'radio' ,'Plasma composition'        ,0,'Composition parameters'};
form{length(form)+1} = {col1};

col1 = {'radio_source'    ,'radio' ,'Sources'   ,0,'source term functions'};
form{length(form)+1} = {col1};

col1 = {'radio_gene'    ,'radio' ,'Global simulation parameters'      ,0,'Global simulation parameters of CRONOS'};
form{length(form)+1} = {col1};

col1 = {'radio_equili'    ,'radio' ,'Equilibrium'      ,0,'Equilibrium parameters'};
form{length(form)+1} = {col1};

col1 = {'radio_neo'    ,'radio' ,'Neoclassical'      ,0,'Neoclassical parameters'};
form{length(form)+1} = {col1};

col1 = {'radio_mhdyna'    ,'radio' ,'MHD'      ,0,'MHD parameters'};
form{length(form)+1} = {col1};


%col1 = {'radio_fonction','radio' ,'External modules'   ,0,'External modules parameters'};
%form{length(form)+1} = {col1};

%col1 = {'radio_split'   ,'radio' ,'Time splitting' ,0,'Time splitting'};
%form{length(form)+1} = {col1};

%col1 = {'radio_from'    ,'radio' ,'Information'    ,0,'Information'};
%form{length(form)+1} = {col1};

%col1 = {'radio_intv'    ,'radio' ,'Time intervals'  ,0,'Define various time intervals for the calculation of external modules'};
%form{length(form)+1} = {col1};


% separation
% ----------
%sepa = {'separation_comm','frame','',-5,''};
%form{length(form)+1} = {sepa};
%sepa = {'separation_comm','frame','',3,''};
%form{length(form)+1} = {sepa};
%sepa = {'separation_comm','frame','',-5,''};
%form{length(form)+1} = {sepa};

% 2emes frame : Commandes de data
% -------------------------------
%col1 = {'libel','text','DATA',[],''};
%form{length(form)+1} = {col1};

% ---------------------------------------------
col1 = {'radio_mode'   ,'radio' ,'Transport equations',0 ,'Transport equations parameters'};
form{length(form)+1} = {col1};

%col1 = {'radio_coef'   ,'radio' ,'Transport coefficients' ,0 ,'Edit transport coefficients'};
%form{length(form)+1} = {col1};

%col1 = {'radio_geo'    ,'radio' ,'Geometry'                 ,0 ,'plasma equilibrium geometry'};
%form{length(form)+1} = {col1};

%col1 = {'radio_cons'   ,'radio' ,'Plasma references'       ,0 ,'Define palsma references (heat, injection, ...)'};
%form{length(form)+1} = {col1};

col1 = {'radio_asserv' ,'radio' ,'Feedback control'           ,0 ,'Define feedback control'};
form{length(form)+1} = {col1};

col1 = {'radio_plot'    ,'radio' ,'Graphics'          ,0,'Graphics'};
form{length(form)+1} = {col1};

col1 = {'radio_prof'   ,'radio' ,'Profiles'                   ,0 ,'Edit profiles '};
form{length(form)+1} = {col1};

%col1 = {'radio_source' ,'radio' ,'Sources'                   ,0 ,'Edit sources '};
%form{length(form)+1} = {col1};

%col1 = {'radio_bord'   ,'radio' ,'Edge and wall'            ,0 ,'Edge and wall',[],'','enable','off'};
%form{length(form)+1} = {col1};

%col1 = {'radio_impur'  ,'radio' ,'Impurities and radiation'  ,0 ,'Impurities and radiation',[],'','enable','off'};
%form{length(form)+1} = {col1};

%col1 = {'radio_neo'    ,'radio' ,'Neoclassical data'     ,0 ,'Neoclassical data',[],'','enable','off'};
%form{length(form)+1} = {col1};

%col1 = {'radio_exp'    ,'radio' ,'Experimental data'    ,0 ,'Experimental data',[],'','enable','off'};
%form{length(form)+1} = {col1};

%col1 = {'radio_equi'   ,'radio' ,'Initial equilibrium',0 ,'Compute initial equilibrium'};
%form{length(form)+1} = {col1};

%col1 = {'radio_scal'   ,'radio' ,'Composition Update'         ,0 ,'Update the plasma composition, ion pressure and density when zeff or gaz composition have been modified ',[],''};
%form{length(form)+1} = {col1};

%col1 = {'radio_mhd'    ,'radio' ,'MHD'                       ,0 ,'MHD',[],'','enable','off'};
%form{length(form)+1} = {col1};

%col1 = {'radio_evx'    ,'radio' ,'Events'                ,0 ,'Events triggering',[],'','enable','off'};
%form{length(form)+1} = {col1};

%col1 = {'radio_simulink'    ,'radio' ,'Simulink'                ,0 ,'Open Simulink @ Cronos dialog',[],''};
%form{length(form)+1} = {col1};
col1 = {'radio_postrait'   ,'radio' ,'Post-treatments'                   ,0 ,' '};
form{length(form)+1} = {col1};

col1 = {'radio_device'   ,'radio' ,'Device'                   ,0 ,' '};
form{length(form)+1} = {col1};

% separation
% ----------
%sepa = {'separation_comm','frame',' ',-5,''};
sepa ={'separation_comm','frame','',1,''};
form{length(form)+1} = {sepa};

% 9eme ligne : Bouton Quit
% ------------------------
comm{1}={'btn_quit','radio','Close',0,''};
comm{2}={'aide','radio@right','Help',0,''};

hout=zuicreeform('Edition','edit','zuiedit_fct','',form,comm);

