% ZUICREATE  fonction de creation du formulaire "Creation"
%--------------------------------------------------------------------------
% fichier zuicreate.m ->  
%	zuicreeform	: creation du formulaire
%	zuicreate_fct	: gestion des callbacks des uicontrols
%		   			  du formulaire
%
% fonction Matlab 5 :
%	Cette fonction cree la feuille principale/creation  (interface graphique) sous 
%	la forme d'un formulaire. Le formulaire est constitue de lignes.
%	Chaque ligne peut avoir un nombre de colonnes quelconques. Les lignes
%	qui se suivent et qui ont un nombre de colonnes egal, voient leurs colonnes
%	alignees verticalement.
%
% syntaxe  :
%	zuicreate ;
%
% entrees :
% 
% sorties :
%  				
% fonction ecrite par C. Passeron, poste 61 19
% version 2.0, du 10/12/2002.
% 
% liste des modifications : 
% * 10/12/2002 : interface en anglais
%--------------------------------------------------------------

function zuicreate

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('create');
if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

% 1eres lignes : Rappel du nom du fichier de travail
% --------------------------------------------------
form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel_loadfile','text@full',' Cronos input file creation',[],''};
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% ---------------------------------------------
col1 = {'radio_ts'   ,'radio' ,'Tore Supra',0,'simulations for Tore Supra'};
form{length(form)+1} = {col1};

col1 = {'radio_tsg'   ,'radio' ,'Gcronos  ',0,'simulation for Tore Supra, tcronos option'};
form{length(form)+1} = {col1};

col1 = {'radio_jet'  ,'radio' ,'JET'       ,0,'simulation for JET'};
form{length(form)+1} = {col1};

col1 = {'radio_iter' ,'radio' ,'ITER'      ,0,'simulation for ITER'};
form{length(form)+1} = {col1};

col1 = {'radio_sst1' ,'radio' ,'SST1'      ,0,'simulation for SST1'};
form{length(form)+1} = {col1};

col1 = {'radio_east' ,'radio' ,'EAST'      ,0,'simulation for EAST'};
form{length(form)+1} = {col1};

col1 = {'radio_hl2a' ,'radio' ,'HL2A'      ,0,'simulation for HL2A'};
form{length(form)+1} = {col1};

col1 = {'radio_compass' ,'radio' ,'COMPASS'      ,0,'COMPASS-D simulation'};
form{length(form)+1} = {col1};

col1 = {'radio_ftu' ,'radio' ,'FTU'      ,0,'simulation for FTU'};
form{length(form)+1} = {col1};

col1 = {'radio_tcv' ,'radio' ,'TCV'      ,0,'simulation for TCV'};
form{length(form)+1} = {col1};

col1 = {'radio_diiid' ,'radio' ,'DIIID'      ,0,'simulation for DIIID'};
form{length(form)+1} = {col1};

col1 = {'radio_pdb' ,'radio' ,'Profile DB'      ,0,'simulation for ITPA Profile Database'};
form{length(form)+1} = {col1};

col1 = {'radio_0d' ,'radio' ,'Metis'      ,0,'simulation data set created from Metis simulation'};
form{length(form)+1} = {col1};

col1 = {'radio_autre','radio' ,'Other'     ,0,'Virtual tokamak'};
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

col1 = {'radio_plus','radio' ,'extension'     ,0,'extension of a previous simulation'};
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

col1 = {'radio_tcronos','radio' ,'Tcronos'     ,0,'reading of tcronos file + new simulation'};
form{length(form)+1} = {col1};

% separation
% ----------
sepa = {'separation_comm','frame',' ',-5,''};
form{length(form)+1} = {sepa};

% Bouton Quit
% -----------
comm{1}={'btn_quit','radio@center','close',0,'to close the window'};

hout=zuicreeform('Creation','create','zuicreate_fct','',form,comm);
