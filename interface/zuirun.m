% ZUIRUN   formulaire de soumission en interactif de zineb
%--------------------------------------------------------------
% fichier zuirun.m ->  zuicreeform, gmarge
%
% fonction Matlab 5 :
%	Cette fonction cree la feuille (interface graphique) sous 
%	la forme d'un formulaire. Le formulaire est constitue de lignes.
%	Chaque ligne peut avoir un nombre de colonnes quelconques. Les lignes
%	qui se suivent et qui ont un nombre de colonnes egal, voient leurs colonnes
%	alignees verticalement.
%
% syntaxe  : 
%	zuirun ;
%
% entrees :
%
% sorties :
%
% fonction �rite par J-F Artaud , poste 46-78
% version  1.9  du  18/03/2002 
% 
% liste des modifications : 
% * 28/11/2001 -> modification de la presentation
% * 18/03/2002 -> ajout du mode reprise sans reinitilisation de l'equilibre et des sources
%
%--------------------------------------------------------------
%
function zuirun

% initialisation de cr
cr    = 0;

filename = evalin('base','param.edit.currentfile','[]');
if isempty(filename)
	disp(' no input data file')
	warndlg(' no input data file','Zineb run interactif');
	return
end

% formulaire 
% ----------
form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

% 1ere ligne : fichier de simulation
% ---------------------------------- 
col1 = {'libel','text@full',' input data file name',[],''};
form{length(form)+1} = {col1};

col1 = {'text','text@full',filename,35,'input data file name',[],[]};
form{length(form)+1} = {col1};

% separation
% ----------
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% 2eme ligne : reprise ?
% ----------------------
filerep  = evalin('base','param.gene.file','[]');
%
colj=  {'jump1','jump','jump',[],''};
colm1 = {'text_reprise_run','text','Cronos run',[],''};
if ~isempty(filerep)
	col1 = {'radio_reprise_run','popup','Complete          |restart (init)| restart (without init)' ,1,' restart of a partial Cronos run with/without initialization of the equilibrium and source term',{0,1,2}};
else
	col1 = {'radio_reprise_run','jump',' ' ,0,''};
end
form{length(form)+1} = {colm1,colj,col1};

% separation
% ----------
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

hout=zuicreeform('Interface run Zineb','run','zuirun_fonction',' ',form);
