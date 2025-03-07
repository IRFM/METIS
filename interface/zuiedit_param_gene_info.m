% ZUIEDIT_PARAM_GENE_INFO formulaire 'informations' sous 'Parametres generaux - Edition'
%--------------------------------------------------------------
%
% fichier zuiedit_param_gene_info.m ->
%			zuicreeform 	: creation du formulaire
%			zuiedit_param_gene_info_fct : fonction de gestion des callbacks des uicontrols
%                                         du formulaire

% fonction Matlab 5 : 
%	fonction de creation de GUI pour le formulaire
%	d'informations des parametres sous le mode edition du formulaire principal
%
% syntaxe  :
%  
%  zuiedit_param_gene_info ;
%
% entrees
%
% sorties :
%
% fonction ecrite par C. Passeron, poste 61 19
% version 1.3, du 10/04/2001.
%
% 
% liste des modifications : 
%
%--------------------------------------------------------------

function zuiedit_param_gene_info

% si le formulaire existe deja, on l'active
[hfig,h] = zuiformhandle('ed_param_gene_info') ;
if ishandle(hfig)
	zuiformvisible(hfig)
	return
end

info=zinfo;
form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel_loadfile','text@full','Global simulation parameters',[],''};
form{length(form)+1} = {col1};

colj = {'jump','jump',' ',[],''};
col1 = {'libel_loadfile','text@full','Information',[],''};
form{length(form)+1} = {colj,col1};

%  
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% le formulaire
%--------------------------------------------------------------------
val  = evalin('base','param.gene.nbrho');
col1 = {'libel','text',sprintf('%-15s','nbrho'),[],info.param.gene.nbrho,[],[]};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.nbt');
col1 = {'libel','text',sprintf('%-15s','nbt'),[],info.param.gene.nbt,[],[]};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.nbeq');
col1 = {'libel','text',sprintf('%-15s','nbeq'),[],info.param.gene.nbeq,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.dt');
col1 = {'libel','text',sprintf('%-15s','dt'),[],info.param.gene.dt,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.dk');
col1 = {'libel','text',sprintf('%-15s','dk'),[],info.param.gene.dk,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.nbfci');
col1 = {'libel','text',sprintf('%-15s','nbfci'),[],info.param.gene.nbfci,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.nbfce');
col1 = {'libel','text',sprintf('%-15s','nbfce'),[],info.param.gene.nbfce,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.nbhyb');
col1 = {'libel','text',sprintf('%-15s','nbhyb'),[],info.param.gene.nbhyb,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.nbidn');
col1 = {'libel','text',sprintf('%-15s','nbidn'),[],info.param.gene.nbidn,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.nbglacon');
col1 = {'libel','text',sprintf('%-15s','nbglacon'),[],info.param.gene.nbglacon,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.dx');
col1 = {'libel','text',sprintf('%-15s','dx'),[],info.param.gene.dx,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.nbg');
col1 = {'libel','text',sprintf('%-15s','nbg'),[],info.param.gene.nbg,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.nbrhorz');
col1 = {'libel','text',sprintf('%-15s','nbrhorz'),[],info.param.gene.nbrhorz,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.nbthetarz');
col1 = {'libel','text',sprintf('%-15s','nbthetarz'),[],info.param.gene.nbthetarz,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.nbsepa');
col1 = {'libel','text',sprintf('%-15s','nbsepa'),[],info.param.gene.nbsepa,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.origine');
ind  = max(find(val == '/'));
val  = strcat('..',val(ind:length(val)));
col1 = {'libel','text',sprintf('%-15s','file'),[],info.param.gene.origine,''};
col2 = {'libel','text',sprintf('%-10s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.version_zineb');
col1 = {'libel','text',sprintf('%-15s','version_zineb'),[],info.param.gene.version_zineb,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

val  = evalin('base','param.gene.date_zineb');
col1 = {'libel','text',sprintf('%-15s','date'),[],info.param.gene.date_zineb,''};
col2 = {'libel','text',sprintf('%-30s',num2str(val)),[],''};
form{length(form)+1} = {col1,col2};

% separation
% ----------
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% Bouton Quit
% ----------
comm{1}={'btn_quit','radio@center','Close',0,''};

hout=zuicreeform('Edition','ed_param_gene_info','zuiedit_param_gene_info_fct','',form,comm);
