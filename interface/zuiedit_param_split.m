% ZUIEDIT_PARAM_SPLIT   formulaire 'decoupe temporelle' sous 'Edition'
%--------------------------------------------------------------
% fichier zuiedit_param_split.m ->  
%		zuicreeform : creation du formulaire
%		zuiedit_param_plot_fct	: fonction de control des callbacks
%		zuiedit_param_plot_crtl	: fonction de control des donnees
% 
% fonction Matlab 5 :
%	formulaire 'decoupe temporelle' sous 'Edition'
%	fonction de creation de GUI pour le formulaire
%	de decoupe temporelle (split), commandes de parametres  
%	sous le mode edition du formulaire principal%
%
% syntaxe  :
%   	zuiedit_parma_split ;
%
% entrees :
%
% sorties :
%
% fonction ecrite par C. Passeron, poste 61 19
% version 1.3, du 10/04/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------

function zuiedit_param_split

[hfig,h] = zuiformhandle('ed_param_split');
if ishandle(hfig)
	zuiformvisible(hfig) ;
	return
end

info=zinfo ;

%--------------------------------------------------------------------
% le formulaire

form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel_loadfile','text@full','Global Parameters ',[],''};
form{length(form)+1} = {col1};

col1 = {'libel_loadfile','text@full','time splitting',[],''};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% Informations : param.plot
% -------------------------
colj = {'void','jump','void',[],''} ;
form{length(form)+1} = {colj,colj,colj,colj} ;

col1 = {'libel','text','Split',5,info.param.split.onoff,''};
val = evalin('base','param.split.onoff') ;
col2 = {'onoff_split','popup',{'Off','On'},val+1,info.param.split.onoff,{0,1},'param.split.onoff'};
form{length(form)+1} = {colj,col1,col2,colj};

col1 = {'libel','text','split mode',5,info.param.split.mode,''};
val_mode = evalin('base','param.split.mode')  ;
col2 = {'mode_split','popup',{'forced','automatic'},val_mode+1,info.param.split.mode,{0,1},'param.split.mode'};
form{length(form)+1} = {colj,col1,col2,colj};

col1 = {'libel','text','dtmax',5,info.param.split.dtmax,''};
val = evalin('base','param.split.dtmax') ;
col2 = {'dtmax_split','edit',val,5,info.param.split.dtmax,[],'param.split.dtmax'};
%col2 = {'dtmax_split','edit',val,5,info.param.split.dtmax,[],''};
col3 = {'borne','text','[1.e-5  0.1]',5,info.param.split.dtmax,[],''};
form{length(form)+1} = {colj,col1,col2,col3};

col1 = {'libel','text','dtmin',5,info.param.split.dtmin,''};
val = evalin('base','param.split.dtmin') ;
col2 = {'dtmin_split','edit',val,5,info.param.split.dtmin,[],'param.split.dtmin'};
col3 = {'borne','text','[1.e-6  0.01]',5,info.param.split.dtmin,[],''};
form{length(form)+1} = {colj,col1,col2,col3};

val = evalin('base','param.split.nb') ;
if val_mode==0
	col1 = {'libel_nb','text','nb',5,info.param.split.nb,''};
	col2 = {'nb_split','edit',val,5,info.param.split.nb,[],'param.split.nb'};
	col3 = {'borne_nb','text','>1',5,info.param.split.nb,[],''};
else
 	col1 = {'libel_nb','text','nb',5,info.param.split.nb,[],'','Enable','off'};
 	col2 = {'nb_split','edit',val,5,info.param.split.nb,[],'param.split.nb','Enable','off'};
 	col3 = {'borne_nb','text','>1',5,info.param.split.nb,[],'','Enable','off'};
end
form{length(form)+1} = {colj,col1,col2,col3};

col1 = {'libel','text','equi',5,info.param.split.equi,''};
val = evalin('base','param.split.equi') ;
col2 = {'equi_split','edit',val,5,info.param.split.equi,[],'param.split.equi'};
col3 = {'borne','text','[dtmin dtmax]',5,info.param.split.equi,[],''};
form{length(form)+1} = {colj,col1,col2,col3};

hout=zuicreeform('Edition','ed_param_split','zuiedit_param_split_fct','zuiedit_param_split_ctrl',form);
