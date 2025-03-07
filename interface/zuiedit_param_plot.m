% ZUIEDIT_PARAM_PLOT formulaire 'graphique' sous 'Edition'
%--------------------------------------------------------------
% fichier zuiedit_param_plot.m ->  
%	zuicreeform : creation du formulaire
%	zuiedit_param_plot_fct	: fonction de control des callbacks
% 
% fonction Matlab 5 :
%	creation du formulaire 'graphique' sous 'Edition'
%
% syntaxe  :
%	zuiedit_parma_plot;
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

function zuiedit_param_plot

[hfig,h] = zuiformhandle('ed_param_plot');
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

col1 = {'libel_loadfile','text@full','Graphics ',[],''};
form{length(form)+1} = {col1};

col1 = {'libel_loadfile','text@full','plot',[],''};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% Informations : param.plot
% -------------------------
colj = {'void','jump','void',[],''} ;

col1 = {'libel','text','Plot',10,info.param.plot.onoff,''};
val = evalin('base','param.plot.onoff') ;
col2 = {'onoff_plot','popup',{'off','on'},val+1,info.param.plot.onoff,{0,1},'param.plot.onoff'};
form{length(form)+1} = {colj,col1,col2};

col1 = {'libel','text','interval',10,info.param.plot.intervalle,''};
val = evalin('base','param.plot.intervalle') ;
tps = {'    none      ', ...
       '  any time    ', ...
       '     1/2      ', ...
       '     1/5      ', ...
       '    1/10      '} ;
userdata = [0,1,2,5,10] ;
ind = findstr(userdata,val) ;
col2 = {'intervalle_plot','popup',tps,ind,info.param.plot.intervalle,userdata,'param.plot.intervalle'};
form{length(form)+1} = {colj,col1,col2};

col1 = {'libel','text','start running at launch',10,'',''};
val = evalin('base','param.plot.run') ;
col2 = {'run_plot','popup',{'no','yes'},val+1,info.param.plot.run,{0,1},'param.plot.run'};
form{length(form)+1} = {colj,col1,col2};

col1 = {'libel','text','sleep',10,'',''};
val = evalin('base','param.plot.pause') ;
col2 = {'run_plot','edit',val,5,info.param.plot.pause,[],'param.plot.pause'};
form{length(form)+1} = {colj,col1,col2};

hout=zuicreeform('Edition','ed_param_plot','zuiedit_param_plot_fct','',form);
