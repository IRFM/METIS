% ZUIEDIT_PARAM_COMP   formulaire 'composition' sous Edition
%--------------------------------------------------------------
% fichier zuiedit_param_comp.m ->
%		zuicreeform	: creation du formulaire
%
% fonction Matlab 5 :
%	fonction de creation de GUI pour le formulaire
%	de composition sous commandes des parametres
% 	du mode edition du formulaire principal
%
% syntaxe  :
%  zuiedit_param_comp ;
%
% entrees :
%
% sorties :
%	hform : handle du formulaire cree
%
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 30/08/2001.
%
% liste des modifications
%  * 30/08/2001 -> ajout du handle de sortie (J-F Artaud)
%  * 30/08/2001 -> reinitialisation si deja ouvert (J-F Artaud)
%
%--------------------------------------------------------------

function hform = zuiedit_param_comp

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('ed_param_comp');

if ishandle(hform)
	zuiformvisible(hform) ;
	zuiedit_param_comp_fct('init');
	return
end


%--------------------------------------------------------------------
% le formulaire
% nb d'impuretes
nbg = evalin('base','param.gene.nbg') ;

form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel_loadfile','text@full','Plasma composition ',[],''};
form{length(form)+1} = {col1};

col1 = {'libel_loadfile','text@full','ions species ',[],''};
form{length(form)+1} = {col1};

sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% Gaz
% ---
colj = {'void','jump','void',1,''} ;
col1 = {'majoritaire','text','main',20,'main gas name',''} ;
col2 = {'minoritaire1','text','1st minority',20,'first minority species',''} ;
col3 = {'minoritaire2','text','2nd minority',20,'second minority species',''} ;
form{length(form)+1} = {col1,colj,col2,colj,col3} ;

liste = 'H|D|T|He' ;
ind  = evalin('base','param.compo.a(1)') ;
col1 = {'pop_majoritaire','popup',liste,ind,'main species',liste} ;

liste = strcat(liste,'|He3|Li') ;
ind1  = evalin('base','param.compo.z(2)') ;
ind   = evalin('base','param.compo.a(2)') ;
if ind1==2
	if ind==4
		ind = 4 ;
	elseif ind==3
		ind = 5 ;
	elseif ind==7
		ind = 6 ;
	end
end
col2 = {'pop_minoritaire1','popup',liste,ind,'first minority species',liste} ;

ind1 = evalin('base','param.compo.z(3)') ;
ind  = evalin('base','param.compo.a(3)') ;
if ind1==2
	if ind==4
		ind = 4 ;
	elseif ind==3
		ind = 5 ;
	elseif ind==7
		ind = 6 ;
	end
end
col3 = {'pop_minoritaire2','popup',liste,ind,'second minority species',liste} ;
form{length(form)+1} = {col1,colj,col2,colj,col3};

% separation
% ----------
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
% separation
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
% separation
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% Impuretes
% ---------
col1 = {'impuretes','text','impureties',12,'impureties',''} ;
col2=  {'jump1','jump','jump',[],''};
cols={};
for l=4:nbg
	aide   = sprintf('impureties# %d',l);
	col2={strcat('colonne',int2str(l)),'text',sprintf('%5d',l),1,aide,[]};
	cols{l-3}= col2;
end
form{length(form)+1} = {col1,cols{:}};

% Charge
% ------
col1 = {'charge','text','charge',12,'charge',''} ;
col2=  {'jump1','jump','jump',[],''};
cols={};
for l=4:nbg
	aide   = sprintf('impurity charge # %d',l);
	var = strcat('param.compo.z(',int2str(l),')') ;
	chg = evalin('base',var) ;
	col2={strcat('charge',int2str(l)),'edit',sprintf('%5d',chg),1,aide,'',var};
	chg = '' ;
	cols{l-3}= col2;
end
form{length(form)+1} = {col1,cols{:}};

% Masse
% -----
col1 = {'masse','text','mass',12,'mass number',''} ;
col2=  {'jump1','jump','jump',[],''};
cols={};
for l=4:nbg
	aide   = sprintf('impurity mass# %d',l);
	var = strcat('param.compo.a','(',int2str(l),')') ;
	ms  = evalin('base',var) ;
	col2 = {strcat('masse',int2str(l)),'edit',sprintf('%5d',ms),1,aide,'',var};
	ms = '' ;
	cols{l-3}= col2;
end
form{length(form)+1} = {col1,cols{:}};

% separation
% ----------
sepa = {'separation_comm','frame','',1,''};
form{length(form)+1} = {sepa};
colj = {'jump','jump',' ',[],''};
colj = {'jump','jump',' ',[],''};
form{length(form)+1} = {colj,colj};

hform = zuicreeform('Edition','ed_param_comp','zuiedit_param_comp_fct','',form);

