% ZUIPARAM formulaire de parametrage du mode profil
%--------------------------------------------------
% fichier zuiparam.m ->  
%		zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	creation du formulaire de parametrage du mode profil
%
% syntaxe  :
%   	hout=zuiparam(nom,text_modul,var_modul) ;
%
% entree :
%	nom        : nom du profil
%	text_modul : texte affiche dans le bouton commandant la modulation temporelle du profil
%	var_modul  : nom de la variable  du workspace servant de consigne de modulation temporelle du profil
%
% sorties :
%   hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 1.6, du 31/08/2001.
% 
% Modifications
%
%--------------------------------------------------------------
function hout=zuiparam(variable,text_modul,var_modul)

if nargin ~= 3
	warning('zuiparam : Nb d''arguments incorrect')
	return
end
if nargin < 2
	text_modul = '';
	var_modul  = '';
end

% cas specifique de Te,Pe, Ti & Pion
switch variable
case 'data.prof.te'
	evalin('base','param.edit.tepe =''te'';');
case 'data.prof.pe'
	evalin('base','param.edit.tepe =''pe'';');
case 'data.prof.ti'
	evalin('base','param.edit.tipion =''ti'';');
case 'data.prof.pion'
	evalin('base','param.edit.tipion =''pion'';');
end	

% formulaire 
% avec la sous structure from
form={};

% Titre
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'titre','text@full',variable,[],''};
form{length(form)+1} = {col1};

% S�aration
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

colj = {'jump','jump',' ',[],''};
col3 = {'list_plot','list','',70,''};
form{length(form)+1} = {colj,colj,colj,col3} ;

% Centre
col1 = {'text_centre','text','center',5,''} ;
col2 = {'edit_centre','edit','1',5,'center value'} ;
col3 = {'jump_void','jump','texte pour reserver de la place encore et encore',[],''};
form{length(form)+1} = {col1,col2,colj,col3} ;

col3 = {'jump_void','jump','texte pour reserver de la place encore et encore',[],''};
form{length(form)+1} = {colj,colj,colj,col3} ;

% Bord
col1 = {'text_bord','text','edge',5,''} ;
col2 = {'edit_bord','edit','0',5,'edge value'} ;
col3 = {'jump_void','jump','texte pour reserver de la place encore et encore',[],''};
form{length(form)+1} = {col1,col2,colj,col3} ;

col3 = {'jump_void','jump','texte pour reserver de la place encore et encore',[],''};
form{length(form)+1} = {colj,colj,colj,col3} ;

% Alpha
col1 = {'text_alpha','text','alpha',5,''} ;
col2 = {'edit_alpha','edit','2',5,'alpha'} ;
colj = {'jump','jump',' ',[],''};
col3 = {'jump_void','jump','texte pour reserver de la place encore et encore',[],''};
form{length(form)+1} = {col1,col2,colj,col3} ;

col3 = {'jump_void','jump','texte pour reserver de la place encore et encore',[],''};
form{length(form)+1} = {colj,colj,colj,col3} ;

% Beta
col1 = {'text_beta','text','beta',5,''} ;
col2 = {'edit_beta','edit','1',5,'beta'} ;
colj = {'jump','jump',' ',[],''};
col3 = {'jump_void','jump','texte pour reserver de la place encore et encore',[],''};
form{length(form)+1} = {col1,col2,colj,col3} ;

col3 = {'jump_void','jump','texte pour reserver de la place encore et encore',[],''};
form{length(form)+1} = {colj,colj,colj,col3} ;
form{length(form)+1} = {colj,colj,colj,col3} ;
% Modulation
col1 = {'text_module','text','modulation',5,''} ;
if ~isempty(text_modul)
    col2 = {'pop_modul' ,'popup',strcat('    none    |',text_modul),1,'temporal dependence of the profiles (by the reference value normalized to one)'};
else
    col2 = {'text_modul','text','none',1,'temporal dependence of the profiles (by the reference value normalized to one)'};
end
form{length(form)+1} = {col1,col2,colj,col3} ;
form{length(form)+1} = {colj,colj,colj,col3} ;
form{length(form)+1} = {colj,colj,colj,col3} ;
form{length(form)+1} = {colj,colj,colj,col3} ;
form{length(form)+1} = {colj,colj,colj,col3} ;
form{length(form)+1} = {colj,colj,colj,col3} ;


% separation
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
sepa = {'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};

% Formulaire
hout=zuicreeform('Profil parameter','parametrer','zuiparam_fct','',form) ;

setappdata(hout,'variable',variable) ;
setappdata(hout,'var_modul',var_modul) ;

%------------------------------------------------------------------------------
% Creation de la partie graphique

[hform,hui] = zuiformhandle('parametrer') ;

% Partie graphique trace du mode
hui.axes_plot=zuiplotin(hui.list_plot) ;

% memorisation des handles
setappdata(hform,'zhandle',hui) ;

x      = evalin('base','param.gene.x') ;

centre = zuidata(hui.edit_centre) ;
bord   = zuidata(hui.edit_bord) ;
alpha  = zuidata(hui.edit_alpha) ;
beta   = zuidata(hui.edit_beta) ;

% calcul profil = (centre-bord)*(1.-rho^alpha)^beta + bord
profil = (centre - bord) * (1 - x.^alpha).^beta + bord ;

% trace du profil
axes(hui.axes_plot)
plot(x,profil)
title('(center-edge)*(1-x^\alpha)^\beta + edge')
xlabel('x (su)');
zoom on

