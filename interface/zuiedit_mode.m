% ZUIEDIT_MODE   formulaire de modification de la valeur des modes
%--------------------------------------------------------------
% fichier zuiedit_mode.m ->  zuiedit_mode
%		                       zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	creation du formulaire de modification de la valeur des modes 
%
% syntaxe  :
%	hout=zuiedit_mode(nom_mode,valeur_mode,type_mode);
%
% entree :
%	nom_mode : nom de la variable de mode 
%	valeur_mode : valeurs possibles {0,1,2...}
%	type_mode : signification {mis a zero,lu en entree, calcule ...}
%
% sorties :
%	hout : handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 1.3, du 26/06/2001.
% 
% liste des modifications : 
%
%  * 14/05/2002 -> correction redondance tag
%
%--------------------------------------------------------------
function hout=zuiedit_mode(nom_mode,valeur_mode,type_mode)

% argument

if nargin < 1
	error('il faut donner le nom de la variable mode, la liste des valeurs et leur signification ') ;
end
if nargin < 2
	error ('il faut donner la liste des valeurs et leur signification') ;
end
if nargin < 3
	error ('il faut donner la signification') ;
end

% nom du tag
% ----------
%[tagc,reste] = strtok(nom_mode,'.') ;
%while ~isempty(reste)
%	[tagc,reste] = strtok(reste,'.') ;
%end
tagc  = strrep(nom_mode,'data.','') ;
tagc  = strrep(tagc,'.','_') ;
% on supprime les caracteres : , ( ) du tagc
tagc=strrep(tagc,'(','') ;
tagc=strrep(tagc,':','') ;
tagc=strrep(tagc,',','') ;
tagc=strrep(tagc,')','') ;

% si l'interface a deja ete appelee
[hform,hui] = zuiformhandle([tagc '_mode']) ;
if ishandle(hform)
%          zuiformvisible(hform) ;
	zuicloseone(hform) ;
	zuiedit_mode(nom_mode,valeur_mode,type_mode) ;
	return
end

% formulaire d'edition de modes
% avec la sous structure from
form={};

% Titre
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

colj = {'jump','jump',' ',[],''};
col1 = {'nom_mode','text@full',nom_mode,[],''};
form{length(form)+1} = {colj,col1};

% S�aration
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% Liste des intervalles : nom [xdeb - xfin]
col1 = {'list_intervalle','list','texte pour reserver de la place',1,''};
coljg = {'jump_void','jump','',20,''};
col3 = {'list_plot','list','',70,''};
form{length(form)+1} = {colj,col1,coljg,col3,colj};

% Partie graphique
col1 = {'jump_void','jump','texte pour reserver de la place encore et encore et encore',[],''};
col3 = {'jump_void','jump','',70,''};
form{length(form)+1} = {colj,col1,coljg,col3,colj};
form{length(form)+1} = {colj,col1,coljg,col3,colj};
form{length(form)+1} = {colj,col1,coljg,col3,colj};
form{length(form)+1} = {colj,col1,coljg,col3,colj};
form{length(form)+1} = {colj,col1,coljg,col3,colj};
form{length(form)+1} = {colj,col1,coljg,col3,colj};
form{length(form)+1} = {colj,col1,coljg,col3,colj};
form{length(form)+1} = {colj,colj,coljg,col3,colj};
form{length(form)+1} = {colj,colj,coljg,col3,colj};
form{length(form)+1} = {colj,colj,coljg,col3,colj};
form{length(form)+1} = {colj,colj,coljg,col3,colj};
form{length(form)+1} = {colj,colj,coljg,col3,colj};

% Separation
sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

% Radio boutons : modification des intervalles, Importer, Xgrid
col1 = {'radio_modif','radio','interval edition',0,'Open form of interval edition',[],''} ;
col2 = {'radio_import','radio','Import',0,'Copy values from another mode',[],''} ;
col3 = {'radio_xgrid','radio','Xgrid',0,'turn on x grid',[],''} ;
colj  = {'void','jump','v',1,''} ;
form{length(form)+1} = {colj,col1,colj,colj,colj,coljg,col2,colj,col3,colj,colj,colj,colj,colj};

% saut de ligne
form{length(form)+1} = {colj,colj,colj,colj,colj,coljg,colj,colj,colj,colj,colj,colj,colj,colj};

% Legende des boutons de controles
col1a = {'text_tps'     ,'text@merge'  ,'time',10,''} ;
col1b = {'text_tps2'    ,'text@merge' ,'     ',10,''} ;
col1c = {'text_tps3'    ,'text@merge' ,'inter',10,''} ;
col1d = {'text_tps4'    ,'text'       ,'val ',10,''} ;
col2  = {'text_type'   ,'text' ,'type'               ,5,'Values allowed for the mode',''} ;
col3  = {'text_indicie','text' ,'Time or Index'    ,10,'Time step or index at which the mode value is changed',''} ;
col4  = {'text_valeur' ,'text' ,'value'             ,30,'mode value',''} ;
form{length(form)+1} = {colj,col1a,col1b,col1c,col1d,coljg,col2,colj,col3,colj,colj,col4,colj,colj} ;

% Valeurs de l'intervalle et Boutons de controles du mode
x = evalin('base','param.intervalle.calcul') ;
x = x(:) ;

col1a = {'text_tmin'  ,'text','tmin',5,''} ;
col1b = {'valeur_tmin','text',sprintf('%10s',num2str(x(1))),10,''} ;
col1c = {'text_tmax'  ,'text','tmax',5,''} ;
col1d = {'valeur_tmax','text',sprintf('%10s',num2str(x(2))),10,''} ;
type_popup = {'    at the index     ', ...
              '    at the time      ', ...
              '    index step       ', ...
              '    time step        ', ...
              '    all the interval ' } ;
valeur_popup = [1,2,3,4,5] ;
col2 = {'pop_1','popup',type_popup,1,'Values allowed for the mode',valeur_popup,''} ;
col3 = {'edit_1','edit','     ',10,'Time step or index at which the mode value is changed'} ;
col3a = {'text_com','text','',5,'',[],'','visible','off'} ;
col4 = {'pop_2','popup',type_mode,1,'mode type',valeur_mode,''} ;
col5 = {'go','radio','Do',0,'applied the modification',[],''} ;
form{length(form)+1} = {colj,col1a,col1b,col1c,col1d,coljg,col2,colj,col3,col3a,colj,col4,colj,col5} ; 

% recherche des indices kmin et kmax 
tt = evalin('base','data.gene.temps') ;
dt = abs(tt-x(1)) ;
kmin = min(find(dt==min(dt))) ;
dt = abs(tt-x(2)) ;
kmax = min(find(dt==min(dt))) ;

% affichage kmin, kmax
col1a = {'text_kmin'  ,'text','kmin',5,''} ;
col1b = {'valeur_kmin','text',sprintf('%10s',num2str(kmin)),10,''} ;
col1c = {'text_kmax'  ,'text','kmax',5,''} ;
col1d = {'valeur_kmax','text',sprintf('%10s',num2str(kmax)),10,''} ;
col2 = {'retour','radio','undo',0,'Canceled the last action',[],''} ;
form{length(form)+1} = {colj,col1a,col1b,col1c,col1d,coljg,col2,colj,colj,colj,colj,colj,colj,colj} ;

% saut de ligne
form{length(form)+1} = {colj} ;

% Separation
sepa ={'separation_comm','frame','',3,''} ;
form{length(form)+1} = {sepa} ;

% Ligne de commentaires
col3 = {'commentaire','text','ceci est une ligne de commentaires pour les messages d''erreur ',10,'',[],'','visible','off'} ;
form{length(form)+1} = {colj,colj,col3} ;

% Separation
sepa ={'separation_comm','frame','',3,''} ;
form{length(form)+1} = {sepa} ;

% Formulaire
hout=zuicreeform('Mode',[tagc,'_mode'],'zuiedit_mode_fct','',form) ;

%------------------------------------------------------------------------------
% Creation de la partie graphique

[hform,hui] = zuiformhandle([tagc,'_mode']) ;

% Trace du mode
hui.axes_plot=zuiplotin(hui.list_plot) ;
title('') ;
xlabel('temps') ;
ylabel('') ;
xlim = [tt(1) tt(end)] ;
val  = cat(1,valeur_mode{:}) ;
ylim = [val(1)-1 val(end)+1] ;

set(hui.axes_plot,'xlim',xlim,'ylim',ylim) ;

%ytick = get(hui.axes_plot,'ytick')
%YTickLabel{1} = ' ' ;
%for i=1:length(valeur_mode)
%	ind = find(ytick==valeur_mode{i}) ;
%	YTickLabel{ind} = type_mode{i}
%end
%set(hui.axes_plot,'YTickLabel',YTickLabel) ;

set(hui.axes_plot,'xlim',xlim,'ylim',ylim,'YTickLabel','','YTick',[]) ;
zoom on

% om memorise le mode initial
y = evalin('base',nom_mode) ;

mode{1} = y ;

setappdata(hform,'xdata',tt) ;
setappdata(hform,'ydata',mode) ;
setappdata(hform,'valeur',valeur_mode) ;

color = get(0,'defaultaxesColorOrder') ;
color_legende = color(1:length(type_mode),:) ;

for i=1:length(valeur_mode)
	cmd = ['hui.line_' num2str(i) ' = line(''xdata'',tt,''ydata'',[valeur_mode{' num2str(i) '} '] ;
%	cmd = [cmd 'valeur_mode{' num2str(i) '}],''color'',color_legende(' num2str(i) ',:));' ] ;
	cmd = [cmd 'valeur_mode{' num2str(i) '}],''color'',color_legende(' num2str(i) ',:),' ] ;
	cmd = [cmd '''linestyle'',''none'',''visible'',''on'') ;'] ;
	eval(cmd) ;
end

% ecriture de la legende
xtick = get(hui.axes_plot,'Xtick') ;

for i=1:length(type_mode)
% 	x_legende = -0.10*(xlim(2)-xlim(1))
 	x_legende = xlim(1)-0.65*(xtick(2)-xtick(1)) ;
	y_legende =  val(i) ;
	hui.legende = text(x_legende,y_legende,type_mode(i),'color',color_legende(i,:), ...
                           'units','data','hori','center') ;
	units = get(hui.legende,'units') ;
	set(hui.legende,'units','normalized');
	pos=get(hui.legende,'position');
	posy = pos(2);
	posx = -0.22;
	set(hui.legende,'position',[posx posy pos(3)])
	set(hui.legende,'units',units)
end

% Trace de l'intervalle
hui.trait_intrvl = line([0 0],[0 0],'linestyle','-','linewidth',2);
set(hui.trait_intrvl,'xdata',[tt(1) tt(end)], ...
             'ydata',[ylim(2)-0.25 ylim(2)-0.25], ...
             'linewidth',3,'visible','on','linestyle','-', ...
             'color',[0 0 0]) ;

% liste des intervalles
int = evalin('base','param.intervalle') ;
name_int = fieldnames(int) ;
xd = [] ;
xf = [] ;
name = [] ;
for i=1:length(name_int)
	val  = getfield(int,name_int{i}) ;
 	xd   = [xd val(1)] ;
	xf   = [xf val(2)] ;
	name = strvcat(name,name_int{i}) ;
end
setappdata(hui.list_intervalle,'libel',name) ;
setappdata(hui.list_intervalle,'xd',xd) ;
setappdata(hui.list_intervalle,'xf',xf) ;

% Affichage liste des intervalles
zuilisteintrvl([tagc,'_mode']') ;

% memorisation des handles
setappdata(hform,'zhandle',hui) ;

% trace du mode initial
zuiedit_mode_plot(tt,y) ;

