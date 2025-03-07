% ZUIEDIT_PARAM_INTERVALLE   formulaire 'intervalle' sous 'Edition'
%--------------------------------------------------------------
% fichier zuiedit_param_intervalle.m ->  
%	zuicreeform : creation du formulaire
% 
% fonction Matlab 5 :
%	fonction de creation de GUI pour le formulaire
%	d'intervalles , commandes de parametres  
% 	sous le mode edition du formulaire principal
%
% syntaxe  :
%	zuiedit_param_intervalle;
%
% entrees
%
% sorties :
%	hout   ; handle du formulaire
%
% fonction ecrite par C. Passeron, poste 61 19
% version 1.7, du 08/10/2001.
% 
% liste des modifications : 
%	08/10/2001 -> rajout du zoom
%
%--------------------------------------------------------------
function hout=zuiedit_param_intervalle

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle('intervalle');

if ishandle(hform)
        zuiformvisible(hform) ;
	return
end

liste_ref1 = '      Ip    |Vloop|Flux|ICRH|ECRH|LH|NBI|Gas|Pump|Zeffm|Pellet|empty';
var_ref1   = {{'data.gene.temps','data.cons.ip','-'}, ...
             {'data.gene.temps','data.cons.vloop','-'}, ...
             {'data.gene.temps','data.cons.flux','-'}, ...
             {'data.gene.temps','abs(data.cons.fci)','-'}, ...
             {'data.gene.temps','abs(data.cons.fce)','-'}, ...
             {'data.gene.temps','abs(data.cons.hyb)','-'}, ...
             {'data.gene.temps','abs(data.cons.idn)','-'}, ...
             {'data.gene.temps','data.cons.c','-'}, ...
             {'data.gene.temps','data.cons.pomp','-'}, ...
             {'data.gene.temps','data.cons.zeffm','-'}, ...
             {'data.gene.temps','data.cons.glacon','-'}, ...
             {'[]','[]',''}};

liste_ref2 = '    nl0    |ne0|ne1|nemoy|te0|te1|li|q0|c|empty';
var_ref2   = {{'data.gene.temps','data.cons.asser.nl0',':'}, ...
             {'data.gene.temps','data.cons.asser.ne0',':'}, ...
             {'data.gene.temps','data.cons.asser.ne1',':'}, ...
             {'data.gene.temps','abs(data.cons.asser.nemoy)',':'}, ...
             {'data.gene.temps','abs(data.cons.asser.te0)',':'}, ...
             {'data.gene.temps','abs(data.cons.asser.te1)',':'}, ...
             {'data.gene.temps','abs(data.cons.asser.li)',':'}, ...
             {'data.gene.temps','data.cons.asser.q0',':'}, ...
             {'data.gene.temps','data.cons.asser.c',':'}, ...
             {'[]','[]',''}};


x= evalin('base','param.intervalle.calcul') ;
y=[0,0];

% orientation
%if size(x,1) > 1
%	spline_flag =0;
%else
%	spline_flag =1;
%end
x=x(:);
y=y(:);

% choix du multiplicateur
%yy    = abs(y);
%ind   = find(yy > 0);
%multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
%if isempty(multi)
%    multi = 1;
%elseif ~isfinite(multi)
%    multi = 1;
%end

% extraction des points
% x_mem = x;
% y_mem = y;
% curve_init =0;
% if length(x) <31
%     % on garde tout les points
%     % mode lineaire par defaut
% elseif all(~isfinite(y))
%     x = linspace(min(x),max(x),11);
%     y = zeros(size(x));
%     multi = 1;
% else
%
%    % choix des noeux
%    [x,y]= zextraitnoeud(x,y,11 + 20 .* (1 - spline_flag),spline_flag);
%    curve_init = 1;
% end

% Color
color = get(0,'defaultaxesColorOrder') ;

% formulaire de description des intervalles
% avec la sous structure from
form={};
colj = {'jump_void','jump',' ',[],''};

% Titre
col1 = {'titre','text@full',' Intervals ',[],''};
form{1} = {col1};

% S�aration
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% Popup : consignes et asservissements
col1 = {'text','text','   References  ',10,''};
col2 = {'text','text','feedback references',10,''};
col3 = {'list_plot','list','Ceci est un texte tres long pour reserver de la placeplaceplaceplaceplaceplaceplaceplaceplace',1,''};
form{length(form)+1} = {col1,col2,col3};

col1 = {'cons_1' ,'popup',liste_ref1,length(var_ref1),'trace la reference 5',var_ref1,'','BackgroundColor',color(1,:)};
col2 = {'asser_1','popup',liste_ref2,length(var_ref2),'trace la reference 1',var_ref2,'','BackgroundColor',color(1,:)};
col3 = {'jump_void','jump','Ceci est un texte tres long pour reserver de la place',[],''};
form{length(form)+1} = {col1,col2,col3};

col1 = {'cons_2' ,'popup',liste_ref1,length(var_ref1),'trace la reference 1',var_ref1,'','BackgroundColor',color(2,:)};
col2 = {'asser_2','popup',liste_ref2,length(var_ref2),'trace la reference 1',var_ref2,'','BackgroundColor',color(2,:)};
form{length(form)+1} = {col1,col2,col3};

col1 = {'cons_3' ,'popup',liste_ref1,length(var_ref1),'trace la reference 1',var_ref1,'','BackgroundColor',color(3,:)};
col2 = {'asser_3','popup',liste_ref2,length(var_ref2),'trace la reference 1',var_ref2,'','BackgroundColor',color(3,:)};
form{length(form)+1} = {col1,col2,col3};

col1 = {'cons_4' ,'popup',liste_ref1,length(var_ref1),'trace la reference 1',var_ref1,'','BackgroundColor',color(4,:)};
col2 = {'asser_4','popup',liste_ref2,length(var_ref2),'trace la reference 1',var_ref2,'','BackgroundColor',color(4,:)};
form{length(form)+1} = {col1,col2,col3};

col1 = {'cons_5' ,'popup',liste_ref1,length(var_ref1),'trace la reference 1',var_ref1,'','BackgroundColor',color(5,:)};
col2 = {'asser_5','popup',liste_ref2,length(var_ref2),'trace la reference 1',var_ref2,'','BackgroundColor',color(5,:)};
form{length(form)+1} = {col1,col2,col3};

form{length(form)+1} = {colj,colj,col3};
form{length(form)+1} = {colj,colj,col3};
form{length(form)+1} = {colj,colj,col3};
form{length(form)+1} = {colj,colj,col3};
form{length(form)+1} = {colj,colj,col3};
form{length(form)+1} = {colj,colj,col3};
form{length(form)+1} = {colj,colj,col3};
form{length(form)+1} = {colj,colj,col3};
form{length(form)+1} = {colj,colj,col3};

% Separation
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};


% Liste des intervalles : nom [xdeb - xfin]
%col1 = {'list_intervalles','list','Ceci est un texte pour reserver de la place',1,''};
col1 = {'list_intervalle','list','                                            ',1,''};
colj = {'jump_void','jump',' ',[],''};
coljg = {'jump_void','jump','ceci est un grand jump',[],''};
form{length(form)+1} = {col1,colj,coljg,coljg,colj};

% Bouton edit des saisie du temps de debut et de fin + popup enchaine 
col1 = {'jump_void','jump','Ceci est un texte pour reserver de la place',[],''};
col2 = {'text_X','text','initial time',5,'initial time of the interval'};
col3 = {'text_Y','text','final time  ',5,'final time of the interval'};
form{length(form)+1} = {col1,colj,col2,col3,colj};

tdeb_calcul = evalin('base','param.intervalle.calcul(1)') ;
col2 = {'valeur_xdeb','edit',tdeb_calcul,5,'initial time of the interval'};
tfin_calcul = evalin('base','param.intervalle.calcul(2)') ;
col3 = {'valeur_xfin','edit',tfin_calcul,5,'final time of the interval'};
col4 = {'radio_enchaine','radio','follow',0,'Intial time = final time of the previous interval',[],''};
form{length(form)+1} = {col1,colj,col2,col3,col4};

% popup de saisie de xdeb et xfin avec ginput sur le graphique
col2 = {'select_xdeb','radio','select',0,'Choose on the graph the initial time'};
col3 = {'select_xfin','radio','select',0,'Choose on the graph the final time'};
form{length(form)+1} = {col1,colj,col2,col3,colj};

form{length(form)+1} = {col1,colj,coljg,coljg,colj};

% Saisie du mon de l'intervalle_ + popup Ajout et Supprime
col2 = {'text_nomintervalle','text','interval name ',15,''};
form{length(form)+1} = {col1,colj,col2,coljg,colj};

col2 = {'nom_intervalle','edit',' ',15,'interval name'};
form{length(form)+1} = {col1,colj,col2,coljg,colj};

col2 = {'ajouter','radio','add',0,'Add a new interval'};
col3 = {'supprimer','radio','suppress',0,'Suppress an interval'};
form{length(form)+1} = {col1,colj,col2,col3,colj};

form{length(form)+1} = {col1,colj,coljg,coljg,colj};
form{length(form)+1} = {col1,colj,coljg,coljg,colj};

% Separation
sepa ={'separation_comm','frame','',10,''};
form{length(form)+1} = {sepa};

% Ligne de commentaires
col3 = {'commentaire','text','ceci est une ligne de commentaires pour les messages d''erreur ',10,'',[],'','visible','off'} ;
form{length(form)+1} = {col1,colj,col3};

% Separation
sepa ={'separation_comm','frame','',5,''};
form{length(form)+1} = {sepa};

%hout=zuicreeform('intervalle','intervalle','zuiedit_param_intervalle_fct','zuiedit_param_intervalle_ctrl',form) ;
hout=zuicreeform('interval','intervalle','zuiedit_param_intervalle_fct','zuiedit_param_intervalle_fct',form) ;

[hform,hui] = zuiformhandle('intervalle');

% on sauve les valeurs de xdeb et xfin
setappdata(hui.valeur_xdeb,'xdprec',tdeb_calcul) ;
setappdata(hui.valeur_xfin,'xfprec',tfin_calcul) ;

% Creation de la partie graphique
hui.axes_plot=zuiplotin(hui.list_plot);
title('');
xlabel('time');
ylabel('');
hui.cons = line(x,y,'marker','o','color',[1 0 0],'linestyle','none');

hui.trait = line([0 0],[0 0],'linestyle','-','linewidth',2);

int = evalin('base','param.intervalle') ;
name_int = fieldnames(int) ;
xd=[] ;
xf=[] ;
name=[] ;
for i=1:length(name_int)
	val = getfield(int,name_int{i}) ;
 	xd = [xd val(1)] ;
	xf = [xf val(2)] ;
	name = strvcat(name,name_int{i}) ;
end
setappdata(hui.list_intervalle,'libel',name) ;
setappdata(hui.list_intervalle,'xd',xd) ;
setappdata(hui.list_intervalle,'xf',xf) ;
set(hui.trait,'xdata',[xd(1) xf(1)],'ydata',[0 0], ...
	    'linewidth',3,'visible','on','linestyle','-') ;

set(hui.nom_intervalle,'string',name(1,:)) ;

% les references
hui.line_1 = line(x,y,'linestyle','-','visible','off');
hui.line_2 = line(x,y,'linestyle','-','visible','off');
hui.line_3 = line(x,y,'linestyle','-','visible','off');
hui.line_4 = line(x,y,'linestyle','-','visible','off');
hui.line_5 = line(x,y,'linestyle','-','visible','off');

zoom on
% memorisation des handles
setappdata(hform,'zhandle',hui);

% mise a jour de la liste des inetrvalles
zuilisteintrvl('intervalle');
