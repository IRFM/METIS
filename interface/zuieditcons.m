% ZUIEDITCONS cree la fenetre de l'editeur de consignes ou de profil
%---------------------------------------------------------------------
% fichier zuieditcons.m ->  zuieditcons
%
% fonction Matlab 5 :
%
%	Cette fonction cree la fenetre de l'editeur de consignes ou de profil.
%	L'editeur de consigne permet de dessiner la consigne ou le profil a la
%	souris et au clavier  en placant les points des noeuds. Il est possible
%	d'enlever ou d'ajouter des points. Lors de la validation les
%	variables designe par var_x et var_y sont mises a jour et la fonction
% 	code_retour' est executee.
% 
%	Les noeuds sont des 'o' rouge'.
%	La consignes est trait bleu.
%	La courbe visualiser est celle qui est retournee dans le workspace. 
% 
%	Pour changer la valeur d'un noeud : 
%		1 - cliquer sur le noeud ou selectionner un point dans la liste
%		2 - deplacer le point a la souris ou changer les valeurs dans les champs X et Y
%		   (l'abscisse ndoit rester dans l'intervalle entre les noeuds voisins)
%  
%	Pour supprimer des noeuds :
%		1 - selectionner un ou plusieurs noeuds dans la liste ou cliquer  sur un noeud 
%		2 - appuyer sur le bouton supprimer
%  
%	ou avec le bouton droit de la souris, cliquer sur un noeud et selectionner "supprimer"
%	dans le menu contextuel
%   
%	Les noeuds des des extremites ne peuvent pas etre supprimer
% 
%	Pour ajoutezuieditcons.mr des noeuds :
%		1 - selectionner un ou plusieurs noeuds dans la liste ou cliquer  sur un noeud 
%		2 - appuyer sur le bouton ajouter
%		    si plusieurs noeuds sont selectionnes, un noeud est cree entre chaque noeud selectionne 
%		    (a milieu de chaque intervalle)
%		    si un seul noeud est selectionne, un noeud est cree a droite de ce noeud.
%  
%	ou avec le bouton droit de la souris, cliquer sur un noeud et selectionner "ajouter"
%	dans le menu contextuel,  un noeud est cree a droite de ce noeud.
%    
%	ou avec le bouton droit de la souris, cliquer sur dans les axes et selectionner "ajouter"
%	dans le menu contextuel,  un noeud est cree a la position designee.
%    
%	Il est impossible de cree des noeud en dehors de l'intervalle de temps initial.
%    
%	La representation : 
%		Le bouton spline permet de passer d'une interpolation lineaire a une interpolation spline.
%		La consigne ou le profil retourne est celui tracer en bleu a l'ecran.
% 
%	L'echelle : 
%		Le multiplicateur est appliqu a la consigne ou au profil lors de la validation.
% 
%	Les references : 
%		Les boutons de couleur en haut de l'ecran permettent d'afficher jusqu'a trois courbes de references.
%      
%	Les quadrillage :
%		Le bouton en haut a gauche permet d'affciher un quadrillage.  
%     
%	La modulation des profils : 
%		La valeur d'un profil peut etre modulee (temporellement) par une consigne, si le bouton portant le nom de la consigne est 
%		selectione. Le mdoule de la consigne normalise a 1 est utilise.
%
%
%	syntaxe  :
%	hout=zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal,code_retour,liste_ref,var_ref,texte_prop,var_prop,varargin);
%    
% entree :
%	nom           = nom de la fenetre (mot simple sans pontuation)
%	aide          = tooltip associe a la fenetre
%	x             = abscisse de la consigne [nbt,1] ou du profil [1,nbrho]
%					    (les varaibles de retour sont echantillonees sur cette base)
%	y             = valeur initiale de la consigne [nbt,1] ou du profil [1,nbrho]
%	texte_x       = xlabel des axes de la fenetre d'edition
%	texte_y       = ylabel des axes de la fenetre d'edition
%	var_x         = nom de la variable de retour dans le workspace pour les abscisse (si inutile, mettre a vide '')
%	var_y         = nom de la variable de retour dans le workspace pour les donnees de la consigne ou du profil
%				       (si inutile, mettre a vide '')
%	canal         = numero de la consigne editee pour les groupe (par exemple pour le deuxieme coupleur fci, canal = 2 
%                  et la variable mise a jour est data.cons.fci(:,2) )
%	code_retour   = prend les valeurs :
%		* '' (vide) -> recopie les variables et execute "zuisavenonok" 
%		* 'abs'     -> edition du module de la consigne, recompose la consigne avec sa phase initiale,
%		               puis recopie les variables et execute "zuisavenonok" 
%		* 'angle'   -> edition de la phase de la consigne, recompose la consigne avec son amplitude initiale,
%		               puis recopie les variables et execute "zuisavenonok" 
%		* 'degres'  -> comme angle mais la phase est editee en degres
%		               pour toutes les autres valeur execute la fonction design� par code_retour a la fin du callback 
%		               du bouton "validation", si vide execute "zuisavenonok"  (cette fonction permet des operations 
%		               supplemenetaires sur la consigne).
%	liste_ref     = liste de noms [vecteur de "cell" contenant des chaines de caracteres] des signaux pouvant 
%	                etre affiche avec la consigne comme des references. le dernier elements doit etre la chaine
%	                'vide'.
%	var_ref       = description des variables du workspace prises pour reference : c'est un  vecteur de  "cell". 
%	                Chaque "cell" est de la forme : 
%							{'nom de la variable abscisse dans le workspace', ...
%                     'nom de la variable donnee dans le workspace', ...
%                     'style de la ligne pour le plot'}
%                  Le style 'o' trace de 'o' aux points non nul (utilise pour la consigne des glacons)
%                  le dernier element du vecteur doit etre :{'[]','[]',''}
%                     
%	texte_prop    = texte affiche dans le bouton commandant la modulation temporelle du profil
%	var_prop      = nom de la variable du workspace servant de consigne de modulation temporelle du profil  
%                  (Si pas de modulation doit contenir le nom de la variable de temps pour le profil. 
%                  Doitetre vide pour les consignes)
%       varargin   Options passees a zuicreeform
%                                     
%
% sorties : 
%	hout  =  handle de la fenetre cree.
% 
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.2, du 15/09/2003.
% 
% 
% liste des modifications : 
%	* 28/08/2001 -> ajout du mode basetemps
%	* 30/08/2001 -> correction bug initialisation x et y
%	* 13/09/2001 -> correction bug mode 'degres'
%	* 20/09/2001 -> correction bug references par defaut
%	* 25/09/2001 -> tag fenetre dependant du canal
%	* 08/10/2001 -> correction bug abs(y) dans le cas par defaut
%  * 29/03/2002 -> ajout des codes de retour pour l'IdN
%  * 17/10/2002 -> ajout des codes de retour pour te et ti
%  * 15/09/22003 -> ajout du mode nhnd et ntnd
%  * 15/09/22003 -> ajout du mode toro et  polo pour fce
%
%--------------------------------------------------------------
%
function hout=zuieditcons(nom,aide,x,y,texte_x,texte_y,var_x,var_y,canal,code_retour,liste_ref,var_ref,texte_prop,var_prop,varargin)

% cas specifique de Te,Pe, Ti & Pion
switch var_y
case 'data.prof.te'
	evalin('base','param.edit.tepe =''te'';');
case 'data.prof.pe'
	evalin('base','param.edit.tepe =''pe'';');
case 'data.prof.ti'
	evalin('base','param.edit.tipion =''ti'';');
case 'data.prof.pion'
	evalin('base','param.edit.tipion =''pion'';');
end	

% tag du formulaire
% -----------------
if strmatch(code_retour,{'pfce','toro','polo','nbar','gaspuff','idn','ind2','ntnd','nhe3nd'})
   tagc = strcat(strrep(var_x,'.',''),strrep(var_y,'.',''),code_retour) ;
else
   tagc = strcat(strrep(var_x,'.',''),strrep(var_y,'.','')) ;
end
if ~isempty(canal)
   tagc = sprintf('%s_%d',tagc,canal);
end
   
% si l'interface a deja ete appelee
[hout,hui] = zuiformhandle([tagc '_editeur']) ;
if ishandle(hout)
       zuiformvisible(hout) ;
       return
end

% mode de test si pas d'argument
if nargin <10
	code_retour = '';
end
if nargin <11
    liste_ref = '';
     var_ref   = {};
end
if isempty(liste_ref) | isempty(var_ref)
	liste_ref = '     Ip    |Vloop|Flux|ICRH|ECRH|LH|NBI|Gas|Pumping|Zeffm|Pellet|Empty';
	var_ref   = {{'data.gene.temps','data.cons.ip',':'}, ...
	             {'data.gene.temps','data.cons.vloop',':'}, ...
	             {'data.gene.temps','data.cons.flux',':'}, ...
	             {'data.gene.temps','abs(data.cons.fci)',':'}, ...
	             {'data.gene.temps','abs(data.cons.fce)',':'}, ...
	             {'data.gene.temps','abs(data.cons.hyb)',':'}, ...
	             {'data.gene.temps','real(data.cons.idn)',':'}, ...
	             {'data.gene.temps','data.cons.c',':'}, ...
	             {'data.gene.temps','data.cons.pomp',':'}, ...
	             {'data.gene.temps','data.cons.zeffm',':'}, ...
	             {'data.gene.temps','data.cons.glacon','o'}, ...
	             {'[]','[]',''}};
end
if nargin < 13
	texte_prop = '';
	var_prop   = '';
end

% si length(x) < 2
modessc =0;
if length(x) <2
	x=[0,1];
	y=[0,0];
elseif length(x) ==2
	modessc = 1;
end

% orientation
if size(x,1) > 1
	spline_flag =0;
else
	spline_flag =1;
end
x=x(:);
y=y(:);

% selon le code retour
ymod   = abs(y);
yangle = angle(y);
yreal  = real(y);
yimag  = imag(y);
[ypuiss,ytoro,ypolo] = zdecodefce(y);

switch code_retour
case 'angle'
	y     = yangle;
	multi = 1;
	nom_titre = sprintf('angle(%s) {radian}',nom);
case 'degres'
	multi = 1;
	y     = yangle / pi * 180;
	nom_titre = sprintf('angle(%s) {degrees}',nom);
case 'polo'
	multi = 1;
	y     = ypolo;
	nom_titre = sprintf('poloidal angle (%s) {degrees}',nom);
case 'toro'
	multi = 1;
	y     = ytoro;
	nom_titre = sprintf('toroidal angle (%s) {degrees}',nom);
case {'abs','pfce'}
	y  = ymod;
	yy = ymod;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
	nom_titre = sprintf('abs(%s)',nom);

case 'idn'
	y  = yreal;
	yy = yreal;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
	nom_titre = sprintf('NBI, real(%s)',nom);

case 'idn2'
	y  = yimag;
	yy = yimag;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
	nom_titre = sprintf('NBI2, imag(%s)',nom);

case 'didndt'
	y  = yimag;
	yy = yimag;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
	nom_titre = sprintf('dPnbidt, imag(%s)',nom);

case 'nbar'
	y  = yreal;
	yy = yreal;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
	nom_titre = sprintf('nbar, real(%s)',nom);

case 'gaspuff'
	y  = yimag;
	yy = yimag;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
	nom_titre = sprintf('gaspuff, imag(%s)',nom);

case 'fisrt'
	y  = yreal;
	yy = yreal;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
	nom_titre = sprintf('%s - first',nom);
case 'second'
	y  = yimag;
	yy = yimag;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
	nom_titre = sprintf('%s - second',nom);
case 'nhnd'
	y  = yreal;
	yy = yreal;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
	nom_titre = sprintf('nH/nD, real(%s)',nom);
    
case 'nhe3nd'
	y  = yreal;
	yy = yreal;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
	nom_titre = sprintf('nHe3/nD, real(%s)',nom);
    
case 'ntnd'
	y  = yimag;
	yy = yimag;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
	nom_titre = sprintf('nT/nD, imag(%s)',nom);
otherwise
	%y  = y;
	yy = ymod;
	% choix du multiplicateur
	ind   = find(yy > eps);
	multi = 10.^fix(log(median(abs(yy(ind))))./log(10));
	if isempty(multi)
		multi = 1;
	elseif ~isfinite(multi)
		multi = 1;
	end
	nom_titre =nom;
end

if modessc == 1
	multi =1 ;
end


% extraction des points
x_mem = x;
y_mem = y;
curve_init =0;
if all(~isfinite(y))
    x = linspace(min(x),max(x),11);
    y = zeros(size(x));
    multi = 1;
elseif length(x) <31
    % on garde tout les points
    % mode lineaire par defaut
    x = x';
    y = y';
else

   % choix des noeux
   [x,y]= zextraitnoeud(x,y,11 + 20 .* (1 - (spline_flag~=0)),spline_flag);

   curve_init = 1;
end

% formulaire de test
% avec la sous structure from
form={};
colj = {'void','jump','void',[],''} ;

% 1ere ligne nom de la machine
col1 = {'nom_var','text@full',['Edition of : ',nom_titre],[],aide};
form{1} = {col1};

% 2ieme ligne
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% 3ieme ligne
col1 = {'grid','popup','  grid X  | grid Y | grid XY | None ',4,'adds grid lines to the figure'};
col2 = {'jumpref','jump','12345',[],''};
if isempty(var_ref)
	col3 = {'ref_1','popup',liste_ref,length(var_ref),'plot the reference 1','','','visible','off'};
	col4 = {'ref_2','popup',liste_ref,length(var_ref),'plot the reference 2','','','visible','off'};
	col5 = {'ref_3','popup',liste_ref,length(var_ref),'plot the reference 3','','','visible','off'};
	col6 = {'jumpref2','jump','12345',[],''};
	col7 = {'node_1','popup',liste_ref,length(var_ref),'adds nodes from the reference','','','visible','off'};
else
	col3 = {'ref_1','popup',liste_ref,length(var_ref),'plot the reference 1',var_ref,'','BackgroundColor',[0 0.5 0]};
	col4 = {'ref_2','popup',liste_ref,length(var_ref),'plot the reference 2',var_ref,'','BackgroundColor',[0.75 0 0.75]};
	col5 = {'ref_3','popup',liste_ref,length(var_ref),'plot the reference 3',var_ref,'','BackgroundColor',[0 0.75 0.75]};
	col6 = {'jumpref2','jump','12345',[],''};
	col7 = {'node_1','popup',liste_ref,length(var_ref),'adds nodes from reference',var_ref,'','BackgroundColor',[0.75 0.5 0.25]};
end
% gestion de l'aide
helpurl =strcat('file:',which('zuieditcons'));

% le bouton d'aide
col8 = {'jumphelp','jump','1234581011',[],''};
col9 = {'url_aide','help','Help',[],'Open the help window',helpurl};
%form{length(form)+1} = {colj,col1,col2,col3,col4,col5,col6,col7,col8,col9};
form{length(form)+1} = {colj,col1,col2,col3,col4,col5,col6,col7,col8};

%4ieme ligne
col1 = {'list_plot','list','Ceci est un texte tres long pour reserver de la place & Ceci est un texte tres long pour reserver de la place',1,''};
col2 = {'axes_jump','jump','123',[],''};
col3 = {'liste_valeurs','list','#115 [1.45687,3.141617]',1,'list of reference points'};
form{length(form)+1} = {colj,col1,col2,col3};

% 15 lignes identiques pour faire de la place
col1 = {'jump_void','jump','Ceci est un texte tres long pour reserver de la place',[],''};
col2 = {'jump_void','jump','123',[],''};
col3 = {'jump_void','jump','#115 [1.45687,3.141617]',[],''};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};
form{length(form)+1} = {colj,col1,col2,col3};

% nieme ligne
sepa ={'separation_comm','frame','',10,''};
form{length(form)+1} = {sepa};

% ligne de titre de la ligne de commande
col1 = {'ajouter','radio','Add',0,'add a point to the   curve'};
col2 = {'supprimer','radio','Delete',0,'Delete the previous point'};
col3 = {'jump_cmd','jump','123',[],''};
col4 = {'text_X','text','X = ',[],''};
col5 = {'valeur_X','edit',sprintf('%12g',x(1)),[],'X coordinate of the selected point'};
col6 = {'text_Y','text','  Y = ',[],''};
col7 = {'valeur_Y','edit',sprintf('%12g',y(1)./multi),[],'Y coordinate of the selected point'};
col8 = {'jump_cmd','jump','12',[],''};
col10 = {'jump_cmd','jump','12',[],''};
col11 = {'text_multi','text','multiplier = ',[],''};
col12 = {'valeur_multi','edit',sprintf('%16g',multi),5,'curve multiplier'};
col13 = {'jump_cmd','jump','12',[],''};
col14 = {'spline','popup','linear|spline|pchip',spline_flag+1,'spline fit of the curve',{0,1,2},''} ;
if ~isempty(texte_prop)
    col15 = {'modulation','popup',strcat('constant|',texte_prop),1,'time modulation of the profile  (normalisation = 1)'};
else
    col15 = {'modulation','text','constant',1,'time modulation of the profile  (normalisation = 1)'};
end
form{length(form)+1} = {colj,col1,col2,col3,col4,col5,col6,col7,col8,col10,col11,col12,col13,col14,col15};


%tagc = strcat(strrep(var_x,'.',''),strrep(var_y,'.',''));

hout=zuicreeform('Waveform curve editor',[tagc,'_editeur'],'zuieditcons_action','zuieditcons_action',form ...
                ,{},0,0,varargin{:});

setappdata(hout,'code_retour',code_retour);
setappdata(hout,'multi_mem',multi);
setappdata(hout,'var_x',var_x);
setappdata(hout,'var_y',var_y);
setappdata(hout,'canal',canal);
setappdata(hout,'x_mem',x_mem);
setappdata(hout,'y_mem',y_mem);
setappdata(hout,'ymod',ymod);
setappdata(hout,'yangle',yangle);
setappdata(hout,'var_prop',var_prop);

% autres objets
[hform,hui] = zuiformhandle([tagc,'_editeur']);
hui.axes_plot=zuiplotin(hui.list_plot);
if modessc == 1
	set(hui.axes_plot,'ylim',[1e-7,10],'yscale','log');
end

title('');
xlabel(texte_x);
ylabel(texte_y);
hui.cons = line(x,y./multi,'marker','o','color',[1 0 0],'linestyle','none','hittest','on');
xx=linspace(min(x_mem),max(x_mem),length(x_mem));
if spline_flag == 0
    hui.interp = line(x,y./multi,'color',[0 0 1],'hittest','off');
elseif spline_flag == 1
    hui.interp = line(xx,zspline(cat(2,x(1)-mean(diff(x)),x),cat(2,y(1),y)./multi,xx),'color',[0 0 1],'hittest','off');
else
    hui.interp = line(xx,zpchip(cat(2,x(1)-mean(diff(x)),x),cat(2,y(1),y)./multi,xx),'color',[0 0 1],'hittest','off');
end
hui.marque = line(x(1),y(1)./multi,'color',[0 0 0],'hittest','off','linestyle','none','marker','*','userdata',1);

if curve_init == 1
    hui.curve = line(x_mem,y_mem./multi,'color',[0.25 0.25 0.25],'hittest','off','linestyle',':');
end  

% les references
hui.line_1 = line(x,y./multi,'color',[0 0.5 0],'linestyle',':','visible','off');
hui.line_2 = line(x,y./multi,'color',[0.75 0 0.75],'linestyle',':','visible','off');
hui.line_3 = line(x,y./multi,'color',[0 0.75 0.75],'linestyle',':','visible','off');


% le bouton raz n'a pas de sens
zuicloseone(hui.raz)
hui.raz=[];

% memorisation des handles
setappdata(hform,'zhandle',hui);

% cration des menu contextuels
hui.context = uicontextmenu;
set(hui.cons,'uicontextmenu',hui.context);
hui.item1 = uimenu(hui.context, 'Label', 'add', 'Callback', 'zuieditcons_action(''ajouter_souris'')');
hui.item2 = uimenu(hui.context, 'Label', 'delete', 'Callback','zuieditcons_action(''supprimer_souris'')');
hui.context_axes = uicontextmenu;
set(hui.axes_plot,'uicontextmenu',hui.context_axes);
hui.item3 = uimenu(hui.context_axes, 'Label', 'add', 'Callback', 'zuieditcons_action(''ajouter_axes'')');

% creation des callbacks de selection
%set(hui.cons,'ButtonDownFcn','zuieditcons_action(''souris'')');
% mise a jour de la liste
zuilistepoint([tagc,'_editeur']);

