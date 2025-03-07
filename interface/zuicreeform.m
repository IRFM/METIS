% ZUICREEFORM fonction de creation de GUI a base de formulaire
%--------------------------------------------------------------
% fichier zuicreeform.m ->  zuicreeform, gmarge
%
%
% fonction Matlab 5 :
% 
% Cette fonction cree une feuille (interface graphique) sous 
% la forme d'un formulaire. Le formulaire est constitue de lignes.
% Chaque ligne peut avoir un nombre de colonnes quelconques. Les lignes
% qui se suivent et qui ont un nombre de colonnes egal, voient leurs colonnes
% alignees verticalement.
%
% syntaxe  :
%  
%  hout=zuicreeform(titre,tag,fonction,control,form,{comm,nb,menubar,p1,v1, ...pn,vn});
%
% entrees :
% 
%  titre       =  titre de la fenetre (chaine caracteres)
%  tag         =  immatriculation de la fenetre (chaine de caracteres)
%  fonction    =  nom de la fonction de gestion des callback (elle est associe a 
%                 chaque uicontrol qui n'est pas lui meme associe a une variable)
%  control     =  nom de la fonction qui assure le control du contenu du formulaire 
%                 (elle est associe a chaquee uicontrol qui est lui meme associe
%                  a une variable)
%  form        =  vecteur de cellule qui decrit le formulaire
%  comm        =  vecteur de cellule qui decrit le formulaire commun [si non standart, optionnel].
%                 Si comm n'est pas un vecteur de cell, le formulaire n'a pas de partie commune.
%  nb          =  1 -> fenetre avec numero, 0 -> pas de numero [defaut 0, optionnel]
%  menubar     =  1 -> menu d'edition, 0 -> pas de menu d'edition [defaut 0, optionnel]
%  p1,v1 ... 
%  pn,vn       = autres proprietes passees a la commande figure (optionnel)
% 
%
% sorties :
% 
%  hout       = handle de la fenetre graphique
%
% description du formulaire :
% 
%  la varibale form est un vecteur de "cell". Chaque element du vecteur decrit une ligne.
%  Chaque ligne est decrite par un vecteur de "cell". Chaque element du vecteur (ligne)
%  decrit un case du formulaire.
%  Chaque case est drecite par un vecteur de "cell" constituer comme suit :
%  
%   {tag,style,string,value,tooltip,userdata,var,p1,v1, ... pn,vn}
%   
%   tag      = immatriculation du uicontrol, utilise pour le callback
%   style    = style de l'objet (cf. plus loin)
%   string   = chaine de caractere de l'objet
%   value    = valeur de l'objet
%   tooltip  = aide contextuelle
%   userdata = donnees associees a l'objet (utiliser pour les popup)
%   var      = variable asociee a l'objet dans le wokspace (pour les mise a jour automatique)
%   p1,v1, ...
%   pn,vn    = autres proprietes passe au uicontrol (optionnel)
% 
%   remarque 1 : longueur des zones d'edition 
%   Les objet edit sont initialiser avec la chaine "string" et ont une longueur
%   resever pour l'edition "value" en nombre de caracteres
%   
%   remraque 2 : la gestion des table de translation (par exemple pour les popupmenu)
%   Si le userdata de l'objet n'est pas vide, il sert de table de translation. Le user data
%   contient un vecteur (numerique ou de cell) ou une matrice de text.
%   Chaque element du vecteur correspond a un choix dans l'objet (dans le meme ordre).
%   La table de translation est appliquee sur le champ "value".
%   
%   remarque 3 : la variable style
%   tous les styles admis par les uicontrol sont psossible, plus le style "jump"
%   qui provoque un saut de colonne. La largeur du saut est donne par le nombre 
%   de caracteres de la chaine "string". Il est possible d'ajouter un bouton de 
%   style "help" par formulaire qui ouvre une url ou donne l'aide de la fonction. 
%   L'argument doit etre passe dans le champ userdata. 
%   
%   remarque 4 : les frame
%   Si une frame a une donnee "value" non nulle, alors
%   elle est traitee comme un ligne de separation de "value" pixels de haut. elle est
%   visible si la valeur est >0, sinon elle est invisible
%   
%   remarque 5 : modification du style
%   La variable style peut etre complete par une commande. le format est style@commande.
%   Les commandes possibles sont : 
%         right  -> impose que l'objet soit aligner a droite de la fenetre 
%         center -> impose que l'objet soit au centre de la fenetre (horizontalement)
%         left   -> impose que l'objet soit aligner a gauche de la fenetre 
%         full   -> impose que l'objet occupe toute la fenetre
%         merge  -> impose que l'objet soit etendu a la colonne suivante (de droite), 
%                   l'operateur est transitif a gauche.
%   
%   remarque 6 : objet "listbox" "edit" et  "frame"
%   Les objets "listbox" occupe plusieurs lignes si les cases du dessous (du meme groupe) contiennent des 
%   objets invisibles (jump)
%   
% les fonctions de callback :
% 
%  Elles recoivent comme argument de retour le tag de l'objet qui a declenche le callback.
%  La fonction "fonction", si elle recoit l'argument 'init' elle doit creer le fenetre 
%  si necesssaire, mettre a jour les donnees du formulaire et la rendre visible.
%  Si elle recoit l'argument 'close' elle doit rendre la fenetre inivisible
%  Si elle recoit l'argument 'annulation' elle doit rendre la fenetre invisible
%  Si elle recoit l'argument 'validation' elle doit modifier les donnees dans l'espace de 
%  travail, sauver le fichier et rendre le fenetre invisible.
%  Si elle recoit l'argument 'raz' elle rinitialiser le formulaire.
%  
% variables d'environements :
% 
% Elles sont stockees dans les appdata de l'objet root. Si elles ne sont
% pas definies, la fonction utilise des valeurs par defauts.
% 
%  La variable presentation : 
%  c'est une structure. Elle contient les informtions de formatage : 
% 
%    presentation.dx               = espace horizontal entre deux boutons
%    presentation.dy               = espace verticale entre deux boutons
%    presentation.ym               = marge verticale autour du texte
%    presentation.xm.default       = marge horizontal par defaut autour du texte
%    presentation.xm.pushbutton    = marge selon style
%    presentation.xm.togglebutton  = marge selon style
%    presentation.xm.radiobutton   = marge selon style
%    presentation.xm.checkbox      = marge selon style
%    presentation.xm.edit          = marge selon style
%    presentation.xm.text          = marge selon style
%    presentation.xm.slider        = marge selon style
%    presentation.xm.frame         = marge selon style
%    presentation.xm.listbox       = marge selon style
%    presentation.xm.popupmenu     = marge selon style
%
%  La variable bouton_commun : 
%  elle contient un formulaire (meme syntaxe que form) qui est ajouter en bas
%  de chaque fenetre. Il est commun a toutes les fenetres.
%
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.2, du 05/04/2004.
% 
% 
% liste des modifications : 
%
%  * 09/07/2001 -> ajout de l'operateur 'merge'
%  * 20/07/2001 -> correction de la taille des fenetres
%  * 20/07/2001 -> option multilangue
%  * 27/07/2001 -> rajout de defaut.style    : on garde le style �la cr�tion des uicontrol (C. Passseron)
%                            dafaut.etat     :  "  "    l'�at  �la cr�tion     "
%  * 11/10/2001 -> ajout des references croisees pour les fenetres
%  * 02/04/2004 -> ajout de la protection des fenetres de l'interface
%  * 05/04/2004 -> remise du patch qui a disparu pour le popup et qui empeche le plantage de matlab sous linux
%  * 12/04/2004 -> extention du multiligne a "edit"
%  * 18/07/2005 -> tentative de diminution du debit reseau pour l'affichage matlab7
%
%-----------------------------------------------------------------------------------------------------------
%
function hout=zuicreeform(titre,tag,fonction,control,form,comm,nb,menubar,varargin)


%dbstack

% info pour les referenes croisees
info_0   = get(0);
info_gcf = get(info_0.CurrentFigure);

% creation de la fenetre
% choix des parametres
if nargin < 6
   comm = {};
end
if nargin < 7
   nb = 0;
end
if nargin <8
   menubar = 0;
end
if nargin <9
   varargin = [];
end
if isempty(titre)
   nb = 1;
elseif isempty(nb)
   nb = 0;
end
if isempty(tag)
   tag = titre;
end
if isempty(menubar)
   menubar = 0;
end
if isempty(fonction)
   error('il faut une fonction de callback')
end
if isempty(control)
   control ='';
end
if nb  == 0
   nb  = 'off';
else
   nb  = 'on';
end
if menubar == 0
   menubar = 'none';
else
   menubar = 'figure';
end

% creation de la fenetre
hout = figure('name',titre, ...
	      'numbertitle',nb, ...
	      'menubar',menubar, ...
	      'tag',tag, ...
	      'units','pixels', ...
	      'visible','off', ...
	      'resize','off', ...
	      'defaultuicontrolHorizontalAlignment','center', ...
              'HandleVisibility','callback',...
              'CloseRequestFcn',strcat(fonction,'(''close'');'));

% memorisation des handles des formulaires
if isappdata(0,'formulaire_handle')
	fh  =  getappdata(0,'formulaire_handle');
else
	fh = [];
end
fh  =  zsetfield(fh,tag,hout);
setappdata(0,'formulaire_handle',fh);


% autres donnees attache a la fenetre
setappdata(hout,'control',control);
setappdata(hout,'fonction',fonction);
setappdata(hout,'form',form);

% gestion de la fermeture
%set(hout,'CloseRequestFcn',strcat(fonction,'(''close'');'));

%parametres de control de la presentation
if isappdata(0,'presentation')
	par = getappdata(0,'presentation');
else
	par = [];
end

if isempty(par)
   par.dx               = 1;    % espace horizontal entre deux boutons
   par.dy               = 1;    % espace verticale entre deux boutons
   %par.ym               = 1;    % marge verticale autour du texte
   if ismac
	  par.ym               = 10;    % marge verticale autour du texte
   else
	  par.ym               = 1;    % marge verticale autour du texte    
   end
   par.xm.default       = 1;    % marge horizontal par defaut autour du texte
   par.xm.pushbutton    = 3;    % selon style
   par.xm.togglebutton  = 1;    % selon style
   par.xm.radiobutton   = 30;   % selon style
   par.xm.checkbox      = 1;    % selon style
   par.xm.edit          = 3;    % selon style
   par.xm.text          = 1;    % selon style
   par.xm.slider        = 1;    % selon style
   par.xm.frame         = 1;    % selon style
   par.xm.listbox       = 30;    % selon style
   par.xm.popupmenu     = 30;    % selon style
end
   
% ajout des bouttons communs a toutes les fenetres
if isempty(comm)
   if isappdata(0,'bouton_commun')
   	comm = getappdata(0,'bouton_commun');
   end
end
if isempty(comm)
   % description des boutons communs
   % bouton annulation /precedent
   comm{1}= {'annulation', ...
	     'radio', ...
	     'Cancel', ...
	     0, ...
	     ''};

   % bouton Raz
   comm{2}= {'raz', ...
	     'radio', ...
	     'reset', ...
	     0, ...
	     ''};

   % espace
   comm{3}= {'espace_cmd', ...
	     'jump', ...
	     'unpeuplusadroite', ...
	     0, ...
	     ''};

   % bouton validation/suivant
   comm{4}= {'validation', ...
	     'radio@right', ...
	     'Ok', ...
	     0, ...
	     ''};
end
if iscell(comm)
   sepa ={'separation_comm','frame','',-5,''};
   form{length(form)+1} = {sepa};
   form{length(form)+1} = comm;
end

% creation des objets
harray = [];    % tableau des handles

% calcul de la taille de la fenetre
taille_x = 0;
taille_y = 0;

% boucle sur les lignes
nbligne = length(form);
for k = 1:nbligne
   % ligne courante
   ligne = form{k};
   nbcol = length(ligne);
   % boucle sur les colonnes
   tx =0;
   ty =0;
   for l = 1:nbcol
      % description de l'objet
      col = ligne{l};
      % proprietes des objets 
      tag     = col{1};
      style   = col{2};
      % decodage style
      [style,reste] = strtok(style,'@');
      modifieur     = strtok(reste,'@');
		
      string  = col{3};
      value   = col{4};
      if isempty(col{5})
            tooltip = '';
      else
            tooltip = col{5};
      end
      if length(col) >= 6
         userdata = col{6};
      else
         userdata = [];
      end
      if length(col) >= 7
         var = col{7};
      else
         var = '';
      end
      if length(col) > 7
         prop = col(8:end);	 
      else
         prop = {};
      end
      if strcmp(style,'jump')
         style ='text';
         visible = 'off';
      else
         visible ='on';
      end
      flag_help =0;  
      if strcmp(style,'help')
         style ='push';
         flag_help =1;
      end

      % protection anti crash pour linux
      if strcmp(style,'popup')
	  if isempty(value)
               value =1;
          elseif ~isfinite(value)
               value =1;
          else
               value = max(1,value);
          end
      end
	
      % uicontrol
      harray(k,l) = uicontrol(hout, ...
                   'units','pixels', ...
		   'interruptible','off', ...
		   'busyaction','queue', ...
		   'style',style, ...
		   'tag',tag, ...
		   'value',value, ...
		   'string',string, ...
		   'tooltip',tooltip, ...
		   'visible',visible, ...
		   'userdata',userdata, ...
		    prop{:});
		                       
      if flag_help == 1    
         if ~isempty(findstr(userdata,'.m'))
            set(harray(k,l),'callback',strcat('zne(''',userdata,''');'));
         else
            set(harray(k,l),'callback',strcat('netscape(''',userdata,''');'));
         end
      elseif isempty(var)
         set(harray(k,l),'callback',strcat(fonction,'(''',tag,''');'));
      else
         if ~isempty(control)
            set(harray(k,l),'callback',strcat(control,'(''',tag,''');'));
         end
         setappdata(harray(k,l),'variable',var);
      end
		 
      % memorisation des valeurs par defaut
      defaut.string = string ;
      defaut.value  = value ;
      defaut.style  = style ;
      defaut.etat   = get(harray(k,l),'Enable') ;

      setappdata(harray(k,l),'init_defaut',defaut);

      % memorisation du modifieur de style
      setappdata(harray(k,l),'modifieur',modifieur);

      % calcul des tailles                                  
      ttc         = get(harray(k,l),'extent');
      if ~isempty(value)
         if (strcmp(style,'edit') |strcmp(style,'text'))& (value > 1)
            ttc(3) = max(ttc(3),ceil(ttc(3) .* value ./ max(length(string),1)));
         end
      end
      %tx          = tx + ttc(3) + par.dx + gmarge(par,style);
      %ty          = max(ty,ttc(4));
      %taille{k,l} = ttc+[0,0,gmarge(par,style),par.ym];
      aux_gmarge=gmarge(par,style);
      gmarge_aux  = gmarge(par,style);
      tx          = tx + ttc(3) + par.dx + gmarge_aux;
      ty          = max(ty,ttc(4));
      taille{k,l} = ttc+[0,0,gmarge_aux,par.ym];
   end
   taille_x = max(taille_x,tx);
   taille_y = taille_y + ty + par.dy + par.ym; % -1;
   longueur(l) = tx;
end

% taille de la fenetre a peu pres
taille_x = taille_x - par.dx;
taille_y = taille_y - par.dy;
setappdata(hout,'taille_x',taille_x);
setappdata(hout,'taille_y',taille_y);

% sauvegarde des handles
setappdata(hout,'harray',harray);
% tableau des dimensions
setappdata(hout,'taille',taille);
% vecteur des longueurs de ligne
setappdata(hout,'longueur',longueur);

% calcul des position
% 1 - nombre de colonnes :
for l =1:size(harray,1)
%   ind      = find(harray(l,:)~=0);
   nbcol(l) = length(find(harray(l,:)~=0));
end
setappdata(hout,'nbcol',nbcol);

% 2 - calcul des positions
% initialisation de la boucle  
l =1;
% position initialle
ypos  = taille_y;
% tableau des hauteurs
hauteur = zeros(1,size(harray,1));

% boucle sur les lignes		
while(l <=size(harray,1))
    
   % nb de colonnes de la section
   nbc = nbcol(l);
    
   % tableau taille de chaque colonne
   tcol = zeros(1,nbc);
    
   % recherche changement d'ordonnancement suivant
   ind = find(nbcol ~= nbc);
   ind = min(ind(find(ind>l)));
   if isempty(ind) 
      lf =size(harray,1);
   else
      lf = ind -1;
   end	
   index =l:lf;
    
   % taille des colonnes
   % boucle sur les colonnes
   for k=1:nbc
      % allignement des colonnes
      for m =l:lf
         ttc =taille{m,k};
         tcol(k) = max(tcol(k),ttc(3));
         hauteur(m) = max(hauteur(m),ttc(4));
      end
   end
   % mise a jour de la taille x
   tx = sum(tcol);
   taille_x =max(taille_x,tx);
    
   % creation des positions
   % boucle sur les lignes
   for m =l:lf
      % boucle sur les colonnes
      for k=1:nbc
         % calcul de la position	
         ttc =taille{m,k};
         if k == 1
            x0 =0;
         else
            x0 =sum(tcol(1:(k-1)))+ (k-1) .* par.dx;
         end 
         dx =tcol(k);
         y0 = ypos - hauteur(m);
         style = get(harray(m,k),'style');
         if strcmp(style,'frame')
            val = get(harray(m,k),'value');
            if val >0
               dy =val;
            elseif val <0
               dy =abs(val);
               set(harray(m,k),'visible','off');
            else
               dy = hauteur(m);
            end
         else
            dy = hauteur(m);
         end
         y0 = ypos - dy;
	    
         % placement
         pos = [x0,y0,dx,dy];
         set(harray(m,k),'position',pos);
      end
      ypos = ypos - dy - par.dy;
   end
    
   % traitement des objets qui occupent plusieurs lignes
   % boucle sur les colonnes
   for k=1:nbc
      style=get(harray(l,k),'style');
      stringl=get(harray(l,k),'string');
      
      % cas du listobx
      if (strcmp(style,'listbox') | strcmp(style,'frame')| ...
              (strcmp(style,'edit') & (length(stringl) > 16))) & (lf > l)
         dl =l;
         % recherche des objets invisbles
         multi = 1;
         for m =(l+1):lf
            if strcmp(get(harray(m,k),'visible'),'off') && (multi ==1)
               dl =m;
            else
               multi = 0;
            end
         end
         % si plusieur lignes modification des position
         if dl >l
            posini = get(harray(l,k),'position');
            posfin = get(harray(dl,k),'position');
            pos(1) = posini(1);
            pos(3) = posini(3);
            pos(4) = posini(2) + posini(4) - posfin(2);
            pos(2) = posfin(2);
            set(harray(l,k),'position',pos);
         end	
      end
   end
    
   % nouvelle serie de ligne
   l = lf + 1;
% fin du while
end    

% redimensionnement final de la fenetre
taille_x = taille_x - par.dx;
taille_y = taille_y - par.dy;
posf = get(hout,'position');
if ispc
    
    set(hout,'position',[100,100,taille_x,taille_y]);
else
    set(hout,'position',[posf(1:2),taille_x,taille_y]);
end
setappdata(hout,'taille_x',taille_x);
setappdata(hout,'taille_y',taille_y);

% formatage des lignes de separation
hfr = findobj(hout,'type','uicontrol','style','frame','visible','on');
for k=1:length(hfr)
	pos =get(hfr(k),'position');
	pos(3) = taille_x -pos(1)+1;
	set(hfr(k),'position',pos);
end

% application des modification de position
posf = get(hout,'position');
for k = 1:size(harray,1)
	for l=1:size(harray,2)
		hc = harray(k,l);
		if ishandle(hc)
			if isappdata(hc,'modifieur')
				modifieur = getappdata(hc,'modifieur');
			else
				modifieur = [];
			end
			if ~isempty(modifieur)
				switch modifieur,
				case 'right'
					pos =get(hc,'position');
					pos(1) = posf(3)-pos(3)+1;
					set(hc,'position',pos);
				case 'left'
					pos =get(hc,'position');
					pos(1) = 0;
					set(hc,'position',pos);
				case 'center'
					pos =get(hc,'position');
					pos(1) = fix((posf(3) - pos(3))/2);
					set(hc,'position',pos);
				case 'full'
					pos =get(hc,'position');
					pos(1) = 0;
					pos(3) =posf(3);
					set(hc,'position',pos);	
			   end
			end
		end
	end
end
% allongement de l'objet avec merge
for k = 1:size(harray,1)
	for l=(size(harray,2)-1):-1:1
		hc = harray(k,l);
		hcs = harray(k,l+1);
		if ishandle(hc)&ishandle(hcs)
			if isappdata(hc,'modifieur')
				modifieur = getappdata(hc,'modifieur');
			else
				modifieur = [];
			end
			if ~isempty(modifieur)
				switch modifieur,
				case 'merge'
					pos1 =get(hc,'position');
					pos2 =get(hcs,'position');
               				txt1 = get(hc,'string');
               				txt2 = get(hcs,'string');
               				pos  = [pos1(1:2),pos2(1) - pos1(1) + pos2(3),pos1(4)];
					set(hc,'position',pos,'string',[txt1,txt2]);
					set(hcs,'visible','off');
			   	end
			end
		end
	end
end


% appdata des handle sous forme de structure tag
hfr = findobj(hout,'type','uicontrol');
zhandle =[];
for k=1:length(hfr)
	zhandle=zsetfield(zhandle,get(hfr(k),'tag'),hfr(k));
end
setappdata(hout,'zhandle',zhandle);

% correction du decalage vertical
% calcul de l'offset
posend = get(harray(end,1),'pos');
yoff   = posend(2) - 1 - par.dy;
% boucle de changement des offset
hh= harray(:);
for k =hh'
	if k ~=0
		pos = get(k,'position');
		pos(2) = pos(2) - yoff;
		set(k,'position',pos);
	end
end
% reduction de la fenetre
pos1 =get(harray(1,1),'position');
ymax = pos1(2) + pos1(4) + par.dy;
posf =get(hout,'position');
posf(4) = ymax;
set(hout,'position',posf);


% passe en unites normalisees
%drawnow  
%zuiplace(hout);
%drawnow
child  = findobj(hout,'type','uicontrol');
set(child,'unit','normalized');
set(hout,'unit','normalized');
%set(child,'unit','pixel');
%set(hout,'unit','pixel');

if isempty(varargin)
	set(hout,'resize','off');
	set(hout,'visible','on');
	drawnow
else
   	set(hout,varargin{:});
end

%drawnow
% appel de la fonction des references croisees
%try
%  zinterface_crossref(info_0,info_gcf,hout);
%catch
%  disp('Erreur dans zinterface_crossref !');
%end
% fin de la fonction

  
    
% fonction qui recherche les parametres du style    
function marge=gmarge(par,style)

   data  = par.xm; 
   field = fieldnames(data);   
   ind   = strmatch(style,field);
   marge = zgetfield(data,field{ind});
   
   
