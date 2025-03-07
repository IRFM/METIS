% ZUIDEDIT_IMPORT_MODE_FCT   gestion des callbacks du formulaire d'importation de donnees
%--------------------------------------------------------------
% fichier zuiedit_import_mode_fct.m  
%
% fonction Matlab 5 :
%	fonction de gestion des callbacks des uicontrols du formulaire
%	d'importation de donnees
%
% syntaxe :
%	zuiedit_import_mode_fct(action)
%
% entrees :
%  action =  tag du uicontrol active
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 2.2, du 18/08/2004.
% 
% liste des modifications : 
% 18/08/2004 -> choix de l'occurence 0 rendu possible (1 a 9 avant seulement)
% 18/08/2004 -> securites sur nom de donnee, espace, et temps enlevees dans le cas ou on importe un profil de la base
%--------------------------------------------------------------
function zuiedit_import_mode_fct(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end

% initialisation des varibles
name ='';

% recupere le handle de la fenetre concernee
%[hfig,h] = zuiformhandle('import_mode') ;
hfig = gcf;
h    = getappdata(hfig,'zhandle');

% information pour l'assistant
zuicr(hfig,action) ;

mode   = getappdata(hfig,'mode') ;
source = zuidata(h.pop_srce) ;
type   = zuidata(h.pop_file) ;

% Fichier
if source==1
	set(h.text_nom,'String','nom du fichier') ;
	set(h.text_nom,'TooltipString','Nom du fichier') ;
	set(h.edit_nom,'TooltipString','Nom du fichier') ;
	zuienable(h.text_nom) ;
	zuienable(h.edit_nom) ;
	zuienable(h.text_type) ;
	zuienable(h.pop_file) ;
	zuienable(h.choix) ;

	% mode consigne 
	if mode==0
		% fichier Matlab
		if type==0
			set(h.edit_nom_tps,'style','edit') ;
			set(h.edit_nom_donnee,'style','edit') ;
		% fichier ASCII
		elseif type==1
			set(h.edit_nom_tps,'style','text') ;
			set(h.edit_nom_tps,'String','1ère colonne') ;
			set(h.edit_nom_donnee,'style','text') ;
			val = zuidata(h.edit_nom_espace) ;
			if strcmp(val,' ') val=2 ; end
			st = sprintf('%sè colonne',num2str(val)) ;
			set(h.edit_nom_donnee,'String',st) ;
		end
	% mode profil 
	elseif  mode==1
		% fichier Matlab
		if type==0
			set(h.edit_nom_tps,'style','edit') ;
			set(h.edit_nom_espace,'style','edit') ;
			set(h.edit_nom_donnee,'style','edit') ;
		% fichier ASCII
		elseif type==1
			set(h.edit_nom_tps,'style','text') ;
			set(h.edit_nom_tps,'String','1ère colonne') ;
			set(h.edit_nom_espace,'style','text') ;
			set(h.edit_nom_espace,'String','1ère ligne') ;
			set(h.edit_nom_donnee,'style','text') ;
			set(h.edit_nom_donnee,'String','(2:end,2:end)') ;
		end
	
	end

% Workspace
elseif source==0
	zuidisable(h.text_nom) ;
	zuidisable(h.edit_nom) ;
	zuidisable(h.text_type) ;
	zuidisable(h.pop_file) ;
	zuidisable(h.choix) ;
	set(h.edit_nom_tps,'style','edit') ;
	set(h.edit_nom_espace,'style','edit') ;
	set(h.edit_nom_donnee,'style','edit') ;

% Arcad
else source==-1
	set(h.text_nom,'String','objet@choc.occ') ;
	set(h.text_nom,'TooltipString', ...
        'nom_producteur@numchoc.occurence ou nom_signal@numchoc.occurence') ;
	zuienable(h.text_nom) ;

	set(h.edit_nom,'TooltipString', ...
        'nom_producteur@numchoc.occurence ou nom_signal@numchoc.occurence') ;
	zuienable(h.edit_nom) ;

	zuidisable(h.text_type) ;
	zuidisable(h.pop_file) ;
	zuidisable(h.choix) ;

	set(h.edit_nom_tps,'style','edit') ;
	set(h.edit_nom_espace,'style','edit') ;
	set(h.edit_nom_donnee,'style','edit') ;
end		

% selon ation
switch lower(action)

% choisir
case 'choix'
	% dialogue
	name = '*' ;
	[file,path]=uigetfile(name,'Nom du fichier à lire ?') ;
	drawnow
	if ~ischar(file) | file==0
		warning('chargement annulé')
	else
		filename = strcat(path,file) ;
		[pf,nf,ext,ver]=fileparts(filename) ;
		if strcmp(ext,'.gz')
			[pf,nf,ext,ver]=fileparts(nf) ;
		end
		if strcmp(ext,'.mat') 
			set(h.pop_file,'value',1) ;
		else
			set(h.pop_file,'value',2) ;
		end
		zuidata(h.edit_nom,filename) ;
	end

	zuireset(h.choix) ;

case {'annulation','close'}
	zuicloseone(hfig) ;	
	
case {'init','raz'}
	zuiformvisible(hfig) ;
	
	zuiformreset(hfig) ;
	zuiuploadform(hfig) ;
	
	zuireset(h.raz) ;
	
case {'validation'}
% test sur nom fichier
	% fichier
	if source==1
		name = zuidata(h.edit_nom) ;
		if strcmp(name,' ')
			herror = warndlg('le champ nom de fichier est vide ','Probleme','modal') ;
			zuireset(h.validation) ;
			return
		end			
	end
	% Arcad
	if source==-1
		% validation du format nom@numchoc.occurence
		name = zuidata(h.edit_nom) ;
		if strcmp(name,' ')
			herror = warndlg('le champ objet@choc.occ est vide ','Probleme','modal') ;
			zuireset(h.validation) ;
			return
		end			
		% on cherche '@'
		ind = findstr(name,'@') ;
		if isempty(ind)
			st = sprintf(' Attention au format  => @ \n  nom_producteur@numchoc.occurence ou nom_signal@numchoc.occurence') ;
			herror = warndlg(st,'Probleme','modal') ;
			zuireset(h.validation) ;
			return
		else
			% on extrait le nom et on recupere le numchoc et l'occurence
			[nom,rest] = strtok(name,'@') ;
			choc       = strtok(rest,'@') ;
			numchoc    = str2num(choc) ;

			% on teste si le numchoc.occurence est numerique
			if isempty(numchoc)
				st = sprintf(' Attention au format dans \n  nom_producteur@numchoc.occurence ou nom_signal@numchoc.occurence \n  numchoc.occurence  est un champ numérique') ;
				herror = warndlg(st,'Probleme','modal') ;
				zuireset(h.validation) ;
				return
			else
				% on cherche le '.' entre numchoc et occurence
				ind = findstr(choc,'.') ;
				if isempty(ind)
					st = sprintf(' Attention au format  => . \n  numchoc.occurence ') ;
					herror = warndlg(st,'Probleme','modal') ;
					zuireset(h.validation) ;
					return
				else
					% on recupere l'occurence
					[choc,occ] = strtok(choc,'.') ;
					occ        = strtok(occ,'.') ;
					if ~isempty(occ)
						occ        = str2num(strtok(occ,'.')) ;
						% l'occurence doit etre comprise entre 0 et 9
						if occ<0 | occ>9 
							st = sprintf(' Attention au format dans \n nom_producteur@numchoc.occurence ou nom_signal@numchoc.occurence \n  0 <= occurence => 9') ;
							herror = warndlg(st,'Probleme','modal') ;
							zuireset(h.validation) ;
							return
						end
					else
						st = sprintf(' Attention au format  => . \n  numchoc.occurence ') ;
						herror = warndlg(st,'Probleme','modal') ;
						zuireset(h.validation) ;
						return
					end
				end
			end
		end						
	end

	% variable temps
	val = zuidata(h.edit_nom_tps) ;
	if strcmp(val,' ') & source > -1
		herror = warndlg('le nom de la variable temps doit etre un champ caractères non vide ','Probleme','modal') ;
		zuireset(h.validation) ;
		return
	end

	% variable espace
	val = zuidata(h.edit_nom_espace) ;
	% mode consigne
	if mode==0
		if ischar(val)
			herror = warndlg('la voie doit etre un integer','Probleme','modal') ;
			zuireset(h.validation) ;
			return
		else
			val = floor(val) ;
			zuidata(h.edit_nom_espace,val) ;
			if val<= 0
				herror = warndlg('la voie doit etre >0 ','Probleme','modal') ;
				zuireset(h.validation) ;
				return
			end
		end
		type = zuidata(h.pop_file) ;
		if type==1
			st = sprintf('%sè colonne',num2str(val)) ;
			set(h.edit_nom_donnee,'String',st) ;
		end
		
	% mode profil
	elseif mode==1 & source > -1
		if strcmp(val,' ')
			herror = warndlg('le nom de la variable espace doit etre un champ caractères non vide ','Probleme','modal') ;
			zuireset(h.validation) ;
			return
		end
	end		

	% variable donnee
	val = zuidata(h.edit_nom_donnee) ;
	if strcmp(val,' ')  & source > -1
		herror = warndlg('le nom de la variable donnée doit etre un champ caractères non vide ','Probleme','modal') ;
		zuireset(h.validation) ;
		return
	end

	% valeur hors intervalle
	val = zuidata(h.edit_valeur) ;
	if ~isnan(val) & ~isinf(val) & ischar(val)
		herror = warndlg('Attention la valeur hors intervalle doit être numérique, Nan ou Inf ','Probleme','modal') ;
		zuireset(h.validation) ;
		return
	end

	nom      = getappdata(hfig,'nom_mode') ;
	type     = zuidata(h.pop_file) ;
	temps    = zuidata(h.edit_nom_tps) ;
	espace   = get(h.edit_nom_espace,'string') ;
	data     = zuidata(h.edit_nom_donnee) ;
	% coord    = zuidata(h.pop_espace)
	coord    = 2 ;
	defaut   = zuidata(h.edit_valeur) ; 
	positif  = zuidata(h.pop_defini) ;

	% appel de zimport
	cmd = strcat('cr = zimport(''',nom,''',',num2str(mode),',',num2str(source), ...
                     ',''',name,''',',num2str(type),',''',temps,''',''',espace,''',', ...
                     '''',data,''',',num2str(coord),',',num2str(defaut),',',num2str(positif),');');
	eval(cmd);	

	zuireset(h.validation) ;
	zuicache(hfig) ;	

otherwise
	
	   
end

