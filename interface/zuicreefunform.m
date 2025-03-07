% ZUICREEFUNFORM cree les formulaires associes aux fonctions auodeclaratives
%-------------------------------------------------------------------------------------
% fichier zuicreefunform.m ->  zuicreefunform
%
%
% fonction Matlab 5 :
% 
% Cette fonction cree les formulaires associes aux fonctionx auodeclaratives. 
%
% syntaxe  :
%  
%  hform=zuicreefunform(fonction,racine,nombre,{inerte,code_retour})
%
% entrees :
% 
%  fonction      = nom de la fonction
%  racine        = racine de la structure de parametre de la fonction a modifier
%  nombre        = nombre d'elements (coupleurs, injecteurs) a parametrer
%  inerte        = si a 1 mode affichage info seul
%  code_retour   = callback appele a la fin de 'validation'
%
% sortie : 
%
%  hform     = handle du formulaire 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 28/02/2005.
% 
% 
% liste des modifications : 
% 
% * 28/08/2001 -> ajout du mode inerte
% * 12/09/2001 -> corrcetion formatage des edits
% * 27/11/2002 -> ajout du code_retour
% * 12/12/2002 -> interface en anglais
% * 13/11/2003 -> modification de la politique vis a vis des variables inexistantes
% * 19/01/2005 -> ajout des labels specifique pour ITC
% * 28/02/2005 -> ajout du mode multicolonnes pour les modules avec grand nombre de parametre
%
%--------------------------------------------------------------
%
function hform=zuicreefunform(fonction,racine,nombre,inerte,code_retour,update)

% test des entrees
if nargin < 2
	warning('Il faut donner le nom de la fonction et la racine de la structure concernee')
	return
elseif isempty(fonction)
	warning('Il faut donner le nom de la fonction')
	return
elseif isempty(racine)
	warning('Il faut donner la racine de la structure concernee')
	return
end
if nargin  < 3
	nombre =1;
elseif isempty(nombre)
	nombre =1
end
if nargin < 4
	inerte = 0;
elseif isempty(inerte)
	inerte = 0;
end
if nargin < 5
	code_retour  = '';
end

% global switch for advanced 
if isappdata(0,'GUI_GLOBAL_BASIC_MODE')
   global_basic_mode = 1;
else
   global_basic_mode = 0;

end
% with backward compatibility
basic_mode = 0;
if isappdata(0,'GUI_GLOBAL_BASIC_MODE') &&  (getappdata(0,'GUI_GLOBAL_BASIC_MODE') == 1)
	advanced = 0;
        basic_mode = 1;
        if iscell(fonction)
	  fonction = fonction{1};
	end  
elseif isappdata(0,'GUI_GLOBAL_BASIC_MODE') &&  (getappdata(0,'GUI_GLOBAL_BASIC_MODE') == 0)
	advanced = 0;
        basic_mode = 0;
        if iscell(fonction)
	  fonction = fonction{1};
	end  
% flag du mode advanced
elseif iscell(fonction)
	advanced = 1;
	fonction = fonction{1};
else
	advanced = 0;
end

% decodage du nom de la fonction pour la gestion des sections
if any(fonction =='#')
    % il s'agit de n'afficher que une sections 
    [fonction,reste] = strtok(fonction,'#');
    section = strtok(reste,'#');
else 
    section = [];
end

% recuperation de la structure info
try
   if nombre > 1
   	declaration = feval(fonction,nombre);
   else
     	declaration = feval(fonction);
   end 	
catch
	warning(['wrong function name or path : ',fonction])
	return
end
if isempty(declaration)
	warning(['fonction non auto declarante : ',fonction])
	return
end

% this section have been put here due to multimode implementation (advanced/basic)
% check if varaible exist in workspace
liste = fieldnames(declaration.valeur);
for k = 1:length(liste)
	variable = '';
        try
		nom     = liste{k};
		variable = strcat(racine,'.',nom);
		data    = zgetfield(declaration.valeur,nom);
                % test d'existance dans le workspace de la variable
		if ~existbase(variable)
                        zassignin('base',variable,data);
                        fprintf('ZUICREEFUNFORM : Initialization of %s in Matlab workspace\n',variable);
                end
         catch
	    if isempty(variable)
	      vairable = lasterr;
	    end
	    fprintf('ZUICREEFUNFORM : error creating variable %s in workspace\n',variable); 
         end
end

% traitement des sections
if isfield(declaration,'section')  && (advanced == 0)
  if isempty(section)
      hf = findobj(0,'type','figure','tag',sprintf('section_form_%s',fonction));
      if ~isempty(hf)
	  figure(hf)
	  return
      end
       % creation de l'interface pour les sections
       sections_liste = {};
       noms = fieldnames(declaration.section);
       for k=1:length(noms)
	  if isfield(declaration.section,noms{k})
                if isempty(strmatch( declaration.section.(noms{k}),sections_liste,'exact')) && ...
                   ((isfield(declaration,'mode') && ~isfield(declaration.mode,noms{k})) || (basic_mode == 0))
		      sections_liste{end+1} = declaration.section.(noms{k});
                end
          end
       end
       if isempty(sections_liste)
	    return
       end
       for k=1:length(sections_liste)
            if nargin > 5
		commande_liste{k} = sprintf('zuicreefunform(''%s#%s'',''%s'',%d,%d,''%s'',%d);',fonction,sections_liste{k},racine,nombre,inerte,code_retour,update);
            else
 		commande_liste{k} = sprintf('zuicreefunform(''%s#%s'',''%s'',%d,%d,''%s'');',fonction,sections_liste{k},racine,nombre,inerte,code_retour);
          end
       end
       % ajout bouton interface complete
       sections_liste{end+1} = 'All in one interface';
       if nargin > 5
	  commande_liste{end+1} = sprintf('zuicreefunform(''%s#%s'',''%s'',%d,%d,''%s'',%d);',fonction,sections_liste{end},racine,nombre,inerte,code_retour,update);
       else
	  commande_liste{end+1} = sprintf('zuicreefunform(''%s#%s'',''%s'',%d,%d,''%s'');',fonction,sections_liste{end},racine,nombre,inerte,code_retour);
       end
       switch fonction
       case {'metis4imas','metis4itm','zerod'}
	      fonction_label = 'METIS parameters';
       otherwise
	      fonction_label = fonction;
       end
       hf = zuimenuliste(sprintf('%s (sections list)',fonction_label),sections_liste,commande_liste);
       set(hf,'tag',sprintf('section_form_%s',fonction));

	% memorisation des handles des formulaires
	if isappdata(0,'formulaire_handle')
		fh  =  getappdata(0,'formulaire_handle');
	else
		fh = [];
	end
	fh  =  zsetfield(fh,get(hf,'tag'),hf);
	setappdata(0,'formulaire_handle',fh);

       drawnow
       return
  elseif ~strcmp('All in one interface',section)
       % affichage que de la section 
      % decodage et identification 
      try
		if ~isfield(declaration,'mode')
			declaration.mode = struct([]);
		end
		if ~isfield(declaration,'label')
			declaration.label = struct([]);
		end
		noms = fieldnames(declaration.section);
                nb = 0;
		for k=1:length(noms)	
                        if isfield(declaration.section,noms{k})
			  if ~strcmp(declaration.section.(noms{k}),section) 
			      declaration.valeur      = rmfield(declaration.valeur,noms{k});
			      declaration.type        = rmfield(declaration.type,noms{k});
			      declaration.borne       = rmfield(declaration.borne,noms{k});
			      declaration.defaut      = rmfield(declaration.defaut,noms{k});
			      declaration.info  	    = rmfield(declaration.info,noms{k});
			      if isfield(declaration.mode,noms{k})
				  declaration.mode  	    = rmfield(declaration.mode,noms{k});
			      end
			      if isfield(declaration.label,noms{k})
				  declaration.label  	    = rmfield(declaration.label,noms{k});
			      end
   			      declaration.section  	    = rmfield(declaration.section,noms{k});
                          else
                               nb = nb +1;
			  end
                        else
                               nb = nb +1;
			end
			
		end
                if nb <= 10
		    declaration.multicol = 0;
                end

      catch
	      warning(['incorrect declaration inside : ',fonction])
	      return
      end
       
  end
end


% global advanced/basic mode
if (isfield(declaration,'mode') && ~isempty(declaration.mode) && (basic_mode == 1))
    
    try
        if ~isfield(declaration,'mode')
            declaration.mode = struct([]);
        end
        if ~isfield(declaration,'label')
            declaration.label = struct([]);
        end
        noms = fieldnames(declaration.mode);
        for k=1:length(noms)
            declaration.valeur      = rmfield(declaration.valeur,noms{k});
            declaration.type        = rmfield(declaration.type,noms{k});
            declaration.borne       = rmfield(declaration.borne,noms{k});
            declaration.defaut      = rmfield(declaration.defaut,noms{k});
            declaration.info        = rmfield(declaration.info,noms{k});
            if isfield(declaration.mode,noms{k})
                declaration.mode  	    = rmfield(declaration.mode,noms{k});
            end
            if isfield(declaration.label,noms{k})
                declaration.label  	    = rmfield(declaration.label,noms{k});
            end
            if isfield(declaration,'section') && isfield(declaration.section,noms{k})
                declaration.section  	    = rmfield(declaration.section,noms{k});
            end
        end
        if length(fieldnames(declaration.valeur)) <= 10
            declaration.multicol = 0;
        end
        
    catch
        warning(['incorrect declaration inside : ',fonction])
        return
    end
    
end
if (global_basic_mode == 1) && isfield(declaration,'mode')
                declaration = rmfield(declaration,'mode');
end
% decodage et identification 
try
        valeur      = declaration.valeur;
        type        = declaration.type;
        borne       = declaration.borne;
        defaut      = declaration.defaut;
        info        = declaration.info;
        interface   = declaration.interface;
        description = declaration.description;
        helpurl     = declaration.help;
        gui         = declaration.gui;
        controle    = declaration.controle;
	if isfield(declaration,'mode')
		mode = declaration.mode;
	else
		mode = struct([]);
	end
	if isfield(declaration,'label')
		label = declaration.label;
	else
		label = struct([]);
	end
	
	if (~isempty(mode) && (advanced == 0))
		commande_mem.fonction       = fonction;
		commande_mem.racine         = racine;
		commande_mem.nombre         = nombre;
		commande_mem.inerte         = inerte;
		commande_mem.code_retour    = code_retour;
		if nargin >5
			commande_mem.update = update;
		end
		noms = fieldnames(mode);
		for k=1:length(noms)	
			valeur      = rmfield(valeur,noms{k});
        		type        = rmfield(type,noms{k});
        		borne       = rmfield(borne,noms{k});
        		defaut      = rmfield(defaut,noms{k});
        		info  	    = rmfield(info,noms{k});			
		end
	else
		commande_mem = [];
	end
		
catch
        warning(['incorrect declaration inside : ',fonction])
        return
end


% gestion des formatages
noms = fieldnames(info);
for k=1:length(noms)
      info.(noms{k}) = sprintf(info.(noms{k}));
end




% flag de passage en mode multicolonnes pour les scalaire vrai
if isfield(declaration,'multicol')
      multicol = declaration.multicol;
elseif isempty(valeur)
      multicol = 0;	
elseif length(fieldnames(valeur)) >= 25
      multicol = 1;
else
      multicol = 0;
end

% si la fonction a un gui specifique le travail est presque terminer
if ~isempty(gui)
	% ouverture du gui
	try
        	feval(gui,'init');
	catch
	   warning(sprintf('le gui (%s) defini dans la fonction %s produit une erreur (%s)',gui,fonction,lasterr))
	end
	return
end

% si l'interface a deja ete appele
[hform,hui] = zuiformhandle(fonction);
if ishandle(hform)
        zuicloseone(hform)
end

% creation de l'interface
% les fonctions de callback
callback  = 'zuifaitfunaction';
if isempty(controle)
	controle = 'zuicontrolfun';
end
tag      = fonction;
switch fonction
case {'metis4imas','metis4itm','zerod'}
      fonction_label = 'METIS parameters';
otherwise
      fonction_label = sprintf('module interface %s',fonction);
end
if inerte == 0
        if ~isempty(section)
	    titre    = sprintf('%s (section #%s)',fonction_label,section);
        else
	    titre    = sprintf('%s',fonction_label);
        end
else
        if ~isempty(section)
	    titre    = sprintf('%s (section #%s, mode consultation)',fonction_label,section);
        else
	    titre    = sprintf('%s (mode consultation)',fonction_label);
        end
end

% gestion de l'aide
% removed option for METIS (only useful for CRONOS)
%if isempty(helpurl)
%	helpurl =strcat('file:',which(fonction));
%end

% le formulaire
% titre et aide simplifiee + le bouton d'aide
col1 = {'description','text@left',titre,[],description};
if isempty(helpurl)
      form{1} = {col1};
else
      col2 = {'url_aide','help@right','Help',[],'to open help function',helpurl};
      form{1} = {col1,col2};
end
% separation
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
% separation
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};
% separation
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};


% liste des variables
if ~isempty(valeur) && (isstruct(valeur) && ~isempty(fieldnames(valeur)))
	liste = fieldnames(valeur);
        poscol = 0;

	% 1ere partie : les variables ne dependant pas du nombre d'objets (injecteur, antenne ...)
	for k =1:length(liste)
		nom     = liste{k};
		data    = zgetfield(valeur,nom);
		limite  = zgetfield(borne,nom);
		def     = zgetfield(defaut,nom);
		aide    = zgetfield(info,nom);
		format  = zgetfield(type,nom);
		variable = strcat(racine,'.',nom);
		comment = sprintf('(%s)',format);
		if isfield(mode,nom)
			mode_baba = zgetfield(mode,nom);
		else
			mode_baba = 'advanced';
		end
		if isfield(label,nom)
			label_nom = zgetfield(label,nom);
		else
			label_nom = nom;
		end
		
		
                % test d'existance dans le workspace de la variable
		if ~existbase(variable)
                        zassignin('base',variable,data);
                        fprintf('ZUICREEFUNFORM : Initialization of %s in Matlab workspace\n',variable);
                end
		% si la variable est une cell
		dxl    = evalin('base',variable);
		if iscell(dxl)
			variable = sprintf('%s{:}',variable);
		end
		%
		if (length(data) <= 1)|(ischar(data) & any(size(data)<=1))
			% c'est une variable unique
			% le texte
			col1 = {strcat('text_',nom),'text',label_nom,[],aide};
			% le bouton/popup
			if iscell(limite)
				% c'est un popup
				value =[];
				chaine ='';
				%comment = [comment,' {'];
				% boucle sur les valeurs possibles
				for lk = 1:length(limite)
					val = limite{lk};
					if isempty(val)
					    val = '<vide>';
					end
					% cas de la chaines de caracteres
					if ischar(val)
						% string du popup
						chaine = strvcat(chaine,val);
						%comment = [comment,val,','];
						if strcmp(def,val)
							% valeur par defaut
							value = lk;
						end
					else
						% cas numerique
						% string du popup
						chaine = strvcat(chaine,sprintf('%g',val));
						%comment = [comment,sprintf('%g',val),','];
						if val == def
							value = lk;
						end
					end					
				end
                                % securite popup
                                if isempty(value)
                                    value = 1;
                                elseif value <= 0
                                    value = 1;
                                end
				%comment(end)='}';
				% declaration de l'objet
				col2={nom,'popup',chaine,value,aide,limite,variable};
			else
				% c'est un edit
				if ischar(def)
					chaine = def;
				else
					chaine = sprintf('%g',def);
				end	
				if length(chaine) < 24
					chaine = sprintf('%24.24s',chaine);
				end
				col2 ={nom,'edit',chaine,[],aide,[],variable};
				comment = [comment,' ',mat2str(limite)];
			end	
			% le type de la donnee				
			col3 = {strcat('type_',nom),'text',comment,[],''};

                        % selon le cas 
                        if multicol == 1
                           % posol
                           poscol = poscol +1;
                           if poscol > 3
                                 poscol = 1;
                           end
                           % selon la position
                           if poscol == 1
			      % ajout au formulaire
			      form{length(form)+1} = {col1,col2,col3};
                           else
                              colj     = {'jump','jump',' ',[],''};
                              listcol  = form{length(form)};
			      form{length(form)} = cat(2,listcol,{colj,col1,col2,col3});
                           end
                        else
			   % ajout au formulaire
			   form{length(form)+1} = {col1,col2,col3};
                        end
		end
	end
	
        % il faux remplir de case vide
        if multicol == 1
            if poscol == 1
                   coltx     = {'txjump','text',' ',[],''};
                   listcol  = form{length(form)};
		   form{length(form)} = cat(2,listcol,{colj,coltx,coltx,coltx},{colj,coltx,coltx,coltx});

            elseif poscol == 2
                   coltx     = {'txjump','text',' ',[],''};
                   listcol  = form{length(form)};
		   form{length(form)} = cat(2,listcol,{colj,coltx,coltx,coltx});
            end
        end

	% separation
	sepa ={'separation_comm','frame','',-5,''};
	form{length(form)+1} = {sepa};
	% separation
	sepa ={'separation_comm','frame','',3,''};
	form{length(form)+1} = {sepa};
	% separation
	sepa ={'separation_comm','frame','',-5,''};
	form{length(form)+1} = {sepa};
	
	% 2ieme partie : les variables dependant du nombre d'objets (injecteur, antenne ...)
	if nombre >1
           % cas particulier de itc
           spe = 0;
           if strcmp(fonction,'zitc');
               try
                  compo=evalin('base','param.compo');
                  spe = 1;
              end
           end
  	   % ligne de titre
	   col1=  {'jump1','jump','jump',[],''};
	   col3=  {'jump3','jump','jump',[],''};
	   cols={};
	   for l=1:nombre
                if spe == 0
		  aide   = sprintf('gaz, injector or antenna # %d',l);
		  col2={strcat('colonne',int2str(l)),'text',int2str(l),[],aide,[]};
		  cols{l}= col2;
                else
                     % cas itc
                     imp_alias = zitc_equiv(compo.z(l),compo.a(l));
                     if isempty(imp_alias)
                           imp_alias = sprintf('(%d,%d)',compo.z(l),compo.a(l));
                     else
                           imp_alias = sprintf('%s(%d,%d)',imp_alias,compo.z(l),compo.a(l));
                     end
		    aide   = sprintf('gaz # %d',l);
		    col2={strcat('colonne',int2str(l)),'text',imp_alias,[],aide,[]};
		    cols{l}= col2;
                end
	   end
	   % ajout au formulaire
	   form{length(form)+1} = {col1,cols{:},col3};
	end
	
	% boucle sur les champs
	for k =1:length(liste)
		nom     = liste{k};
		data    = zgetfield(valeur,nom);
		limite  = zgetfield(borne,nom);
		def     = zgetfield(defaut,nom);
		aide    = zgetfield(info,nom);
		format  = zgetfield(type,nom);
		variable = strcat(racine,'.',nom);
		if isfield(mode,nom)
			mode_baba = zgetfield(mode,nom);
		else
			mode_baba = 'advanced';
		end
		if isfield(label,nom)
			label_nom = zgetfield(label,nom);
		else
			label_nom = nom;
		end
		%
		if (length(data) > 1) & (~(ischar(data) & any(size(data)==1)))
			% c'est une variable multiple
			% le texte
			col1 = {strcat('text_',nom),'text',label_nom,[],aide};
			% boucle sur le nombre d'antennes/injecteurs
			cols ={};
			for l =1:length(data)
				% initialisation de  comment
		               comment = sprintf('(%s)',format);
			       cas_cell=0;
				% le bouton/popup
				if iscell(limite)
					% c'est un popup
					value =[];
					chaine ='';
				   % comment = [comment,' {'];
					% boucle sur les valeurs possibles
					for lk = 1:length(limite)
						val = limite{lk};
						% cas de la chaines de caracteres
						if ischar(val)
							% string du popup
							chaine = strvcat(chaine,val);
						       %comment = [comment,val,','];
						        cas_cell  = 1;
							if strcmp(def,val)
								% valeur par defaut
								value = lk;
							end
						else
							% cas numerique
							% string du popup
							chaine = strvcat(chaine,sprintf('%g',val));
						        %comment = [comment,sprintf('%g',val),','];
							if val == def
								value = lk;
							end
						end					
					end
					% declaration de l'objet
				        %comment(end)='}';
					if cas_cell == 1
					   col2 = {strcat(nom,int2str(l)),'popup',chaine,value, ...
					           sprintf('%s (# %d)',aide,l),limite, ...
					           strcat(variable,'{',int2str(l),'}')};
					else
					   col2 = {strcat(nom,int2str(l)),'popup',chaine,value, ...
					           sprintf('%s (# %d)',aide,l),limite, ...
					           strcat(variable,'(',int2str(l),')')};
				        end
				else
					% c'est un edit
					if ischar(def)
						chaine = def;
					else
						chaine = sprintf('%g',def);
					end	
					if length(chaine) < 12
						chaine = sprintf('%12.12s',chaine);
					end
					col2 = {strcat(nom,int2str(l)),'edit',chaine,[], ...
					sprintf('%s (# %d)',aide,l),[], ...
					strcat(variable,'(',int2str(l),')')};
				   comment = [comment,' ',mat2str(limite)];
				end
				cols{l} = col2;
			end
			
			% le type de la donnee
			col3 = {strcat('type_',nom),'text',comment,[],''};
			% ajout au formulaire
			form{length(form)+1} = {col1,cols{:},col3};
		end
	end
elseif ~isempty(section)
    return
else
texte ={'texte_vide','text@center','no parameter for this function !',[],''};
form{length(form)+1} = {texte};
	
end

% separation
sepa ={'separation_comm','frame','',-5,''};
form{length(form)+1} = {sepa};
% separation
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};


if nargin >5
	if isempty(commande_mem)
		% description des boutons communs
		% bouton annulation /precedent
		comm{1}= {'annulation', ...
			'radio', ...
			'Cancel', ...
			0, ...
			'canceled in progress action'};
		
		% bouton Raz
		comm{2}= {'raz', ...
			'radio', ...
			'reset', ...
			0, ...
			'get initial values'};
		
		% espace
		comm{3}= {'espace_cmd', ...
			'jump', ...
			'unpeuplusadroite', ...
			0, ...
			''};
			
		% bouton update	
		comm{4}= {'update', ...
			'radio', ...
			'Update', ...
			0, ...
			'update data in matlab workspace without closing the window'};
		
		
		% bouton validation/suivant
		comm{5}= {'validation', ...
			'radio@right', ...
			'Ok', ...
			0, ...
			'data formular validation'};
	else
		% description des boutons communs
		% bouton annulation /precedent
		comm{1}= {'annulation', ...
			'radio', ...
			'Cancel', ...
			0, ...
			'canceled in progress action'};
		
		% bouton Raz
		comm{2}= {'raz', ...
			'radio', ...
			'reset', ...
			0, ...
			'get initial values'};
		
		% espace
		comm{3}= {'espace_cmd', ...
			'jump', ...
			'unpeuplusadroite', ...
			0, ...
			''};
			
		% bouton advanced	
		comm{4}= {'advanced', ...
			'radio', ...
			'Advanced', ...
			0, ...
			'expert mode with all keys'};
		
			
		% espace
		comm{5}= {'espace_cmd2', ...
			'jump', ...
			'unpeuplusadroite', ...
			0, ...
			''};
			
		% bouton update	
		comm{6}= {'update', ...
			'radio', ...
			'Update', ...
			0, ...
			'update data in matlab workspace without closing the window'};
		
		
		% bouton validation/suivant
		comm{7}= {'validation', ...
			'radio@right', ...
			'Ok', ...
			0, ...
			'data formular validation'};
	
	end
else
	if isempty(commande_mem)
		% description des boutons communs
		% bouton annulation /precedent
		comm{1}= {'annulation', ...
			'radio', ...
			'Cancel', ...
			0, ...
			'canceled in progress action'};
		
		% bouton Raz
		comm{2}= {'raz', ...
			'radio', ...
			'reset', ...
			0, ...
			'get initial values'};
		
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
			'data formular validation'};
		
	else
		% description des boutons communs
		% bouton annulation /precedent
		comm{1}= {'annulation', ...
			'radio', ...
			'Cancel', ...
			0, ...
			'canceled in progress action'};
		
		% bouton Raz
		comm{2}= {'raz', ...
			'radio', ...
			'reset', ...
			0, ...
			'get initial values'};
		
		% espace
		comm{3}= {'espace_cmd', ...
			'jump', ...
			'unpeuplusadroite', ...
			0, ...
			''};
			
		% bouton advanced	
		comm{4}= {'advanced', ...
			'radio', ...
			'Advanced', ...
			0, ...
			'expert mode with all keys'};
		
			
		% espace
		comm{5}= {'espace_cmd2', ...
			'jump', ...
			'unpeuplusadroite', ...
			0, ...
			''};
			
		% bouton validation/suivant
		comm{6}= {'validation', ...
			'radio@right', ...
			'Ok', ...
			0, ...
			'data formular validation'};
		
	
	end
end

% creation du formulaire 
hform=zuicreeform(titre,tag,callback,controle,form,comm,0,0);
if ~(verLessThan('matlab', '8.2') & ~verLessThan('matlab', '7.9'))
  set(hform,'visible','off');
end
%zuiformreset(hform);
%zuiuploadform(hform);
feval(callback,'init',hform);

% mode inerte
if inerte == 1
	%drawnow
	[hfig,hui] = zuiformhandle(tag);
	set(hui.validation,'visible','off');
	set(hui.raz,'visible','off');
	set(hui.annulation,'string','Quit');
end

% ajout code retour
setappdata(hform,'code_retour',code_retour);

% pour le mode complet
if ~isempty(commande_mem)
	setappdata(hform,'commande_mem',commande_mem);
end
% fenetre visible
if ~(verLessThan('matlab', '8.2') & ~verLessThan('matlab', '7.9'))
  set(hform,'visible','on');
end
%drawnow

