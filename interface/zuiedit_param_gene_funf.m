% ZUIEDIT_PARAM_GENE_FUNF creation formulaires de convergence, sauvegarde, execution des parametres generaux
%--------------------------------------------------------------
%
% fichier zuiedid_param_gene_funf.m ->  
%		zuicreeform : creation du formulaire
%
% fonction Matlab 5 :
%	fonction de definition des uicontrols pour les formulaires
%	de convergence, sauvegarde, execution des parametres generaux
% 	sous le mode edition du formulaire principal
%
% syntaxe  : 
%	zuiedid_param_gene_funf(fonction,titre,callback,controle) ;
%
% entrees
%		fonction : nom de la fonction concernee : 
%			zconv pour paramatres de convergence
%			zsauv pour les parametres relatifs a la sauvegarde des resultats
%			zequa pour les parametres  de configuration des equations
%			zexec pour les parametres commandant l'execution
%			zprofile pour la configuration du profiler
%		titre	: titre du formulaire
%		callback : nom de la fonction des gestion des callbacks
%                si vide on utilise 'zuifaitfunaction'
%		control : nom de la fonction de controls des donnees
%		           si vide on utilise 'zuicontrolfun'
%
% sorties :
% 
% fonction ecrite par C. Passeron, poste 61 19
% version 2.0, du 12/11/2002.
%
% 
% liste des modifications : 
%  * 07/09/2001 -> ajout de l'argument racine
%  * 07/09/2001 -> modification de la longueur de la chaine de caracteres
%                   pour le cas d'un edit
%  * 12/11/2002 -> correction bug matlab6
%
%--------------------------------------------------------------

function zuiedit_param_gene_funf(fonction,racine,titre,callback,controle)

if nargin < 1
	warning('You must give the function name ')
	return
elseif isempty(fonction)
	warning('You must give the function name ')
	return
end

if nargin < 2
	racine = 'param.gene';
end

if nargin < 3
	titre='' ;
end

if nargin < 4
	callback  = 'zuifaitfunaction';
elseif isempty(callback)
	callback  = 'zuifaitfunaction';
end

if nargin < 5
	controle = 'zuicontrolfun';
elseif isempty(controle)
	controle  = 'zuicontrolfun';
end


% si l'interface a deja ete appele
[hform,hui] = zuiformhandle(fonction) ;
if ishandle(hform)
        zuiformvisible(hform)
	return
end

% recuperation de la structure info
declaration = feval(fonction);

if isempty(declaration)
	warning(['unfind function  : ',fonction])
	return
end

% decodage et identification 
try
        valeur      = declaration.valeur;
        type        = declaration.type;
        borne       = declaration.borne;
        defaut      = declaration.defaut;
        info        = declaration.info;
catch
        warning(['wrong function declaration : ',fonction])
        return
end

% creation de l'interface
tag      = fonction;

form = {};
sepa ={'separation_comm','frame','',3,''};
form{1} = {sepa};

col1 = {'libel_loadfile','text@full','Global simulation parameters',[],''};
form{length(form)+1} = {col1};

col1 = {'libel_loadfile','text@full',titre,[],''};
form{length(form)+1} = {col1};

% le formulaire
% separation
sepa ={'separation_comm','frame','',3,''};
form{length(form)+1} = {sepa};

% liste des variables
if ~isempty(valeur)
	liste = fieldnames(valeur);

	colj=  {'jump1','jump','jump',[],''};
	% tester les champs

	% 1ere partie : les variables ne dependant pas du nombre d'objets (injecteur, antenne ...)
	for k =1:length(liste)
		nom     = liste{k};

		champ = strcat('declaration.valeur','.',nom);
		eval(strcat('teststruct=isstruct(',champ,');'))
		if teststruct
			eval(strcat('champnew=fieldnames(',champ,');')) ;  	
 			for l=1:length(champnew)
				nomnew=strcat(nom,'.',champnew{l}) ;
				data    = zgetfield(valeur,nomnew) ;
				limite  = zgetfield(borne,nomnew) ;
				def     = zgetfield(defaut,nomnew) ;
				aide    = zgetfield(info,nomnew) ;
				format  = zgetfield(type,nomnew) ;
				variable = strcat(racine,'.',nomnew) ;
				comment='';
				if ~isempty(format)
 				     comment = sprintf('(%s)',format) ;
				end
		%
				if (length(data) == 1)|(ischar(data) & any(size(data)==1))
					% c'est une variable unique
					% le texte
					col1 = {strcat('text_',nomnew),'text',nomnew,[],aide};
					% le bouton/popup
					if iscell(limite)
						% c'est un popup
						value =[];
						chaine ='';
						%comment = [comment,' {'];
						% boucle sur les valeurs possibles
						for lk = 1:length(limite)
							val = limite{lk};
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
						%comment(end)='}';
						% declaration de l'objet
						col2={nomnew,'popup',chaine,value,aide,limite,variable};
					else
						% c'est un edit
						if ischar(def)
							chaine = def ;
						else
	 						chaine = num2str(def) ;
						end	
						col2 ={nomnew,'edit',chaine,1,aide,[],variable} ;
						comment = [comment,' ',mat2str(limite)] ;
					end	
					% le type de la donnee				
					col3 = {strcat('type_',nomnew),'text',comment,[],''} ;
					% ajout au formulaire
					form{length(form)+1} = {colj,col1,col2,col3};
				end
			end
		else
			data    = zgetfield(valeur,nom);
			limite  = zgetfield(borne,nom);
			def     = zgetfield(defaut,nom);
			aide    = zgetfield(info,nom);
			format  = zgetfield(type,nom);
			variable = strcat(racine,'.',nom);
			comment='';
			if ~isempty(format)
			comment = sprintf('(%s)',format);
			end

			if (length(data) == 1)|(ischar(data) & any(size(data)==1))
				% c'est une variable unique
				% le texte
				col1 = {strcat('text_',nom),'text',nom,[],aide};
				% le bouton/popup
				if iscell(limite)
					% c'est un popup
					value =[];
					chaine ='';
					%comment = [comment,' {'];
					% boucle sur les valeurs possibles
					for lk = 1:length(limite)
						val = limite{lk};
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
					%comment(end)='}';
					% declaration de l'objet
					col2={nom,'popup',chaine,value,aide,limite,variable};
				else
					% c'est un edit
					if ischar(def)
						chaine = def ;
					else
						chaine = num2str(def) ;
					end	
% modification du 07/09/2001
					col2 ={nom,'edit',chaine,length(chaine)+2,aide,[],variable};
					comment = [comment,' ',mat2str(limite)];
				end	
				% le type de la donnee				
				col3 = {strcat('type_',nom),'text',comment,[],''};
				% ajout au formulaire
				form{length(form)+1} = {colj,col1,col2,col3};
			end
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
	
end

% creation du formulaire 
hform=zuicreeform('Edition',tag,callback,controle,form);
zuiformreset(hform);
%zuiuploadform(hform);
