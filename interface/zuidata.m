% ZUIDATA manipule les donnees des uicontrol
%-----------------------------------------------------------------
% fichier zuidata.m ->  zuidata
%
%
% fonction Matlab 5 :
% 
% Cette fonction manipule les donnees des uicontrol. Elle
% permet de mettre ou de recuperer une donnee d'un uicontrol.
%
% syntaxe  :
%  
%  pour lire :
%      donnees=zuidata(hc);
%      
% pour ecrire :
%      zuidata(hc,donnees,repli);
%
% entrees :
% 
%  hc     = handle du uicontrol
%  donnee = donnee a mettre dans le uicontrol 
%           (elle doit etre compatible avec le uicontrol)
%  repli  = si la donnee est incompatible, valeur  ou chaine 
%           mis a la place (selon type, mis dans value ou string 
%           sans utilisation de la table de translation)
%
% sortie : 
% 
%  donnee = donnee contenue dans le uicontrol 
%  
%  remarque : 
%    
%    si repli = 'zuimax' -> repli prend la valeur du champ max de l'objet
%    si repli = 'zuimin' -> repli prend la valeur du champ min de l'objet
%    
%  
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.1, du 20/05/2003.
% 
% 
% liste des modifications : 
%
%   * 20/05/2003 -> protection des popup
%
%--------------------------------------------------------------
%
function donnees=zuidata(hc,donnees,repli)

% test des entrees
if nargin  == 0 
	warning('zuidata: pas de handle')
	return
elseif ~ishandle(hc)
	warning('zuidata: handle invalide')
	return
end

if nargin <3
	repli = [];
elseif strcmp(repli,'zuimax')
	repli = get(hc,'max');
elseif strcmp(repli,'zuimin')
	repli = get(hc,'min');
end

%style
st  = get(hc,'style');
% recupere le user data dans 
user = get(hc,'userdata');
% valeur
value = get(hc,'value');
% string
string = get(hc,'string');

% selon les entrees
if nargin  < 2 
	% mode get
	% selon l'objet
	% cas de la table de translation
	if ~isempty(user)
		% selon le type de donnees
		if ischar(user)
			donnees = deblank(user(value,:));
		elseif iscell(user)
			donnees = user{value};
		elseif length(user) == prod(size(user))
			donnees = user(value) ;
		elseif ndims(user) == 2
			donnees =user(value,:);
		elseif ndims(user) == 3
			donnees =user(value,:,:);
		else
			donnees = user(value);
		end	
	elseif  strcmp(st,'edit')
                if (exist(string) >= 2) && (exist(string) <= 6) && ...
                   isempty(strmatch(lower(string),{'ans','eps','realmax','realmin','pi','i','inf','nan','j'},'exact'))                   

			% c'est le nom d'une fonction
			donnees =string;
		else
			try
				donnees = eval(string);
	%chan
				if isstruct(donnees) 
					donnees = string;
				end
	%chan
			catch
				donnees = string;
			end
                end
	   	if isempty(donnees)
		    donnees =string;
	   	end
   	else
	   	% tous les autres objets
	   	donnees = value;
	end        
else
	% cas set
	% selon l'objet charge la donnee 
	% cas de la table de translation
	if ~isempty(user)
		% donnees format texte
		if ischar(donnees)
			value = strmatch(donnees,user);
			if isempty(value)
				value = repli;
			elseif length(value) >1
				value = min(strmatch(donnees,user,'exact'));
			end
			% securite popupmenu
		 	if strcmp(st,'popupmenu')
                           if isempty(repli)
                              repli = 1;
                           elseif ~isfinite(repli)
                              repli = 1;
                           end 
                           if ~isfinite(value)
                                 value = repli;
                           end
			end
			if ~isempty(value)
				set(hc,'value',value);
			else
				set(hc,'value',repli);
			end
		else
			% donnee numerique
			value = min(find(user == donnees));
			if isempty(value)
				set(hc,'value',repli);
			else
				set(hc,'value',value);
			end
		end
	elseif  strcmp(st,'edit')
		if ischar(donnees)
			set(hc,'string',donnees);
		else
			set(hc,'string',sprintf('%g',donnees));
		end
	else
		% tous les autres objets
		if ischar(donnees) & ~strcmp(st,'popupmenu')
			set(hc,'string',donnees);
		else
			% securite popupmenu
		 	if strcmp(st,'popupmenu')
                           if isempty(repli)
                              repli = 1;
                           elseif ~isfinite(repli)
                              repli = 1;
                           end 
                           if ~isfinite(donnees)
                                 donnees = repli;
                           end
			end
			set(hc,'value',donnees);
		end
	end            	
end
