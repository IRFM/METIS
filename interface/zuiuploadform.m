% ZUIUPLOADFORM   charge dans un formulaires les donnees provenant du workspace
%------------------------------------------------------------------------------
% fichier zuiuploadform.m ->  zuiuploadform
%
% fonction Matlab 5 :
%	Cette fonction charge dans un formulaires les donnees provenant du workspace.
%	Les variables sont decrites dans l'appdata "variable" des uicontrol concernes. 
%
% syntaxe  :
%  zuiuploadform(hfig);
%
% entrees :
%  hfig : handle de la fenetre formulaire
%
% sortie : aucune
%
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 3.0, du 04/04/2005.
% 
% 
% liste des modifications : 
%
%  * 20/05/2003 -> protection des popup
%  * 20/05/2003 -> creation automatique des variables manquantes dans le workspace
%  * 04/04/2005 -> utilise la fonction same a la place == pour les cells
%
%--------------------------------------------------------------
%
function zuiuploadform(hfig)

% recherche des uicontrols
hui = findobj(hfig,'type','uicontrol');
% boucle sur les objets
for k =1:length(hui)
	hc  = hui(k);
	st  = get(hc,'style');
	% recupere le nom de la variable
	if isappdata(hc,'variable')
		var = getappdata(hc,'variable');
	else
		var = [];
	end
	% recupere le user data dans 
	user = get(hc,'userdata');
	% lit la donnees dans le workspace
	if ~isempty(var) & ~strcmp(var,'void');
	   try
	     donnees = evalin('base',var);
	     % selon l'objet charge la donnee
	     % cas de la table de translation
	     if ~isempty(user)
	     	% donnees format texte
	     	if ischar(donnees)
	     		value = strmatch(donnees,user);
	     		if isempty(value)
	     			stw = sprintf('Donnee non valide pour le popupmnu  (%s)',var);
	     			warning(stw)
	     		elseif length(value) >1
	     			value = min(strmatch(donnees,user,'exact'));
	     		end
	     		if ~isempty(value)
 			   % securite popupmenu
		 	   if strcmp(st,'popupmenu')
                              if isfinite(value)
  	     			 set(hc,'value',value);
                              end
                           else
 	     			set(hc,'value',value);
                           end
	     		end
	     	else
	     		% donnee numerique
	     		if iscell(user)
	     			for klk =1:length(user)
	     				if same(user{klk},donnees)
	     					value =klk;
	     				end
	     			end
	     		else
	     			value = min(find(user == donnees));
	     		end
	     		
	     		if isempty(value)
	     			stw = sprintf('Donnee non valide pour le popupmnu  (%s)',var);
	     			warning(stw)
	     		else
		 	   if strcmp(st,'popupmenu')
                              if isfinite(value)
  	     			 set(hc,'value',value);
                              end
                           else
	     			set(hc,'value',value);
                           end
	     		end
	     	end
	     elseif  strcmp(st,'edit')
	     	if ischar(donnees)
	     		set(hc,'string',donnees);
            elseif iscell(donnees)
                if ischar(donnees{:})
                    set(hc,'string', donnees{:});      
               else    
                    set(hc,'string',sprintf('%g',donnees{:}));      
                end
	     	else
	     		set(hc,'string',sprintf('%g',donnees));
	     	end
	     else
	     	 % tous les autres objets
	     	if ischar(donnees) & ~strcmp(st,'popupmenu')
	     		set(hc,'string',donnees);
	     	else
                     if strcmp(st,'text')
                           set(hc,'string',sprintf('%g',donnees));
		     elseif strcmp(st,'popupmenu')
                         if isfinite(value)
  	     	              set(hc,'value',value);
                          end
                     else
	     		set(hc,'value',donnees);
                     end
	     	end
	     end            	
	  catch
	    stw = sprintf('Variable introuvable dans le workspace (%s)',var);     
	    warning(stw);
            zassignin('base',var,[]);
	  end
	end	  
end
