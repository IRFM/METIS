% ZUIDOWNLOADFORM  mise a jour variables du workspace a partir des donnees du formulaire
%--------------------------------------------------------------------------------------
% fichier zuidownloadform.m 
%
% fonction Matlab 5 :
%	Cette fonction met a jour les variables du workspace a partir des donnees du formulaire.
%	Les variables sont decrites dans l'appdata "variable" des uicontrol concernes. 
%
% syntaxe  :
%	zuidownloadform(hfig);
%
% entrees :
%	hfig       = handle de la fenetre formulaire
%
% sortie : aucune
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 2.1, du 20/05/2003.
% 
% liste des modifications : 
%  * 07/09/2001 -> modifie le flag savenonok que si une donnee des structures
%                 data, param, post est modifiee.
%  * 20/05/2003 -> creation automatique des variables manquantes dans le workspace
%
%--------------------------------------------------------------
%
function zuidownloadform(hfig)

% recherche des uicontrols
hui = findobj(hfig,'type','uicontrol');
% variable de modififcation a sauver
modif = 0;
% boucle sur les objets
for k =1:length(hui)
    hc  = hui(k);
    % recupere le nom de la variable
    if isappdata(hc,'variable')
        var = getappdata(hc,'variable') ;
    else
        var = [];
    end
    if ~isempty(var)
        % est-ce data, param ou post
        modif = zmodifsave(var,modif);
        %style
        st  = get(hc,'style');
        % recupere le user data dans
        user = get(hc,'userdata');
        % valeur
        value = get(hc,'value');
        % string
        string = get(hc,'string');
        % selon l'objet
        % cas de la table de translation
        if ~isempty(user)
            % selon le type de donnees
            if ischar(user)
                donnees = deblank(user(value,:));
            elseif iscell(user)
                donnees = user{value} ;
            elseif length(user) == prod(size(user))
                donnees = user(value);
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
                donnees = string;
            else
                try
                    %string
                    donnees = eval(string);
                    if isstruct(donnees)
                        donnees = string;
                    end
                catch
                    donnees = string;
                end
            end
            if isempty(donnees)
                if isempty(string)
                    donnees =[];
                else
                    donnees =string;
                end
            end
        else
            % tous les autres objets
            donnees = value;
        end
        % ecrit la donnees dans le workspace
        try
            % verification des types
            model = evalin('base',var) ;
            if strcmp(st,'text')
                %rien
            elseif ischar(model) 
                if ischar(donnees)
                    zassignin('base',var,donnees);
                elseif isnumeric(donnees)
                    zassignin('base',var,num2str(donnees));
                elseif iscell(donnees)
                    zassignin('base',var,cat(1,donnees{:}));
                end
            elseif isnumeric(model)
                if ischar(donnees)
                    zassignin('base',var,num2str(donnees));
                elseif isnumeric(donnees)
                    if isnumeric(model) && (all(size(donnees)==1)) &&  ~isempty(model)
                        % Chan  debut modif
                        zassignin('base',var,donnees*ones(size(model)));
                    else
                        zassignin('base',var,donnees);
                    end
                elseif iscell(donnees)
                    zassignin('base',var,cat(1,donnees{:}));
                end
            elseif iscell(model)
                if ischar(donnees)
                    for ll=1:size(donnees,1)
                        v{ll}=donnees(ll,:);
                    end
                    zassignin('base',var,v);
                elseif isnumeric(donnees)
                    zassignin('base',var,donnees);
                elseif iscell(donnees)
                    zassignin('base',var,cat(1,donnees{:}));
                end
            else
                %	     			keyboard
                stw   = sprintf('incompatible type (%s)',var);
                warning(stw);
            end
        catch
            stw   = sprintf('unknown variable in the workspace (%s)',var);
            warning(stw);
            zassignin('base',var,[]);
        end
    end
end
if modif == 1
	zuisavenonok;
end

function modif = zmodifsave(var,modif)

tok =strtok(var,'.');
if strmatch(tok,{'data','param','post'})
	modif = 1;
end


