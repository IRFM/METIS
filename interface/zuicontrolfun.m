% ZUICONTROLFUN fonction de control des formulaires crees avec zuicreefunform
%-------------------------------------------------------------------------------------
% fichier zuicontrolfun.m ->  zuicontrolfun
%
% fonction Matlab 5 :
% Cette fonction controle les valeurs des formulaires crees avec zuicreefunform . 
%
% syntaxe  :
%  zuicontrolfun(action)
%
% entrees :
%  action      = action associee au callback
%
% sortie : 
%
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.6, du 28/08/2001.
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function zuicontrolfun(action)

% cas special
hfig = get(0,'currentfigure');
if ~ishandle(hfig)
	return
end

h   = getappdata(hfig,'zhandle');
try
    hoc = getfield(h,action);
catch
    [action1,reste] = strtok(action,'.');
    [action2,reste] = strtok(reste,'.');
    hoc = h.(action1).(action2);
end
% si l'objet n'est pas un edit
if ~strcmp(get(hoc,'style'),'edit')
	% pas de control de valeur (aucun risque d'erreur)
	return
else
	data = zuidata(hoc);
	if ischar(data)
		% pas de control possible sur les chaines de caracteres
		return
	end
end
% nom de la fonction associee
fonction =get(hfig,'tag') ;

% supression des numeros dans action
% recupere le nom de la variable associee
var = getappdata(hoc,'variable');
% recupere le nombre
ind =findstr(var,'(');
if ~isempty(ind);
	nombre = var((ind+1):(end-1));
	actionc = strrep(action,nombre,'');
	nombre =int2str(nombre);
else
	actionc = action;
	nombre  = 1;
end


% recupere les infos
declaration = feval(fonction);
type        = zgetfield(declaration.type,actionc);
borne       = zgetfield(declaration.borne,actionc);
defaut      = zgetfield(declaration.defaut,actionc);

% cas des entiers
if ~isempty(findstr(lower(type),'integer'))|~isempty(findstr(lower(type),'entier'))
	data =fix(data);
end

% les bornes
if iscell(borne)
	flag = 0;
	for k =1:length(borne)
		if borne{k} == data
			flag =1;
		end
	end
	if flag == 0
		data =defaut;
    end
% pour les Inf ou NaN
elseif isempty(type)
	if isempty(data)
		data = defaut;
        elseif ischar(borne)
		return
	elseif data < min(borne)
		data = min(borne);
	elseif data > max(borne)
		data =max(borne);
	end

elseif strmatch(type,{'char','string','text'})
    return
else
	if isempty(data)
		data = defaut;
    elseif ischar(borne)
		return
	elseif data < min(borne)
		data = min(borne);
	elseif data > max(borne)
		data =max(borne);
	end
end

zuidata(hoc,data,defaut);

