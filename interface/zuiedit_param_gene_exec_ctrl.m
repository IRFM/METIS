% ZUIEDIT_PARAM_GENE_EXEC_CTRL   control parametres relatifs a l'execution des parametres généraux
%--------------------------------------------------------------
% fichier zuiedit_param_gene_exec_ctrl.m ->
%
% fonction Matlab 5 :
%	fonction de controle des parametres relatifs a l'execution
%	des parametres  généraux
%	sous le mode edition du formulaire principal
%
% syntaxe :
%	zuiedit_param_gene_exec_ctrl(action)
%
% entrees
%	action       =  tag du uicontrol active
%
% sorties	
%	
% fonction ecrite par C. Passeron, poste 61 19
% version 1.3, du 10/04/2001.
% 
% liste des modifications : 
%
%--------------------------------------------------------------
function zuiedit_param_gene_exec_ctrl(action)

%  
[hfig,h] = zuiformhandle('zexec') ;
if ~ishandle(hfig)
	return
end
hoc = getfield(h,action) ;

tt = evalin('base','data.gene.temps') ;

% Partie zuicontrolfun 
% si l'objet n'est pas un edit
if ~strcmp(get(hoc,'style'),'edit')
	% pas de control de valeur (aucun risque d'erreur)
	return
else
	zdata = zuidata(hoc);
	if ischar(zdata)
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
type        = getfield(declaration.type,actionc);
borne       = getfield(declaration.borne,actionc);
defaut      = getfield(declaration.defaut,actionc);

% cas des entiers
if ~isempty(findstr(lower(type),'integer'))|~isempty(findstr(lower(type),'entier'))
	zdata =fix(zdata);
end

% les bornes
if iscell(borne)
	flag = 0;
	for k =1:length(borne)
		if borne{k} == zdata
			flag =1;
		end
	end
	if flag == 0
		zdata =defaut;
	end
else
	if isempty(zdata)
		zdata = defaut;
	elseif zdata < min(borne)
		zdata = min(borne);
	elseif zdata > max(borne)
		zdata =max(borne);
	end
end

zuidata(hoc,zdata,defaut);

%
if strcmp(get(hoc,'tag'),'tdeb')
	% temps de debut -> tdeb
	zdata = zuidata(hoc) ;
	dt = abs(tt-zdata) ;
	kmin = min(find(dt==min(dt))) ;
	zuidata(h.kmin,kmin) ;

elseif strcmp(get(hoc,'tag'),'tfin')
	% temps de fin -> kmax
	zdata = zuidata(hoc) ;
	dt = abs(tt-zdata) ;
	kmax = min(find(dt==min(dt))) ;
	zuidata(h.kmax,kmax) ;
	
elseif strcmp(get(hoc,'tag'),'t')
	% temps -> indice k 
	zdata  = zuidata(hoc) ;
	dt = abs(tt-zdata) ;
	k  = min(find(dt==min(dt))) ;
	zuidata(h.k,k) ;
	
elseif strcmp(get(hoc,'tag'),'k')
	% indice k -> temps
	zdata  = zuidata(hoc) ;
	t  = tt(zdata) ;
	zuidata(h.t,t) ;
	
elseif strcmp(get(hoc,'tag'),'kmin')
	% indice k -> temps
	zdata  = zuidata(hoc) ;
	t  = tt(zdata) ;
	zuidata(h.tdeb,t) ;
	
elseif strcmp(get(hoc,'tag'),'kmax')
	% indice k -> temps
	zdata  = zuidata(hoc) ;
	t  = tt(zdata) ;
	zuidata(h.tfin,t) ;
	
else
	zdata = zuidata(hoc) ;
	
end
