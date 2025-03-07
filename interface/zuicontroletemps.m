function zuicontroletemps(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('tempsTS');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.commentaire,'string','');
tdeb=getappdata(hfig,'tdeb');
tfin=getappdata(hfig,'tfin');

vdeb  = zuidata(h.valeur_xdeb);
vfin  = zuidata(h.valeur_xfin);

% selon ation
switch lower(action)

case 'valeur_xdeb'
	if vdeb < tdeb
		zuidata(h.commentaire,'temps de debut trop petit ...');
		vdeb = tdeb;
	elseif vdeb > vfin
		zuidata(h.commentaire,'temps de debut plus grand que le temps de fin ...');
		vdeb = tdeb;
	end	
case 'valeur_xfin'
	if tfin < vfin
		zuidata(h.commentaire,'temps de fin trop grand ...');
		vfin = tfin;
	elseif vdeb > vfin
		zuidata(h.commentaire,'temps de fin plus petit que le temps de debut ...');
		vfin = tfin;
	end	
	
otherwise
		warning('ation non prise en compte')
	   
end

zuidata(h.valeur_xdeb,vdeb);
zuidata(h.valeur_xfin,vfin);
