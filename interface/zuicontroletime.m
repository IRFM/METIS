
function zuicontroletime(action)

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

vdeb  = zuidata(h.valeur_xdeb);

% selon ation
switch lower(action)

case 'valeur_xdeb'
	if vdeb < tdeb
		zuidata(h.commentaire,'two small input time ...');
		vdeb = tdeb;
	end	
	
otherwise
		warning('action not taken into account')
	   
end

zuidata(h.valeur_xdeb,vdeb);