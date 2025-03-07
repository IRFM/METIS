function zuifaittempsluke(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('tempsTS');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.commentaire,'string','')

tdeb=getappdata(hfig,'tdeb');

vdeb  = zuidata(h.valeur_xdeb);

% selon ation
switch lower(action)
	
	case 'select_xdeb'
		[x,y] = ginput(1);

	    zuidata(h.valeur_xdeb,max(x,tdeb));
		zuireset(h.select_xdeb);
		return
	case {'annulation','close'}
		zuicloseone(hfig);	
	
	case {'init','raz'}
		zuiformvisible(hfig);
		zuiformreset(hfig);
		zuireset(h.raz);
		zuidata(h.valeur_xdeb,vdeb);
		setappdata(hfig,'tdeb',vdeb);
	case 'validation'
		if vdeb < tdeb
			zuidata(h.commentaire,'too small input time ...');
			vdeb = tdeb;
			zuidata(h.valeur_xdeb,vdeb);
			zuireset(h.validation);
			return
		end	
		zuiformcache(hfig) ;
		zuidownloadform(hfig);
		zuicloseone(hfig);	
		
	otherwise
		warning('action not taken into account')
	   
end
