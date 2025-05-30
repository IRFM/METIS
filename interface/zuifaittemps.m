function zuifaitnum(action)

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
tfin=getappdata(hfig,'tfin');

vdeb  = zuidata(h.valeur_xdeb);
vfin  = zuidata(h.valeur_xfin);


% selon ation
switch lower(action)
	
	case 'select_xdeb'
		[x,y] = ginput(1);
		if x < vfin
			zuidata(h.valeur_xdeb,max(x,tdeb));
		else
			zuidata(h.commentaire,'temps de debut plus grand que le temps de fin ...');
		end
		zuireset(h.select_xdeb);
		return
	case 'select_xfin'	
		[x,y] = ginput(1);
		if x > vdeb
			zuidata(h.valeur_xfin,min(x,tfin));
		else
			zuidata(h.commentaire,'temps de fin plus petit que le temps de debut ...');
		end
		zuireset(h.select_xfin);
		return
	case {'annulation','close'}
		zuicloseone(hfig);	
	
	case {'init','raz'}
		zuiformvisible(hfig);
		zuiformreset(hfig);
		zuireset(h.raz);
		zuidata(h.valeur_xdeb,vdeb);
		zuidata(h.valeur_xfin,vfin);
		setappdata(hfig,'tdeb',vdeb);
		setappdata(hfig,'tfin',vfin);

	case 'validation'
		if vdeb < tdeb
			zuidata(h.commentaire,'temps de debut trop petit ...');
			vdeb = tdeb;
			zuidata(h.valeur_xdeb,vdeb);
			zuidata(h.valeur_xfin,vfin);
			zuireset(h.validation);
			return
		elseif vdeb > vfin
			zuidata(h.commentaire,'temps de debut plus grand que le temps de fin ...');
			vdeb = tdeb;
			zuidata(h.valeur_xdeb,vdeb);
			zuidata(h.valeur_xfin,vfin);
			zuireset(h.validation);
			return
		end	
		if tfin < vfin
			zuidata(h.commentaire,'temps de fin trop grand ...');
			vfin = tfin;
			zuidata(h.valeur_xdeb,vdeb);
			zuidata(h.valeur_xfin,vfin);
			zuireset(h.validation);
			return
		elseif vdeb > vfin
			zuidata(h.commentaire,'temps de fin plus petit que le temps de debut ...');
			vfin = tfin;
			zuidata(h.valeur_xdeb,vdeb);
			zuidata(h.valeur_xfin,vfin);
			zuireset(h.validation);
			return
		end	

		zuiformcache(hfig) ;
		zuidownloadform(hfig);
		zuicloseone(hfig);	
		
	otherwise
		warning('ation non prise en compte')
	   
end

