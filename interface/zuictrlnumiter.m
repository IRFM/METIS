function zuicontrolnum(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('accesITER');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.etat,'string','')

% selon ation
switch lower(action)

	case 'numchoc'
		numchoc = zuidata(h.numchoc);
		if isempty(numchoc)
			set(h.etat,'string','Le numero  du choc ne doit pas etre vide ...')
			zuidata(h.numchoc,1,[]);
		   zuireset(h.validation);
			return
		end
		if numchoc < 0
			set(h.etat,'string','numchoc < 0 ...')
			zuidata(h.numchoc,1,[]);
		   zuireset(h.validation);
			return
		end
		
	case 'chemin'
		chemin = zuidata(h.chemin);
		[s,t]  = unix(sprintf('ls %s',chemin));
		if s ~=0
			set(h.etat,'string','Chemin non valide ...')
			zuidata(h.chemin,'','');
		end
		
	
	otherwise
		warning('ation non prise en compte')
	   
end

