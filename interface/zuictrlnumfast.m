function zuicontrolnum(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('accesFAST');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.etat,'string','')

% selon ation
switch lower(action)

	case 'numchoc'
		numchoc = zuidata(h.numchoc);
		if isempty(numchoc)
			set(h.etat,'string','the shot number can not be empty ...')
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
			set(h.etat,'string','false path ...')
			zuidata(h.chemin,'','');
		end
		
	
	otherwise
		warning('action not taken into account')
	   
end

