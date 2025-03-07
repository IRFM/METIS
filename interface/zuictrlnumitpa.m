function zuictrlnumitpa(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('accesITPADB');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.etat,'string','')

% selon action
switch lower(action)

	case 'numchoc'
		numchoc = fix(zuidata(h.numchoc));
		if numchoc < 0
			set(h.etat,'string','Le numero de choc doit etre > 0.')
			zuidata(h.numchoc,[],[]);
		end
		
	case 'chemin'
		chemin = zuidata(h.chemin);
		[s,t]  = unix(sprintf('ls %s',chemin));
		if s ~=0
			set(h.etat,'string','Chemin non valide ...')
			zuidata(h.chemin,'','');
		end
		
	
	otherwise
		warning('action non prise en compte')
	   
end

