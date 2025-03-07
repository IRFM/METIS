function zuicontrolnum(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('accesJET');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.etat,'string','')

% selon ation
switch lower(action)

	case 'numchoc'
		numchoc = fix(zuidata(h.numchoc));
		if numchoc < 0
			set(h.etat,'string','the shot number must be > 0.')
			zuidata(h.numchoc,[],[]);
		else
		   try
			  racine = zuidata(h.racine);
		     datajet  = load(sprintf('%s/%d/temp%d',racine,numchoc,numchoc));
		     zuidata(h.numchoc,numchoc,[]);
		   catch
			  try
		     		datajet  = load(sprintf('%s/temp%d',racine,numchoc));
		     		zuidata(h.numchoc,numchoc,[]);
			  catch
					set(h.etat,'string','no data for this shot ...')
					zuidata(h.numchoc,[],[]);
				end
		   end
		end
		
	case 'chemin'
		chemin = zuidata(h.chemin);
		[s,t]  = unix(sprintf('ls %s',chemin));
		if s ~=0
			set(h.etat,'string','non valid path ...')
			zuidata(h.chemin,'','');
		end
		
	case 'racine'
		racine = zuidata(h.racine);
		[s,t]  = unix(sprintf('ls %s',racine));
		if s ~=0
			set(h.etat,'string','Chemin non valide ...')
			zuidata(h.racine,'','');
		end
	
	otherwise
		warning('action not taken into account')
	   
end

