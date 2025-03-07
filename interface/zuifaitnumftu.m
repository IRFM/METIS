function zuifaitnum(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('accesftu');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.etat,'string','')

% selon ation
switch lower(action)

	case {'annulation','close'}
		zuicloseone(hfig);	
	
	case {'init','raz'}
		zuiformvisible(hfig);
		zuiformreset(hfig);
		zuiuploadform(hfig);
		zuireset(h.raz);
	
	case 'validation'
		numchoc = fix(zuidata(h.numchoc));
		if numchoc < 0
			set(h.etat,'string','the shot number have to be > 0.')
			zuidata(h.numchoc,[],[]);
			return
		else
		
		   racine = zuidata(h.racine);
		   if isempty(racine)
			set(h.etat,'string','the path should pointed towards FTU data')
			zuidata(h.racine,'','');
			return
		   elseif exist(racine,'dir') ~= 7
			set(h.etat,'string','the path should pointed towards FTU data')
			zuidata(h.racine,'','');
			return
		   end
		end
		
		try
		     dataftu  = load(sprintf('%s/%d/ftutemp',racine,numchoc));
		     zuidata(h.numchoc,numchoc,[]);
	   catch
			try
		      dataftu  = load(sprintf('%s/ftutemp',racine));
		      zuidata(h.numchoc,numchoc,[]);
			catch
				set(h.etat,'string','No data for this shot ...')
				zuidata(h.numchoc,[],[]);
				helpdlg('you have to use zlectftu to generate the datas ', ...
			        'No data for this shot ...');
				return
			end
		end
		
		chemin = zuidata(h.chemin);
		[s,t]  = unix(sprintf('ls %s',chemin));
		if s ~=0
			set(h.etat,'string','unknown path ...')
			zuidata(h.chemin,'','');
		        zuireset(h.validation);
			return
		end

		zuiformcache(hfig) ;
		zuidownloadform(hfig);
		zuicloseone(hfig);	
		
	otherwise
		warning('action not taken into account')
	   
end

