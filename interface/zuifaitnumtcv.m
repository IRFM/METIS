function zuifaitnum(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('accesTCV');
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
			set(h.etat,'string','shot number greater than 0.')
			zuidata(h.numchoc,[],[]);
			return
		else
		
		   racine = zuidata(h.racine);
		   if isempty(racine)
			set(h.etat,'string','the directory must be a path towards TCV data file')
			zuidata(h.racine,'','');
			return
		   elseif exist(racine,'dir') ~= 7
			set(h.etat,'string','the directory must be a valid path towards TCV data file')
			zuidata(h.racine,'','');
			return
		   end
		end
		
		try
		     datatcv  = load(sprintf('%s/%d/tcvtemp',racine,numchoc));
		     zuidata(h.numchoc,numchoc,[]);
	   catch
			try
		      datatcv  = load(sprintf('%s/tcvtemp',racine));
		      zuidata(h.numchoc,numchoc,[]);
			catch
				set(h.etat,'string','no data for this shot number ...')
				zuidata(h.numchoc,[],[]);
				helpdlg('You should use zlecttcv to prepare the data', ...
			        'no data for this shot ...');
				return
			end
		end
		
		chemin = zuidata(h.chemin);
		[s,t]  = unix(sprintf('ls %s',chemin));
		if s ~=0
			set(h.etat,'string','non valid path ...')
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

