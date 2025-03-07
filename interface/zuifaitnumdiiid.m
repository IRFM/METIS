function zuifaitnum(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('accesDIIID');
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
			set(h.etat,'string','the shot bnumber must be > 0.')
			zuidata(h.numchoc,[],[]);
			return
		else
		
		   racine = zuidata(h.racine);
		   if isempty(racine)
			set(h.etat,'string','wrong path for data access')
			zuidata(h.racine,'','');
			return
		   elseif exist(racine,'dir') ~= 7
			set(h.etat,'string','wrong path for DIIID data access')
			zuidata(h.racine,'','');
			return
		   end
		end
		
		try
		     datadiiid  = load(sprintf('%s/%d/diiidtemp',racine,numchoc));
		     zuidata(h.numchoc,numchoc,[]);
	   catch
			try
		      datadiiid  = load(sprintf('%s/diiidtemp',racine));
		      zuidata(h.numchoc,numchoc,[]);
			catch
				set(h.etat,'string','No data for this shot ...')
				zuidata(h.numchoc,[],[]);
				helpdlg('you should use zlectdiiid to prepare the input CRONOS data', ...
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

