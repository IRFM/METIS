function zuifaitnum(action)

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
			set(h.etat,'string','shot number must be > 0.')
			zuidata(h.numchoc,[],[]);
			return
		else
		
		   racine = zuidata(h.racine);
		   if isempty(racine)
			set(h.etat,'string','the root name must be a path towards JET data')
			zuidata(h.racine,'','');
			return
		   elseif exist(racine,'dir') ~= 7
			set(h.etat,'string','the root name must be a path towards JET data')
			zuidata(h.racine,'','');
			return
		   end
		end
		
		try
		     datajet  = load(sprintf('%s/%d/temp%d',racine,numchoc,numchoc));
		     zuidata(h.numchoc,numchoc,[]);
	   catch
			try
		      datajet  = load(sprintf('%s/temp%d',racine,numchoc));
		      zuidata(h.numchoc,numchoc,[]);
			catch
				set(h.etat,'string','no data for this shot number ...')
				zuidata(h.numchoc,[],[]);
				helpdlg('You must prepare the data (use zjet)', ...
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
		
      % appel de la suite de l'assistant
      evalin('base','zuijetacces2(2);');
      
	otherwise
		warning('action not taken into account')
	   
end

