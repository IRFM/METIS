function zuicontrolnum(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('accesTS');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.etat,'string','')

% selon ation
switch lower(action)

	case 'numchoc'
		numchoc = zuidata(h.numchoc);
		try
		    [text,void1,void2,void3]    = tsbase(numchoc,'scprof');
		catch
		    text='';
		end
		if isempty(text)
			set(h.etat,'string','Pas de donnees pour ce choc et cette occurence ...')
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
		warning('ation non prise en compte')
	   
end

