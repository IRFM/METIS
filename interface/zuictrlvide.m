function zuicontrolnum(action)
%
%
%
% fonction ecrite par C. Passeron, poste 6119
% version  2.0 , du 11/12/2002.
% 
% liste des modifications : 
% * 11/12/2002 : interface en anglais
%--------------------------------------------------------------

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('accesvide');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.etat,'string','')

numchoc  = fix(zuidata(h.numchoc));
chemin   = zuidata(h.chemin);
tdeb     = zuidata(h.tdeb);
tfin     = zuidata(h.tfin);
pas      = zuidata(h.pas);
machine  = zuidata(h.machine);


% selon ation
switch lower(action)

case 'numchoc'
	if isempty(numchoc)
		zuidata(h.numchoc,getappdata(hfig,'pid'),[]);
		set(h.etat,'string','You must give a shot number')
	elseif numchoc < 0
		set(h.etat,'string','The shot number must be > 0.')
		zuidata(h.numchoc,getappdata(hfig,'pid'),[]);
	end
	
case 'chemin'
	[s,t]  = unix(sprintf('ls %s',chemin));
	if s ~=0
		set(h.etat,'string','non valid path ...')
		zuidata(h.chemin,'','');
	end
	
case 'tdeb'
	if isempty(tdeb)
		set(h.etat,'string','You must give an initial time')
	elseif tfin < tdeb
		set(h.etat,'string','Initial time must be smallest than the final time')
		zuidata(h.tdeb,0,[]);
	end
	
case 'tfin'
	if isempty(tfin)
		set(h.etat,'string','You must give a final time')
	elseif tfin < tdeb
		set(h.etat,'string','final time must be greatest than the initial time')
		zuidata(h.tfin,tdeb + 30 ,[]);
	end
	
case 'pas'
	if isempty(pas)
		set(h.etat,'string','You must give dt')
	elseif pas > (abs(tfin - tdeb)/2)
		set(h.etat,'string','dt too high ( > (tfin-tdeb)/2)')
		zuidata(h.pas,abs(tfin - tdeb)/2,[]);
	end
	
case 'machine'
	if isempty(pas)
		set(h.etat,'string','Tokamak name ?')
		zuidata(h.machine,'DEMO',[]);
	end
	
otherwise
	warning('action not taken into account')
	
end

