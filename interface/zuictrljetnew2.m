function zuicontroletemps(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('tempsJET');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.commentaire,'string','');
tdeb    = getappdata(hfig,'tdeb');
tfin    = getappdata(hfig,'tfin');
tmse    = getappdata(hfig,'tmse');
tpol    = getappdata(hfig,'tpol');

vdeb    = zuidata(h.valeur_xdeb);
vfin    = zuidata(h.valeur_xfin);
vdia    = zuidata(h.select_listediag);
vtemps  = zuidata(h.select_listetemps);
% selon ation
switch lower(action)

case 'valeur_xdeb'
	if vdeb < tdeb
		zuidata(h.commentaire,'first time too small  ...');
		vdeb = tdeb;
	elseif vdeb > vfin
		zuidata(h.commentaire,'first time greater than final time...');
		vdeb = tdeb;
	end	
case 'valeur_xfin'
	if tfin < vfin
		zuidata(h.commentaire,'final time too large ...');
		vfin = tfin;
	elseif vdeb > vfin
		zuidata(h.commentaire,'final time slower than first time ...');
		vfin = tfin;
	end
case 'select_listediag'
  switch vdia
    case 1
      zuicache(h.select_listetemps)
      zuienable(h.select_xdeb);
	   zuienable(h.valeur_xdeb);

    case 2
      zuicache(h.select_listetemps)
      zuienable(h.select_xdeb);
	   zuienable(h.valeur_xdeb);

    case 3
      zuidisable(h.select_xdeb)
	   zuidisable(h.valeur_xdeb);
      sliste =sprintf('|%d',1:length(tmse));
	   sliste(1)=[];
      set(h.select_listetemps,'string',sliste,'value',1,'userdata',tmse, ...
	       'tooltip','time for MSE reconstruction  (intitial current profile)', ...
		    'visible','on');
      vdeb  = tmse(1);

    case 4
      zuidisable(h.select_xdeb)
	   zuidisable(h.valeur_xdeb);
      sliste =sprintf('|%d',1:length(tpol));
	   sliste(1)=[];
      set(h.select_listetemps,'string',sliste,'value',1,'userdata',tpol, ...
	       'tooltip','time for Polarimetry reconstruction  (intitial current profile)', ...
		    'visible','on');
      vdeb  = tpol(1);
   end
case 'select_listetemps'
	   
       vdeb = zuidata(h.select_listetemps);
		
otherwise

       warning('action not taken into account')

end


zuidata(h.valeur_xdeb,vdeb);
zuidata(h.valeur_xfin,vfin);
