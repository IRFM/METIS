function zuicontroletemps(action)
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
[hfig,h] = zuiformhandle('tempsJET');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.commentaire,'string','');
tdeb=getappdata(hfig,'tdeb');
tfin=getappdata(hfig,'tfin');
tmse=getappdata(hfig,'tmse');
tpol=getappdata(hfig,'tpol');

vdeb  = zuidata(h.valeur_xdeb);
vfin  = zuidata(h.valeur_xfin);
vdia  = zuidata(h.select_listediag);
% selon ation
switch lower(action)

case 'valeur_xdeb'
	if vdeb < tdeb
		zuidata(h.commentaire,'initial time too small ...');
		vdeb = tdeb;
	elseif vdeb > vfin
		zuidata(h.commentaire,'initial time greater than final time ...');
		vdeb = tdeb;
	end	
case 'valeur_xfin'
	if tfin < vfin
		zuidata(h.commentaire,'final time too high ...');
		vfin = tfin;
	elseif vdeb > vfin
		zuidata(h.commentaire,'final time smaller than initial time ...');
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
		     'tooltip','MSE profile for the initial plasma current', ...
			  'visible','on');
	  case 4
	    zuidisable(h.select_xdeb)
		 zuidisable(h.valeur_xdeb);
	    sliste =sprintf('|%d',1:length(tpol));
		 sliste(1)=[];
	    set(h.select_listetemps,'string',sliste,'value',1,'userdata',tpol, ...
		     'tooltip','Polarimetry profile for the initial plasma current', ...
			  'visible','on');
	  
	end
case 'select_listetemps'
	vdeb = zuidata(h.select_listetemps);
		
otherwise
		warning('action not taken into account')
	   
end
zuidata(h.valeur_xdeb,vdeb);
zuidata(h.valeur_xfin,vfin);
