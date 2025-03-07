function zuifaitnum(action)

if nargin ==0
	action = 'init';
elseif isempty(action)
	action = 'init';
end


% recupere le handle de la fenetre concernee
[hfig,h] = zuiformhandle('tempsJET');
% information pour l'assistant
zuicr(hfig,action) ;
set(h.commentaire,'string','')

tdeb=getappdata(hfig,'tdeb');
tfin=getappdata(hfig,'tfin');

vdeb  = zuidata(h.valeur_xdeb);
vfin  = zuidata(h.valeur_xfin);


% selon ation
switch lower(action)
	
	case 'select_xdeb'
		ind  = zuidata(h.select_listediag);
		if ind > 2
			zuidata(h.commentaire,'initial time is a MSE or Pol time...');
		        zuireset(h.select_xdeb);
		        return
		end
		[x,y] = ginput(1);
		if x < vfin
			zuidata(h.valeur_xdeb,max(x,tdeb));
		else
			zuidata(h.commentaire,'initial time greater than final time ...');
		end
		zuireset(h.select_xdeb);
		return
	case 'select_xfin'	
		[x,y] = ginput(1);
		if x > vdeb
			zuidata(h.valeur_xfin,min(x,tfin));
		else
			zuidata(h.commentaire,'final time smallest than the initial time ...');
		end
		zuireset(h.select_xfin);
		return
		
	case {'annulation','close'}
		zuicloseone(hfig);	
	
	case {'init','raz'}
		zuiformvisible(hfig);
		zuiformreset(hfig);
		zuireset(h.raz);
		zuidata(h.valeur_xdeb,vdeb);
		zuidata(h.valeur_xfin,vfin);
		setappdata(hfig,'tdeb',vdeb);
		setappdata(hfig,'tfin',vfin);

	case 'validation'
		if vdeb < tdeb
			zuidata(h.commentaire,'initial time too small ...');
			vdeb = tdeb;
			zuidata(h.valeur_xdeb,vdeb);
			zuidata(h.valeur_xfin,vfin);
			zuireset(h.validation);
			return
		elseif vdeb > vfin
			zuidata(h.commentaire,'initial time greater than final time ...');
			vdeb = tdeb;
			zuidata(h.valeur_xdeb,vdeb);
			zuidata(h.valeur_xfin,vfin);
			zuireset(h.validation);
			return
		end	
		if tfin < vfin
			zuidata(h.commentaire,'final time too high ...');
			vfin = tfin;
			zuidata(h.valeur_xdeb,vdeb);
			zuidata(h.valeur_xfin,vfin);
			zuireset(h.validation);
			return
		elseif vdeb > vfin
			zuidata(h.commentaire,'final time smallest than the initial time ...');
			vfin = tfin;
			zuidata(h.valeur_xdeb,vdeb);
			zuidata(h.valeur_xfin,vfin);
			zuireset(h.validation);
			return
		end	
      vdia  = zuidata(h.select_listediag);
      switch vdia
	   case 1
		       zassignin('base','prepare.tmse',NaN);
		       zassignin('base','prepare.tpol',NaN);
             zassignin('base','prepare.mse',-1);	   
             zassignin('base','prepare.pol',-1);	   
      case 2
		       zassignin('base','prepare.tmse',NaN);
		       zassignin('base','prepare.tpol',NaN);
             zassignin('base','prepare.mse',0);	   
             zassignin('base','prepare.pol',0);	   
	   case 3
	          tmse = zuidata(h.select_listetemps);
		       zassignin('base','prepare.tmse',tmse);
		       zassignin('base','prepare.tpol',NaN);
             zassignin('base','prepare.mse',1);	   
             zassignin('base','prepare.pol',0);	   
	   case 4
		       zassignin('base','prepare.tmse',NaN);
	          tpol = zuidata(h.select_listetemps);
		       zassignin('base','prepare.tpol',tpol);
             zassignin('base','prepare.mse',0);	   
             zassignin('base','prepare.pol',1);	   
	   end
		zuiformcache(hfig) ;
		zuidownloadform(hfig);
		%zuicloseone(hfig);	
		 
      
      evalin('base','zuijetacces2(3);');
      
	otherwise
		warning('action not taken into account')
	   
end

