% ZDATAPLOT permet d'afficher et de tracer les donnees zineb
%-----------------------------------------------------------
% fichier zdataplot.m ->  zdataplot
%
%
% fonction Matlab 5 :
%
% Cette fonction genere une figure avec des menus
% permettant d'afficher et de tracer les donnees de zineb.
% Les donnees doivent etre dans l'espace de travail de base.
% Elles sont contenue dans les structure data et param
%
% Les donnees peuvent etre exportees depuis n'importe quel
% point d'arret dans une fonction vers l'espace de travail de
% base a l'aide de la commande 'zexport'.
%
% syntaxe  :
%
%     zdataplot;
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.7, du 22/10/2001.
%
%
% liste des modifications :
%
%   * 08/06/2001 -> affichage des structures
%   * 08/06/2001 -> correction du bug sur 1temps
%   * 08/06/2001 -> modification du style pour 1 et 2 temps
%   * 21/09/2001 -> ajout de la parois
%   * 24/10/2001 -> ajout de l'option de recherche
%   * 21/11/2001 -> ajout des clipboards
%   * 22/11/2001 -> ajout des fonctions pour selectionner un temps
%                   et supperposer des profils
%   * 29/11/2004 -> english version by default (sorry guys !)
%
%--------------------------------------------------------------
%
function zdataplot_metis(cmd,hcb)


% gestion des entrees
if nargin <1
	cmd = 'init';
elseif isempty(cmd)
	cmd = 'init';
end

% mode initialisation
if strcmp(cmd,'init')

	% suppresion du debgug sur error
	status_debug = dbstatus;
	if ~isempty(status_debug)
	      if ~isempty({status_debug.cond})
		      loc_cond = status_debug.cond;
		      if strmatch(loc_cond,'error')
			      dbclear error
		      end
	      end
	end
	% ouverture des clipboard
	%zgclipp;
	%zgclipt;

	% creation de la figure
	hf=figure('tag','zdataplot_metis','toolbar','figure','name','zdataplot_metis');
	set(hf,'defaultaxesfontsize',12,'units','normalized','position',[0.1 0.1 0.7 0.7], ...
	    'color',[1 1 1],'DeleteFcn','clear_zdataplot_imas_menu');

 	% lecture des informations
        try
            zs  = evalin('base','post');
        catch
            zs = [];
        end
	if isempty(zs)
	     close(hf);
	     warndlg('No METIS data load (needed as template)','METIS data browser');
	     return
	elseif ~isfield(zs,'z0dinput')
                return
        else
		try
	        	zsinfo = zs.z0dinput.zsinfo;
		catch
			zsinfo = zero1t;
		end
		try
	        	profinfo = zs.z0dinput.profinfo;
		catch
			profinfo = z0dprofinfo;
		end
	end

        % merge des structures
        zsinfo2 = zero1t;
        noms = fieldnames(zsinfo2);
        for k = 1:length(noms)
	      if ~isfield(zsinfo,noms{k})
		  zsinfo.(noms{k}) = zsinfo2.(noms{k});
	      end
	end
 	profinfo2 = z0dprofinfo;
        noms = fieldnames(profinfo2);
        for k = 1:length(noms)
	      if ~isfield(profinfo,noms{k})
		  profinfo.(noms{k}) = profinfo2.(noms{k});
	      end
	end


	fprintf('Creation of menus , please be patient ...\n');

	% creation du uicontrol num choc
	umem =get(hf,'units');
	set(hf,'units','pixels');
	pos =get(hf,'position');
	set(hf,'units',umem);



	% menu
	hracine  = uimenu(hf,'label','METIS 0D','tag','0Dmenu');
	cree0d(hracine,zsinfo)
	hracine  = uimenu(hf,'label','EXP 0D','tag','exp0Dmenu');
	creeexp0d(hracine,zsinfo)
	hracine  = uimenu(hf,'label','METIS Profiles','tag','0Dmenu1D');
	cree0d1d(hracine,profinfo)

	fprintf('\n');
	%delete(hcb);
  


    % substituion des labels
    hso = findobj_local(hf,'type','uimenu');
    for k =1:length(hso)
	lab = get(hso(k),'label');
	set(hso(k),'label',ztr(lab));
    end

    % les axes
    h=subplot(4,1,1);
    set(h,'tag','axes_des_temps');
    xlabel('Time (s)')
    % le bouton pour selectionner le temps
    poscc = get(h,'position');
    uicontrol(hf,'style','checkbox', ...
                 'units','normalized', ...
                 'position',[0,poscc(2)+0.5*poscc(4),0.02,0.02], ...
		 'callback','zdataplot_metis(''pick_t'');', ...
		 'tooltip','Select time slice with the mouse', ...
		 'tag','pick_t');
    % suite des axes
    h=subplot(4,1,2);
    set(h,'tag','axes_des_phases','xticklabel',[]);
    %xlabel('temps (s)')
    %title('angle(f(t))')
    h=subplot(2,3,4);
    set(h,'tag','axes_des_profils');
    if verLessThan('matlab', '8.5')
        xlabel('r','Fontname','Symbol')
    else
        xlabel('$$\rho$$','Interpreter','latex');
    end
    title('Profiles')
    % le bouton pour selectionner le temps
    poscc = get(h,'position');
    uicontrol(hf,'style','checkbox', ...
                 'units','normalized', ...
                 'position',[0,poscc(2)+0.5*poscc(4),0.02,0.02], ...
		 'callback','zdataplot_metis(''supperpose'');', ...
		 'tooltip','Hold on time slice', ...
		 'tag','supperpose');
    % suite des axes
    h=subplot(2,3,5);
    set(h,'tag','axes_des_images');
    if verLessThan('matlab', '8.5')
        xlabel('r','Fontname','Symbol')
    else
        xlabel('$$\rho$$','Interpreter','latex');
    end
    title('2D Profiles')
    h=subplot(2,3,6);
    set(h,'tag','axes_des_surfaces');
    xlabel('R (m)')
    ylabel('Z (m)')
    title('Equilibrium')

    % le menu des commandes
    h=uimenu(hf,'label','Commands','tag','commande');
    uimenu(h,'label','Plot also reference', ...
              'Callback','zdataplot_metis(''plot_ref'')', ...
              'tag','plot_ref');
    uimenu(h,'label','Hold on time traces', ...
              'Callback','zdataplot_metis(''efface_temps'')', ...
              'tag','efface_temps');
    uimenu(h,'label','Hold on profiles', ...
              'Callback','zdataplot_metis(''efface_profils'')', ...
              'tag','efface_profil');
    uimenu(h,'label','LogY time traces', ...
              'Callback','zdataplot_metis(''log_temps'')', ...
              'tag','log_temps');
    uimenu(h,'label','LogY profiles', ...
              'Callback','zdataplot_metis(''log_profils'')', ...
              'tag','log_profil');
    uimenu(h,'label','Log10 2D image', ...
              'Callback','zdataplot_metis(''log_image'')', ...
              'tag','log_image');

    uimenu(h,'label','Plot equilibrium', ...
              'Callback','zdataplot_metis(''plot_equi'')', ...
              'tag','plot_equi');

    uimenu(h,'label','Look for NaN or Inf', ...
              'Callback','zdataplot_metis(''NaNInf'')', ...
              'tag','NaNInf');

    uimenu(h,'label','Extract axes', ...
              'Callback','zdataplot_metis(''extrait'')', ...
              'tag','extrait');

    uimenu(h,'label','Look for variable', ...
              'Callback','zdataplot_metis(''cherche'')', ...
              'tag','cherche');

    uimenu(h,'label','Access to WEST and Tore Supra data', ...
              'Callback','zdataplot_metis(''creeTS'')', ...
              'tag','creeTS');

    uimenu(h,'label','Access to IMAS data', ...
              'Callback','zdataplot_metis(''imas'')', ...
              'tag','cree_imas');

    uimenu(h,'label','Access to IMAS data (simplified menu)', ...
              'Callback','zdataplot_metis(''imas_simple'')', ...
              'tag','cree_imas_simple');

    uimenu(h,'label','Clear IMAS data cache', ...
              'Callback','zdataplot_metis(''clear_imas'')', ...
              'tag','clear_imas');

    uimenu(h,'label','Set IMAS environment', ...
              'Callback','imasdb', ...
              'tag','imasdb');

     uimenu(h,'label','Export IMAS data to matfile', ...
              'Callback','zdataplot_metis(''export_imas'')', ...
              'tag','export_imas');

   uimenu(h,'label','Open Clipboard', ...
              'Callback','zdataplot_metis(''clipboard'')', ...
              'tag','open clipboards for axes');

    uimenu(h,'label','Connect to JET SAL data server', ...
              'Callback','sal_auth;', ...
              'tag','clear_imas');
          
   uimenu(h,'label','JET data', ...
              'Callback','zdataplot_metis(''jet_sal_list'')', ...
              'tag','access to JET data');
     
   uimenu(h,'label','Load CRONOS data', ...
              'Callback','zdataplot_metis(''load_cronos'')', ...
              'tag','Load a CRONOS simulation data file');
          
   uimenu(h,'label','Open CRONOS data browser', ...
              'Callback','zdataplot_metis(''cronos_data_browser'')', ...
              'tag','Open CRONOS dedicated data browser');
          
   % restauration du debgug sur error
   %status_debug = dbstatus;
   if ~isempty(status_debug)
       if ~isempty({status_debug.cond})
           if length(status_debug) == 1
               if strmatch(status_debug.cond,'error')
                   dbstop if error
               end
           else
               for klm =1:length(status_debug)
                   if strmatch(status_debug(klm).cond,'error')
                       dbstop if error
                   end                  
               end
               
           end
       end
   end


% fin de init
end


% objet courant
if nargin < 2
    hcb= gcbo;
    hf=gcf;
    ho=gco;
else
    ho = hcb;
    hf = findobj_local(0,'tag','zdataplot_metis','type','figure');
    if ~isempty(hf)
        hf=hf(1);
    end
end

compare = strcmp(get(findobj_local(hf,'type','uimenu','tag','plot_ref'),'checked'),'on');

if isempty(ho)&isempty(hcb)
	return
end




if strcmp(cmd,'struct')
	nom  = strrep(get(hcb,'tag'),'@','');
	aide = get(hcb,'userdata');
	try
	   nomv = nom;
	   data=evalin('base',nom);
	catch
	    nomv = nom;
	    data =[];
	end
	if isempty(data)
		data='<Empty>';
	end

	disp('-------------------------------------------')
	fprintf('Name of the structure : %s\n',nom);
	disp('contents :')
	disp(' ')
	disp(data)
	disp(' ')
end


% log axes des temps
if strcmp(cmd,'log_temps')
	h=findobj_local(hf,'type','axes','tag','axes_des_temps');
	if ~isempty(h)
		yscale= lower(get(h,'yscale'));
		if strcmp(yscale,'log')
			set(h,'yscale','linear');
			set(hcb,'checked','off');
		else
			set(h,'yscale','log');
			set(hcb,'checked','on');
		end
	end

else
	h=findobj_local(hf,'type','axes','tag','axes_des_temps');
 	hm=findobj_local(hf,'type','uimenu','tag','log_temps');
	if ~isempty(hm)
	   chk = get(hm,'checked');
	   if strcmp(chk,'on')
	   	set(h,'yscale','log');
	   else
		set(h,'yscale','linear');
	   end
        end
end


% log axes des profils
if strcmp(cmd,'log_profils')
	h=findobj_local(hf,'type','axes','tag','axes_des_profils');
	if ~isempty(h)
		yscale= lower(get(h,'yscale'));
		if strcmp(yscale,'log')
			set(h,'yscale','linear');
			set(hcb,'checked','off');
		else
			set(h,'yscale','log');
			set(hcb,'checked','on');
		end
	end

else
	h=findobj_local(hf,'type','axes','tag','axes_des_profils');
 	hm=findobj_local(hf,'type','uimenu','tag','log_profil');
	if ~isempty(hm)
	   chk = get(hm,'checked');
	   if strcmp(chk,'on')
	   	set(h,'yscale','log');
	   else
		set(h,'yscale','linear');
	   end
        end
end

% log axes des profils
if strcmp(cmd,'log_image')
	chk = get(hcb,'checked');
	if strcmp(chk,'on')
		set(hcb,'checked','off');
	else
		set(hcb,'checked','on');
	end
end



% hold axes des temps
if strcmp(cmd,'efface_temps')
	yscale= lower(get(hcb,'checked'));
	if strcmp(yscale,'on')
		set(hcb,'checked','off')
	else
		set(hcb,'checked','on')
	end
end

% hold axes des temps
if strcmp(cmd,'efface_profils')
	yscale= lower(get(hcb,'checked'));
	if strcmp(yscale,'on')
		set(hcb,'checked','off')
	else
		set(hcb,'checked','on')
	end
end

% compare
if strcmp(cmd,'plot_ref')
	yscale= lower(get(hcb,'checked'));
	if strcmp(yscale,'on')
		set(hcb,'checked','off')
	else
		set(hcb,'checked','on')
	end
end

% plot de l'equilibre
if strcmp(cmd,'plot_equi')
	%
	ok   = 1;
	R1   = [];
	mode = 0;
	try
	    post = evalin('base','post');
	    [temps,R,Z,Rext,Zext]=z0make2dsurface(post);
	catch
	    ok=0;
	end
	if ok == 1
	  	if compare
	   		try
			      post = evalin('base','jeux1.post');
			      [temps1,R1,Z1,Rext1,Zext1]=z0make2dsurface(post);
	   		end
   	        end
        end
	if isempty(R)
		R=[];
		Z=[];
		Rext=[];
		Zext=[];
		temps=[];
		R1=[];
	end

	disp('-------------------------------------------')
	fprintf('Equilibrium display \n');
	disp('time dependent')

	fprintf('data size :')
	disp(size(R))
	disp(' ')

	% les axes
	h=subplot(2,3,6);
	hold off
	plot(NaN,NaN);
	set(h,'tag','axes_des_surfaces');

	% limitationdu nombre de surfaces
	pas=1;
	indr = 1:pas:size(R,2);
	indr1 = 1:pas:size(R1,2);

	for kl =1:length(indr)
		rr= squeeze(R(:,indr(kl),:));
		zz= squeeze(Z(:,indr(kl),:));
		if ~isempty(R1)
			rr1= squeeze(R1(:,indr1(kl),:));
			zz1= squeeze(Z1(:,indr1(kl),:));
		else
			rr1=[];
			zz1=[];
		end

		if mode  == 0
			zplotprof(h,temps,rr,zz,'color',[0 0.7 1]);
			if ~isempty(R1)
			   zplotprof(h,temps1,rr1,zz1,'color',[0 0.7 1],'linestyle',':');
			end
		else
			plot(rr(1,:),zz(1,:),'color',[0 0.7 1]);
			hold on
			plot(rr(2,:),zz(2,:),'m:');
		end
		fprintf('.');
	end
	if mode  == 0
		zplotprof(h,temps,Rext,Zext,'color','r');
		if ~isempty(R1)
			zplotprof(h,temps1,Rext1,Zext1,'color','r','linestyle',':');
		end
	else
		plot(Rext(1,:),Zext(1,:),'r');
		plot(Rext(2,:),Zext(2,:),':r');
		hold on
	end
	% plot de la paroi
	try
	    % limiter
	    rwall = [];
	    zwall = [];
	    if isfield(post.z0dinput.option,'first_wall')
		[pfw,ffw] = fileparts(post.z0dinput.option.first_wall);
		filename_wall = fullfile(pfw,sprintf('%s.mat',ffw));
		if exist(filename_wall,'file')
			wall = load(filename_wall);
			if isfield(wall,'R') && isfield(wall,'Z')
			    rwall = wall.R(:);
			    zwall = wall.Z(:);
			end
		end
        end
	    if isempty(rwall)
		switch post.z0dinput.machine
		case 'ITER'
		    if  exist('iterwall')
			[rwall,zwall] = iterwall;
		    else
			rwall = [];
			zwall = [];
		    end    
		otherwise
		    if ~isempty(strfind(upper(post.z0dinput.machine),'WEST')) 
			    [rwall,zwall] = west_limiter(post.z0dinput.shot);
			    try
			      walldata = Get_Paroi_WEST(post.z0dinput.shot,temps);
			      rwall    = walldata.Rparoi;
			      zwall    = walldata.Zparoi;
                  rwall(:,end+1) = rwall(:,1);
                  zwall(:,end+1) = zwall(:,1);                  
			    end
			      
		    elseif ~isempty(strfind(post.z0dinput.machine,'JT-60SA')) || ~isempty(strfind(post.z0dinput.machine,'JT60-SA'))
			    wall = load('jt60sawall');
			    rwall = wall.wall.data(1:end,1);
			    zwall = wall.wall.data(1:end,2);
		    else
			rwall = [];
			zwall = [];
		    end
		end
	    end
	   if ~isempty(rwall) && ~isempty(zwall)
	        hold on
            if all(size(rwall) > 1) && (size(rwall,1) == length(temps))
                zplotprof(h,temps,rwall,zwall,'color','k','linestyle','-');
            else
                plot(rwall,zwall,'k')
            end
	   end

	end
	fprintf('\n');
	hold off
	xlabel('R (m)')
	ylabel('Z (m)')
	axis('equal');
	title('Equilibrium (isoflux)')
	set(h,'tag','axes_des_surfaces');

end


if strcmp(cmd,'NaNInf')
	%
	ok   = 1;
	mode = 0;
	post  = evalin('base','post');
	if ~isempty(post)
	   % 2 - la structure data
	   fprintf('---------------------------------\n');
	   fprintf('Search of NaN and Inf :\n');
	   % liste des champs la structures
	   champ = sort(fieldnames(post));
	   for k=1:length(champ)
		   % nom complet pour acceder a la variable
		   champ{k}=strcat('post.',champ{k});
	   end

	   % jusqu'a ce qu'il n'y ait plus de champ
	   test=0;
	   while(~isempty(champ))
              %premier champ de la liste
	      champc=champ{1};
	      champ(1)=[];
          try
            eval(strcat('test=isstruct(',champc,');'));
          catch
               test = false;
          end
	      if test
	   	   % cas d'une sous structure -> ajout de champs
	   	   eval(strcat('champnew=sort(fieldnames(',champc,'));'));
	   	   for k=1:length(champnew)
	   		   % nom complet pour acceder a la variable
	   		   champnew{k}=strcat(champc,'.',champnew{k});
	   	   end
	   	   % ajout a la liste des champs
	   	   champ=cat(1,champ,champnew);
	      else
	   	   % juste pour les tests
	   	   %disp(champc);
	   	   %fprintf('.');
		   var = eval(champc,[]);
                   ss  = size(var);
                   indnan = find(~isfinite(var(:)));
		   nb     = prod(size(var(:)));
		   if ~isempty(indnan)
		      ll = length(indnan);
		      if ll == nb
		          fprintf('found in %s all NaN or Inf in\n',champc);
		      else
		          fprintf('found in %s  #%d (over % d) NaN or Inf\n',champc,ll,nb);
		      end
		   end
	      end
           end
        end
	disp('---------------')
	disp(' ')
end

if strcmp(cmd,'dataTS')
	% acces au info
	nom   = get(hcb,'tag');
	aide  = get(hcb,'userdata');
	label = get(hcb,'label');
	cb    = get(hcb,'callback');
	producteur = 'pas disponible avec Top 10';
	typep      = 'pas disponible avec Top 10';
	pp    = get(hcb,'parent');
	pp    = get(pp,'parent');
	while ~strcmp(get(pp,'type'),'figure') & (pp ~= 0)
	    stpp = get(pp,'label');
	    if isempty(findstr(stpp,'list'))
	         producteur = stpp;
		 break
	    else
	    	pp    = get(pp,'parent');
	    end
	end
	pp    = get(pp,'parent');
	while ~strcmp(get(pp,'type'),'figure') & (pp ~= 0)
	    stpp = get(pp,'label');
	    if isempty(findstr(stpp,'list'))
	         typep = stpp;
		 break
	    else
	    	pp    = get(pp,'parent');
	    end
	end

	% numreo du choc
	hnum    = findobj_local(hf,'type','uicontrol','tag','numchoc');
	if ishandle(hnum)
		numchoc = str2num(get(hnum,'string'));
	else
		numchoc = [];
	end
	if isempty(numchoc)
		try
		   numchoc = evalin('base','param.from.shot.num');
		   if ishandle(hnum)
		   	set(hnum,'string',sprintf('%g',numchoc));
		   end
		end
	end
	if isempty(numchoc)
		disp('shot number must be provided')
		return

	end


	% lecture de la donnee
	[data,temps,rhofit,certif] = tsbase(fix(numchoc),nom);
	if isempty(temps)
		[data,temps,rhofit,certif] = tsbase(numchoc,nom);
		occ =(numchoc - fix(numchoc)) * 10;
	else
		occ =0;
	end

	if isempty(temps)
		disp('no data')
		return
	elseif isempty(certif)
		certif = rhofit;
		rhofit =[];
	end
	temps  = temps(:,1);
	[cc,cvers,cdate,cheure,cunix,cuniy,cuniz]=tsbase_cert(certif);
	cunix = deblank(cunix);
	cuniy = deblank(cuniy);
	cuniz = deblank(cuniz);

	% creation de la variable dans le workspace
	varout.data    = data;
	varout.temps   = temps;
	varout.rhofit  = rhofit;
	varout.certif  = cc;
	varout.vers    = cvers;
	varout.date    = cdate;
	varout.heure   = cheure;
	varout.unix    = cunix;
	varout.uniy    = cuniy;
	varout.uniz    = cuniz;
	varout.numchoc = fix(numchoc);
	varout.occurrence = occ;

	try
	   zassignin('base',nom,varout);
	end

	% affichage des information
	disp('-------------------------------------------')
	fprintf('Data name : %s@%g\n',nom,numchoc);
	fprintf('Name of the data soure : %s\n',producteur);
	fprintf('source type : %s\n',typep);

	fprintf('Tooltip : %s\n',aide);
	if all(size(data) == 1)
		fprintf('valeur :')
		disp(data)
	else
		fprintf('data size :')
		disp(size(data))
	end
	disp('Information : ')
	fprintf('creation date : the %s at %s\n',cdate,cheure);
	fprintf('version : %d\n',cvers);
	fprintf('certification : ');fprintf('%d ',cc);fprintf('\n');
	fprintf('unities X, Y et Z : %s\t%s\t%s\n',cunix,cuniy,cuniz);
	disp(' ')

	% menu des dernieres donnees
	nbtop =10;
	htop10 = findobj_local(hf,'type','uimenu','tag','top10ts');
	if isempty(htop10)
		htop10 = uimenu(hf,'label','Top10','tag','top10ts');
		deja   = 0;
		hr     = [];
	else
		hr   = get(htop10,'children');
		deja = 0;
		for k= 1:length(hr)
			tt =get(hr(k),'tag');
			if strmatch(tt,nom,'exact')
				deja = 1;
			end
			if get(hr(k),'position') == 1
				hder = hr(k);
			end
		end
	end
	if deja == 0
		if length(hr) < nbtop
			uimenu(htop10,'label',nom,'tag',nom,'userdata',aide,'callback',cb);
		else
			% memorisation des donnees
			for k =1:length(hr)
% 				string{k} = get(hr(k),'label');
% 				tag{k}    = get(hr(k),'tag');
% 				ud{k}     = get(hr(k),'userdata');
% 				cbl{k}    = get(hr(k),'callback');
				pos(k)    = get(hr(k),'position');
			end
			[ppos,ind]  = sort(pos);
			delete(hr(ind(1)))
			uimenu(htop10,'label',nom,'tag',nom,'userdata',aide,'callback',cb);
 		end
	end

	%gestion des unite
	if strcmp(upper(cuniy),'KEV')
	      data =data .* 1e3;
	      cuniy ='eV';
	elseif strcmp(upper(cuniy),'MJ')
	      data =data .* 1e6;
	      cuniy ='J';
	elseif strcmp(upper(cuniy),'MW')
	      data =data .* 1e6;
	      cuniy ='W';
	elseif strcmp(upper(cuniy),'MA')
	      data =data .* 1e6;
	      cuniy ='A';
	elseif strcmp(upper(cuniy),'kPa')
	      data =data .* 1e3;
	      cuniy ='Pa';
	end

	% affichage
	h=findobj_local(hf,'type','axes','tag','axes_des_temps');
	hh=findobj_local(hf,'type','uimenu','tag','efface_temps');
	ck = lower(get(hh,'checked'));
	if isempty(h)|strcmp(ck,'off')
		h=subplot(2,1,1);
		hold off
		set(h,'tag','axes_des_temps');
	else
		axes(h);
	end
	if length(temps) >2
		hp = plot(temps,data,'-');
	else
		hp = plot(temps,data,'o');
	end
	hold on
	xlabel('Time (s)');
	ylabel(sprintf('%s (%s)',nom,cuniy));
	title(sprintf('Shot TS #%d.%d',fix(numchoc),occ))
	co = get(h,'colororder');
	cog=co(1,:);
	set(h,'colororder',co(cat(2,2:size(co,1),1),:));
	set(h,'tag','axes_des_temps');

	if size(data,2) > 1
		if isempty(rhofit)
			rhofit = 1:size(data,2);
		elseif all(abs(rhofit) <1e-38)
			rhofit = 1:size(data,2);
		end

		% les profils
		h=findobj_local(hf,'type','axes','tag','axes_des_profils');
		hh=findobj_local(hf,'type','uimenu','tag','efface_profil');
		ck = lower(get(hh,'checked'));
		if isempty(h)|strcmp(ck,'off')
			h=subplot(2,3,4);
			hold off
			plot(NaN,NaN);
			set(h,'tag','axes_des_profils');
		else
			axes(h);
		end
		if length(temps) > 2
			zplotprof(h,temps,rhofit,data,'color',cog);
		else
			plot(rhofit,data,'linestyle','-','marker','none','color',cog);
		end
		hold on
		xlabel(sprintf('x (%s)',cuniz))
		ylabel(sprintf('%s (%s)',nom,cuniy));
		title('Profiles')
		set(h,'tag','axes_des_profils');

		hm=findobj_local(hf,'type','uimenu','tag','log_image');
		if ~isempty(hm)
			chk = strcmp(get(hm,'checked'),'on');
		else
			chk =0;
		end
		if chk
			datas =log(abs(data))./log(10);
			tnoml = sprintf('log_1_0(%s)',nom);
		else
			datas =data;
			tnoml =nom;
		end

		h = findobj_local(hf,'type','axes','tag','axes_des_images');
		if isempty(h)
		  h  = subplot(2,3,5)
		else
		  axes(h);        
		  if ~verLessThan('matlab','8.5')
		      colorbar(h,'delete');
		  end
		  cla;
		end
		zimagesc(rhofit,temps,datas)
		colormap('default')
        if verLessThan('matlab','8.5')
                colorbar;
        else
                colorbar('peer',h);
        end
		set(gca,'ydir','normal')
		xlabel(sprintf('x (%s)',cuniz))
		ylabel(sprintf('%s (%s)',nom,cuniy));
		title(tnoml)
		set(h,'tag','axes_des_images');
	end

end

if strcmp(cmd,'creeTS')
	% lecture des informations
	fprintf('Access to TS database\n');
	[trait,diag,red]=zgettsinfo(0);
	if isfield(trait,'tprof');
		prof =trait.tprof;
	else
	        prof =[];
	end
	fprintf('Creation of menus , please be patient ...\n');

	% creation du uicontrol num choc
	umem =get(hf,'units');
	set(hf,'units','pixels');
	pos =get(hf,'position');
	set(hf,'units',umem);

	hnum = uicontrol(hf,'style','edit','position',[0,pos(4)-21,75,20],'units','pixels','tag','numchoc');
	set(hnum,'units','normalized');
	try
	    numchoc = evalin('base','param.from.shot.num');
	    set(hnum,'string',sprintf('%g',numchoc));
	end

	% menu puour les diagnostic
	if ~isempty(prof)
		hracine  = uimenu(hf,'label','Tprof','tag','tprof');
		creemenu(hracine,prof,1)
	end

	% menu puour les diagnostic
	hracine  = uimenu(hf,'label','Diagnostics','tag','diagnostic');
	creemenu(hracine,diag)

	% % menu pour les traitements
	hracine  = uimenu(hf,'label','Traitements','tag','traitements');
	creemenu(hracine,trait)

	% menu pour les donnees reduites
	hracine  = uimenu(hf,'label','Donnees reduites','tag','DR');
	creemenu(hracine,red)

	fprintf('\n');
	delete(hcb);
end


if strcmp(cmd,'cherche')
	% dialogue pour selectionner le nom de la varibale
	prompt={'field name  in the structure (ex : ''nequi''), <Empty> -> mots clefs :'};
	def={''};
	dlgTitle='how plot this data ?';
	lineNo=1;
	answer=inputdlg(prompt,dlgTitle,lineNo,def);
	if isempty(answer)
		return
	end

	nomvar   = answer{1};
	if ~isempty(nomvar)
		% liste des mots reservee
		handle_liste = findobj_local(gcf,'type','uimenu');
		reserve = cat(1,get(handle_liste,'tag'),get(handle_liste,'label'));
		liste_nomvar = cat(1,get(handle_liste,'tag'),get(handle_liste,'tag'));
		handle_liste = cat(1,handle_liste,handle_liste);
		clef ={};
		hdl = {};
		for k = 1:length(reserve)
			if isempty(strmatch(reserve{k},clef,'exact'))
				  if ~isempty(strmatch(lower(reserve{k}),lower(nomvar)))
				        cb = get(handle_liste(k),'callback');
				        ua =get(handle_liste(k),'userdata');
				        cp = strrep(cb,')',sprintf(',''%s'')',liste_nomvar{k}));
				        hloc = handle_liste(k);
				        ss3  = sprintf(',findobj_local(gcf,''type'',''uimenu'',''callback'',''%s'',''userdata'',''%s''))',strrep(cb,'''',''''''),ua);
				        hcb = strrep(cb,')',ss3);
				        if ~isempty(findstr(cp,liste_nomvar{k})) && ~isempty(cp)
					    clef{end+1} = cp;
					    hdl{end+1} = hcb;
					end
			          end
			end
		end
		clef =sort(clef);
		if length(clef) >1
		      [s,v] = listdlg('PromptString','Choose a key :',...
				      'SelectionMode','single',...
				      'ListString',clef);
		      if v == 0 | isempty(s)
			      return
		      else
			      nomvar =clef{s};
			      cmd_local = hdl{s};
		      end
		 elseif ~isempty(clef)
		      nomvar =clef{1};
		      cmd_local = hdl{1};	
		 else 
		      nomvar ='';
		      disp('empty result searching key word');
		 end
	end
	if ~isempty(nomvar)
		try
		    evalin('base',cmd_local);
		catch
		    disp(cmd_local)
		end
	end
end

% selection d'un temps
if strcmp(cmd,'pick_t')
    [x,y] = ginput(1);
    if strcmp(get(gca,'tag'),'axes_des_images')
    	x = y;
    elseif strcmp(get(gca,'tag'),'axes_des_temps')
        % rien
    else
   	set(ho,'value',0);
    	return
    end
    hot   = findobj_local(hf,'type','uicontrol','style','edit','tag','temps');
    if ~isempty(ho)
        set(hot,'string',sprintf('%g',x));
	zplotprof('temps');
    end
    set(ho,'value',0);
end


% supperposition
if strcmp(cmd,'supperpose')
      hax  = findobj_local(hf,'type','axes','tag','axes_des_profils');
      if isempty(hax)
         set(ho,'value',0);
         return
      end
      axes(hax);
      dy   = abs(diff(get(hax,'ylim')));
      dx   = abs(diff(get(hax,'xlim')));
      hll  = findobj_local(hax,'type','line','tag','zplotprof');
      if isempty(hll)
          set(ho,'value',0);
          return
      end
      proper = set(hll(1));
      mlist  = proper.Marker(6:(end-1));
      hlm    = findobj_local(hax,'type','line');
      mpris  = get(hlm,'Marker');
      md     =  mlist{1};
      rm     = 0.9;
      for km = 1:length(mlist)
         if isempty(strmatch(mlist{km},mpris))
	    md  = mlist{km};
	    rm  = 0.9 - 0.8 .* km ./length(mlist)
	    break
	 end
      end
      marque = 0;
      for kl = 1:length(hll)
         xd  = get(hll(kl),'xdata');
	 yd  = get(hll(kl),'ydata');
	 ud  = get(hll(kl),'userdata');
	 if ~all(~isfinite(yd)) & ~isempty(ud)
   	    cd  = get(hll(kl),'color');
	    line(xd,yd,'color',cd,'linestyle','-.');
	    lx  = length(xd);
	    d    = abs(xd -rm);
	    indm  = max(find(d == min(d)));
	    ind = [1,indm,lx];
	    line(xd(ind),yd(ind),'color',cd, ...
	          'linestyle','none','marker',md);
            hot   = findobj_local(hf,'type','uicontrol','style','edit','tag','temps');
	    if ~isempty(hot)&(marque == 0)
	          text(xd(indm),yd(indm) + 0.05 .* dy,get(hot,'string'));
	          %line(xd(indm)- 0.05 .* dx,yd(indm) + 0.05 .* dy,'color',[0 0 0], ...
	          %'linestyle','none','marker',md);
		  marque = 1;
	    end
        end
    end
    set(ho,'value',0);
end

if strcmp(cmd,'extrait')
     ha =gca;
     hfnew= figure;
     hnew = gca;
     posnew  = get(hnew,'position');
     unew    = get(hnew,'units');
     delete(hnew);
     hnew = copyobj(ha,hfnew);
     set(hnew,'units',unew,'position',posnew);
     set(hnew,'ButtonDownFcn','');
end

if strcmp(cmd,'data0d')
	% acces au info
	nom   = get(hcb,'tag');
	aide  = get(hcb,'userdata');
	label = get(hcb,'label');
	cb    = get(hcb,'callback');

        % acces a la donnees
        try
               post = evalin('base','post');
        catch
               post =[];
        end
        if isempty(post)
	   disp('no data')
            return
        end
        if ~isfield(post,'zerod')
	   disp('no data')
            return
        end

        data    = getfield(post.zerod,nom);
        temps   = post.zerod.temps;
        numchoc = post.z0dinput.shot;
        machine = post.z0dinput.machine;

	% acces au donnees pour la comparaison
       	data_ref     = [];
        temps_ref    = [];
	if compare
		try
                	post_ref = evalin('base','jeux1.post');
       			data_ref    = getfield(post_ref.zerod,nom);
        		temps_ref   = post_ref.zerod.temps;
		end
	end

	if isempty(temps)
		disp('no data')
		return
	elseif isempty(data)
		disp('no data')
		return
	end

	% acces au donnees pour la comparaison
       	data_ref     = [];
        temps_ref    = [];
	if compare
		try
                	post_ref = evalin('base','jeux1.post');
       			data_ref    = getfield(post_ref.zerod,nom);
        		temps_ref   = post_ref.zerod.temps;
		end
	end

	% affichage des information
	disp('-------------------------------------------')
	fprintf('Name of the data : %s [%s@%g]\n',nom,machine,numchoc);
	fprintf('Source of the data : %s\n','METIS');

	fprintf('Tooltip : %s\n',aide);
	if all(size(data) == 1)
		fprintf('valeur :')
		disp(data)
	else
		fprintf('data size :')
		disp(size(data))
		if isdeployed
		      try
			void = cat(2,temps(:),data(:));
                        fprintf('time\tdata :\n')
			fprintf('%g\t%g\n',void')
		      end
		end
	end
	disp(' ')
	set(hf,'name',sprintf('Help zdataplot @ METIS : %s',aide));


	% affichage
	h=findobj_local(hf,'type','axes','tag','axes_des_temps');
	hh=findobj_local(hf,'type','uimenu','tag','efface_temps');
	ck = lower(get(hh,'checked'));
	if isempty(h)|strcmp(ck,'off')
		h=subplot(2,1,1);
		hold off
		set(h,'tag','axes_des_temps');
	else
		axes(h);
	end
	if length(temps) >2
		hp = plot(temps,real(data),'-');
                if any(imag(data))
		    hold on
		    co = get(h,'colororder');
		    cog=co(1,:);
		    hp = plot(temps,imag(data),'--');
		    set(hp,'color',cog);
		end
	else
		hp = plot(temps,real(data),'o');
                if any(imag(data))
		    hold on
		    co = get(h,'colororder');
		    cog=co(1,:);
		    hp = plot(temps,imag(data),'s');
		    set(hp,'color',cog);
		end
	end
	hold on
	xlabel('Time (s)');
	ylabel(sprintf('%s',formund(nom)));
	title(sprintf('Shot %s %d: %s',machine,fix(numchoc),aide))
	co = get(h,'colororder');
	cog=co(1,:);
	set(h,'colororder',co(cat(2,2:size(co,1),1),:));
	set(h,'tag','axes_des_temps');

	if ~isempty(data_ref)
		if length(temps) >2
			hp = plot(temps_ref,real(data_ref),'-.');
			if any(imag(data_ref))
			    set(hp,'color',cog);
			    hold on
			    hp = plot(temps_ref,imag(data_ref),':');
			end
		else
			hp = plot(temps_ref,real(data_ref),'+');
			if any(imag(data_ref))
			    set(hp,'color',cog);
			    hold on
			    hp = plot(temps_ref,imag(data_ref),'x');
			end
		end
		set(hp,'color',cog);
	end

end

if strcmp(cmd,'exp0d')
	% acces au info
	nom   = get(hcb,'tag');
	aide  = get(hcb,'userdata');
	label = get(hcb,'label');
	cb    = get(hcb,'callback');

        % acces a la donnees
        try
               post = evalin('base','post');
        catch
               post =[];
        end
        if isempty(post)
	   disp('no data')
            return
        end
        if ~isfield(post,'zerod')
	   disp('no data')
            return
        end

        data    = getfield(post.z0dinput.exp0d,nom);
        temps   = post.z0dinput.exp0d.temps;
        numchoc = post.z0dinput.shot;
        machine = post.z0dinput.machine;

	% acces au donnees pour la comparaison
       	data_ref     = [];
        temps_ref    = [];
	if compare
		try
                	post_ref = evalin('base','jeux1.post');
       			data_ref    = getfield(post_ref.z0dinput.exp0d,nom);
        		temps_ref   = post_ref.z0dinput.exp0d.temps;
		end
	end

	if isempty(temps)
		disp('no data')
		return
	elseif isempty(data)
		disp('no data')
		return
	end

	% acces au donnees pour la comparaison
       	data_ref     = [];
        temps_ref    = [];
	if compare
		try
                	post_ref = evalin('base','jeux1.post');
       			data_ref    = getfield(post_ref.z0dinput.exp0d,nom);
        		temps_ref   = post_ref.z0dinput.exp0d.temps;
		end
	end

	% affichage des information
	disp('-------------------------------------------')
	fprintf('Name of the data : %s [%s@%g]\n',nom,machine,numchoc);
	fprintf('Source of the data : %s\n','METIS');

	fprintf('Tooltip : %s\n',aide);
	if all(size(data) == 1)
		fprintf('valeur :')
		disp(data)
	else
		fprintf('data size :')
		disp(size(data))
		if isdeployed
		      try
			void = cat(2,temps(:),data(:));
                        fprintf('time\tdata :\n')
			fprintf('%g\t%g\n',void')
		      end
		end
	end
	disp(' ')
	set(hf,'name',sprintf('Help zdataplot @ METIS : %s',aide));


	% affichage
	h=findobj_local(hf,'type','axes','tag','axes_des_temps');
	hh=findobj_local(hf,'type','uimenu','tag','efface_temps');
	ck = lower(get(hh,'checked'));
	if isempty(h)|strcmp(ck,'off')
		h=subplot(2,1,1);
		hold off
		set(h,'tag','axes_des_temps');
	else
		axes(h);
	end
	if length(temps) >2
		hp = plot(temps,real(data),'-');
                if any(imag(data))
		    hold on
		    co = get(h,'colororder');
		    cog=co(1,:);
		    hp = plot(temps,imag(data),'--');
		    set(hp,'color',cog);
		end
	else
		hp = plot(temps,real(data),'o');
                if any(imag(data))
		    hold on
		    co = get(h,'colororder');
		    cog=co(1,:);
		    hp = plot(temps,imag(data),'s');
		    set(hp,'color',cog);
		end
	end
	hold on
	xlabel('Time (s)');
	ylabel(sprintf('%s',formund(nom)));
	title(sprintf('Shot %s %d: %s',machine,fix(numchoc),aide))
	co = get(h,'colororder');
	cog=co(1,:);
	set(h,'colororder',co(cat(2,2:size(co,1),1),:));
	set(h,'tag','axes_des_temps');

	if ~isempty(data_ref)
		if length(temps) >2
			hp = plot(temps_ref,real(data_ref),'-.');
			if any(imag(data))
			    set(hp,'color',cog);
			    hold on
			    hp = plot(temps,imag(data_ref),':');
			end
		else
			hp = plot(temps_ref,real(data_ref),'+');
			if any(imag(data))
			    set(hp,'color',cog);
			    hold on
			    hp = plot(temps,imag(data_ref),'x');
			end
		end
		set(hp,'color',cog);
	end

end


if strcmp(cmd,'data0d1d')

	% acces au info
	nom   = get(hcb,'tag');
	aide  = get(hcb,'userdata');
	label = get(hcb,'label');
	cb    = get(hcb,'callback');

        % acces a la donnees
        try
               post = evalin('base','post');
        catch
               post =[];
        end
        if isempty(post)
	    disp('no data')
            return
        end
        if ~isfield(post,'profil0d')
	   disp('no data')
            return
        end

        data    = getfield(post.profil0d,nom);
        temps   = post.profil0d.temps;
        rmx     = post.profil0d.rmx;
	xli     = rmx ./ (rmx(:,end) * ones(1,size(rmx,2)));
	if size(xli,2) ~= size(data,2)
		xli =1:size(data,2);
	end
        numchoc = post.z0dinput.shot;
        machine = post.z0dinput.machine;


	if isempty(temps)
		disp('no data')
		return
	elseif isempty(data)
		disp('no data')
		return
	end

	% acces au donnees pour la comparaison
       	data_ref     = [];
        temps_ref    = [];
	rmx_ref      = [];
	xli_ref      = [];
	if compare
		try
                	post_ref = evalin('base','jeux1.post');
       		        data_ref    = getfield(post_ref.profil0d,nom);
        		temps_ref   = post_ref.profil0d.temps;
        		rmx_ref    = post_ref.profil0d.rmx;
			xli_ref     = rmx_ref ./ (rmx_ref(:,end) * ones(1,size(rmx_ref,2)));
			if size(xli_ref,2) ~= size(data_ref,2)
				xli_ref =1:size(data_ref,2);
			end
		end
	end

	% affichage des information
	disp('-------------------------------------------')
	fprintf('Name of the data : %s [%s@%g]\n',nom,machine,numchoc);
	fprintf('Source of the data : %s\n','METIS');

	fprintf('Tooltip : %s\n',aide);
	if all(size(data) == 1)
		fprintf('valeur :')
		disp(data)
	else
		fprintf('data size :')
		disp(size(data))
		if isdeployed
		      try
			void = cat(1,cat(2,NaN,xli(1,:)),cat(2,temps(:),data));
                        fprintf('time\tdata :\n')
                        for k=1:length(temps)
			      fprintf('%s\n',strrep(sprintf('%g\t',void(k,:)'),'NaN','r/a ->'));
                        end
		      catch
			keyboard
		      end
		end
	end
	disp(' ')
	set(hf,'name',sprintf('Help zdataplot  @ METIS: %s',aide));


	% affichage
	h=findobj_local(hf,'type','axes','tag','axes_des_temps');
	hh=findobj_local(hf,'type','uimenu','tag','efface_temps');
	ck = lower(get(hh,'checked'));
	if isempty(h)|strcmp(ck,'off')
		h=subplot(2,1,1);
		hold off
		set(h,'tag','axes_des_temps');
	else
		axes(h);
	end
	if size(data,1) ==length(temps)
		if length(temps) >2
			if size(data,2) > 1
				hp = plot(temps,real(data(:,1)) + imag(data(:,1)),'-');
				hold on
				hp = plot(temps,real(data(:,end)) + imag(data(:,end)),':');
			else
				hp = plot(temps,real(data(:,1)) + imag(data(:,1)),'-');
			end
		else
			hp = plot(temps,real(data) + imag(data),'o');
		end
	else
		disp('Value :')
		data
		return
	end
	hold on
	xlabel('Time (s)');
	ylabel(sprintf('%s',formund(nom)));
	title(sprintf('Shot %s %d: %s',machine,fix(numchoc),aide))
	co = get(h,'colororder');
	cog=co(1,:);
	set(h,'colororder',co(cat(2,2:size(co,1),1),:));
	set(h,'tag','axes_des_temps');


	if ~isempty(data_ref)
		if size(data_ref,1) ==length(temps_ref)
			if length(temps_ref) >2
				if size(data_ref,2) > 1
					hp = plot(temps_ref,real(data_ref(:,1)) + imag(data_ref(:,1)),'-.');
					set(hp,'color',cog);
					hold on
					hp = plot(temps_ref,real(data_ref(:,end)) + imag(data_ref(:,end)),'.');
				else
					hp = plot(temps_ref,real(data_ref(:,1)) + imag(data_ref(:,1)),'-.');
				end
			else
				hp = plot(temps_ref,real(data_ref) + imag(data_ref),'+');
			end
			set(hp,'color',cog);
		else
			disp('Value :')
			data_ref
			return
		end
	end


	if size(data,2) > 1

		% les profils
		h=findobj_local(hf,'type','axes','tag','axes_des_profils');
		hh=findobj_local(hf,'type','uimenu','tag','efface_profil');
		ck = lower(get(hh,'checked'));
		if isempty(h)|strcmp(ck,'off')
			h=subplot(2,3,4);
			hold off
			plot(NaN,NaN);
			set(h,'tag','axes_des_profils');
		else
			axes(h);
		end
		if length(temps) > 2
			zplotprof(h,temps,xli,real(data),'color',cog);
                        if any(imag(data(:)))
			    zplotprof(h,temps,xli,imag(data),'color',cog,'linestyle','--');
                        end
			if ~isempty(data_ref) & (length(temps_ref) > 2)
				zplotprof(h,temps_ref,xli_ref,real(data_ref),'color',cog,'linestyle','-.');
				if any(imag(data_ref(:)))
				    zplotprof(h,temps_ref,xli_ref,imag(data_ref),'color',cog,'linestyle',':');
				end
			end
		else
			plot(xli,real(data),'linestyle','-','marker','none','color',cog);
                        if any(imag(data(:)))
			      plot(xli,imag(data),'linestyle','--','marker','none','color',cog);
                        end
			if ~isempty(data_ref)
				hold on
				plot(xli_ref,real(data_ref),'linestyle','-.','marker','none','color',cog);
				if any(imag(data_ref(:)))
				      plot(xli_ref,imag(data_ref),'linestyle',':','marker','none','color',cog);
				end
			end
		end
		hold on
		xlabel(sprintf('x '))
		ylabel(sprintf('%s',formund(nom)));
		title('Profiles')
		set(h,'tag','axes_des_profils');

		hm=findobj_local(hf,'type','uimenu','tag','log_image');
		if ~isempty(hm)
			chk = strcmp(get(hm,'checked'),'on');
		else
			chk =0;
		end
                data = real(data) + imag(data);
		if chk
			datas =log(abs(data))./log(10);
			tnoml = sprintf('log_1_0(%s)',formund(nom));
		else
			datas =data;
			tnoml =formund(nom);
		end

		h = findobj_local(hf,'type','axes','tag','axes_des_images');
		if isempty(h)
		  h  = subplot(2,3,5)
		else
		  axes(h);        
		  if ~verLessThan('matlab','8.5')
		      colorbar(h,'delete');
		  end
		  cla;
		end
		zimagesc(mean(xli,1),temps,datas)
		colormap('default')
        if verLessThan('matlab','8.5')
                colorbar;
        else
                colorbar('peer',h);
        end
		set(gca,'ydir','normal')
		xlabel(sprintf('x'))
		ylabel(sprintf('time (s)'));
		title(tnoml)
		set(h,'tag','axes_des_images');
	end
end

if strcmp(cmd,'clipboard')
	% ouverture des clipboard
	zgclipp;
	zgclipt;
end

if strcmp(cmd,'jet_sal_list')
    
    % recursivity
    if ishandle(hcb) && ~isempty(strfind(get(hcb,'tag'),'data/pulse'))
       data = sal_list(getappdata(hcb,'JET_DATA_PATH'));       
       if isfield(data,'object') && isfield(data.object,'description')
           fprintf('%s: %s\n',getappdata(hcb,'JET_DATA_PATH'),data.object.description);
       end
       creemenu_jet_sal(hcb,data);
       return
    end
    % try to get shot number
    try
        shot = evalin('base','post.z0dinput.shot;');
    catch
        shot = 87737;
    end
    % menu for one given shot
    prompt={'shot number :','User name','SecurID key'};
	def={sprintf('%d',shot'),'prebut','1234567890'};
	dlgTitle='Access to JET data';

    lineNo=1;
    answer=zinputdlg(prompt,dlgTitle,lineNo,def);
    if isempty(answer)
	    return
    end
    shot  = str2num(answer{1});
    user  = strtrim(answer{2});
    securidkey = strtrim(answer{3});
    setappdata(0,'MDSPLUS_USERNAME',user);
    if (length(securidkey) >= 10) 
        try
            authorisation = sal_auth(user,securidkey);
        catch
            error('SAL identification failed !');
        end
    end
    fprintf('Creation of menus , please be patient ...\n');
    hracine  = uimenu(hf,'label',sprintf('JET@%d',shot),'tag',sprintf('JET@%d',shot));
    setappdata(hracine,'JET_SHOT_DATA',shot);
    setappdata(hracine,'JET_DATA_PATH',sprintf('data/pulse/%d',shot));
    data = sal_list(sprintf('data/pulse/%d',shot));
    creemenu_jet_sal(hracine,data);
	
	fprintf('\n');
	setappdata(0,'ZDATAPLOT_JET_MENU',1);
end

if strcmp(cmd,'imas')
    % lecture des informations
    fprintf('Access to IMAS description \n');
    if isappdata(0,'DESCRIPTION_IMAS')
        description = getappdata(0,'DESCRIPTION_IMAS');
    elseif exist('DESCRIPTION_IMAS_DATA4ZDATAPLOT.mat','file')
        info = load('DESCRIPTION_IMAS_DATA4ZDATAPLOT');
        if ~isempty(strmatch(info.IMAS_VERSION,getenv('IMAS_VERSION'),'exact'))
            description = info.description;
        else
            [output,idss_list,description,mixed] = litidss;
            setappdata(0,'DESCRIPTION_IMAS',description);
            IMAS_VERSION = getenv('IMAS_VERSION');
            save('DESCRIPTION_IMAS_DATA4ZDATAPLOT','description','IMAS_VERSION');
        end
    else
	try
	  [output,idss_list,description,mixed] = litidss;
	  setappdata(0,'DESCRIPTION_IMAS',description);
	  IMAS_VERSION = getenv('IMAS_VERSION');
	  save('DESCRIPTION_IMAS_DATA4ZDATAPLOT','description','IMAS_VERSION');
	catch
	  % for used with imas installed
	  info = load('DATA4IMAS_WITHOUT_INFRASTRUCTURE');
          description  = info.description;
          IMAS_VERSION = NaN;
	end
    end
    fprintf('Creation of menus , please be patient ...\n');

	% creation du uicontrol num choc
	umem =get(hf,'units');
	set(hf,'units','pixels');
	pos =get(hf,'position');
	set(hf,'units',umem);

	hnum = uicontrol(hf,'style','edit','units','pixels','position',[76+225,pos(4)-21,224,20],'tag','numshot','tooltip', ...
	'<shot number>.<run number> or matfile name containing exported imas data (without extension) or directory for HDF5 IMAS files');
	set(hnum,'units','normalized');
		
        hracine  = uimenu(hf,'label','IMAS','tag','IMAS');
        creemenu_imas(hracine,description)
	
	fprintf('\n');
	delete(hcb);
	setappdata(0,'ZDATAPLOT_IMAS_MENU',1);

end

if strcmp(cmd,'jet_sal_plot')
	% acces au info
	nom   = get(hcb,'tag');
	aide  = get(hcb,'userdata');
	label = get(hcb,'label');
	cb    = get(hcb,'callback');
    
    % read data via sal
    ref = getappdata(hcb);
    rep = sal_get(ref.JET_DATA_PATH);
    if isempty(rep)
        disp('no data')
        return
    end
    try
        description = rep.object.description.value;
    catch
        description = '';
    end
    try
       data_units = rep.object.units.value;       
    catch
       data_units = '';
    end
    try
       data = rep.object.data.value.data;      
    catch
       data = [];
    end
    try
       time = rep.object.dimensions.value.x0.value.data.value.data(:);
    catch
       time = [];
    end
    try
       time_units = rep.object.dimensions.value.x0.value.units.value;
    catch
       time_units = 's';
    end
     try
       time_desc = rep.object.dimensions.value.x0.value.description.value;
    catch
       time_desc = 'time';
     end
    try
       x_units = rep.object.dimensions.value.x1.value.units.value;
    catch
       x_units = 'au';
    end
     try
       x_desc = rep.object.dimensions.value.x1.value.description.value;
     catch
       x_desc = 'x';
     end
     try
       x = rep.object.dimensions.value.x1.value.data.value.data(:)';
       if ~isempty(strfind(ref.JET_DATA_PATH,'ppf/signal/tranppf/tra0'))
           if isempty(x) && all(size(data) >1)
               x = linspace(-1,1,size(data,2));
               xloc = x;
               x_desc_loc   = 'r/a'
               x_units_loc  = '';
               
           else
               mem_rrr_tr = sprintf('ROUT_RMAG_RIN_%d',ref.JET_SHOT_DATA); 
               if isappdata(0,mem_rrr_tr)
                   mem = getappdata(0,mem_rrr_tr);
                   rout = mem.rout;
                   rin  = mem.rin;
                   rmag = mem.rmag;
               else
                   rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/trt0/rout',ref.JET_SHOT_DATA ));
                   rout      = rep.object.data.value.data;
                   rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/trt0/rin',ref.JET_SHOT_DATA ));
                   rin       = rep.object.data.value.data;
                   rep       = sal_get(sprintf('data/pulse/%d/ppf/signal/tranppf/trt0/rmag',ref.JET_SHOT_DATA ));
                   rmag      = rep.object.data.value.data;
                   mem.rout = rout;
                   mem.rin  = rin;
                   mem.rmag  = rmag;
                   setappdata(0,mem_rrr_tr,mem);
               end
               vtr       = ones(length(rin),1);
               ver       = ones(1,length(x));
               xli_tr    = (vtr * x - rmag * ver) ./ (rout * ver - rmag * ver) .* ((vtr * x) >= (rmag*ver)) + (vtr * x - rmag * ver) ./ (rmag * ver - rin * ver).* ((vtr * x) < (rmag*ver));
               rmx  = evalin('base','post.profil0d.rmx');
               xli  = evalin('base','post.profil0d.xli');
               temps  = evalin('base','post.profil0d.temps');
               xx        = cat(2,-xli(end:-1:2),xli);
               rmx       = cat(2,-rmx(:,end:-1:2),rmx);
               xloc    = NaN * ones(length(time),length(x));
               for k=1:length(time)
                   indt = find(temps >= time(k),1);
                   if isempty(indt)
                       indt = length(temps);
                   end
                   xloc(k,:) = interp1(xx,rmx(indt,:) ./ max(rmx(indt,:)),xli_tr(k,:),'pchip','extrap');
               end
             
               
%                Raxe = evalin('base','post.profil0d.Raxe');
%                epsi = evalin('base','post.profil0d.epsi');
%                rmx  = evalin('base','post.profil0d.rmx');
%                xli  = evalin('base','post.profil0d.xli');
%                temps  = evalin('base','post.profil0d.temps');
%                Rloc = cat(2,Raxe(:,end:-1:2) .* (1 - epsi(:,end:-1:2)),Raxe .* (1 + epsi));
%                xx   = cat(2,-xli(end:-1:2),xli);
%                rmx  = cat(2,-rmx(:,end:-1:2),rmx);
%                if any(size(x) == 1)
%                    Rtr  = ones(length(time),1) * x(:)';
%                else
%                    Rtr  = x;
%                end
%                xloc    = NaN * ones(length(time),length(x));
%                for k=1:length(time)
%                    indt = find(temps >= time(k),1);
%                    if isempty(indt)
%                        indt = length(temps);
%                    end
%                    xloc(k,:) = interp1(Rloc(indt,:),rmx(indt,:) ./ max(rmx(indt,:)),Rtr(k,:),'pchip','extrap');
%                end
               xloc(abs(xloc)>1.2) = NaN;
               x_desc_loc   = 'r/a'
               x_units_loc  = '';
           end
       else
           xloc = x;
           x_desc_loc   = x_desc;
           x_units_loc  = x_units;
       end
    catch
       x = []; xloc =x;
       x_desc_loc   = x_desc;
       x_units_loc  = x_units;
    end
     try
         certication = rep.object.mask.value.status.value.data;
     catch
         certication = NaN *data;
     end

	% affichage des information
	disp('-------------------------------------------')
	fprintf('Data name : %s\n',ref.JET_DATA_PATH);
	fprintf('Name of the data soure : %s\n','JET');
	fprintf('source type : %s\n','SAL');

	fprintf('Tooltip : %s\n',description);
	if all(size(data) == 1)
		fprintf('valeur :')
		disp(data)
	else
		fprintf('data size :')
		disp(size(data))
	end
	disp('Information : ')
    if ~isempty(certication)
        fprintf('version : %d\n',0);
        if all(certication(:) == certication(1))
            lc =length(certication(:));
            nc = certication(1);
        else
            [lc,nc] = hist(certication(:));
        end
        fprintf('certification : ');
        for k=1:length(nc)
            if lc(k) > 0
                fprintf('%d count: %d, ',nc(k),lc(k));
            end
        end
        fprintf('\n')
    end
	fprintf('unities X, Y et Z : %s\t%s\t%s\n',time_units,x_units,data_units);
    % try to display other field
    noms = fieldnames(rep.object);
    for k=1:length(noms)
       nloc = noms{k};
       try
           if isfield(rep.object.(nloc),'value') && ischar(rep.object.(nloc).value)
              fprintf('%s: %s\n',nloc,rep.object.(nloc).value);
           end
       end
    end
    if isfield(rep.object,'items') && isfield(rep.object.items,'value')
        disp('ITEMS data:')
        % try to display other field
        noms = fieldnames(rep.object.items.value);
        for k=1:length(noms)
            nloc = noms{k};
            try
               fprintf('%s: ',nloc);
               disp(rep.object.items.value.(nloc).value)
            end
        end

    end
	disp(' ')

	% affichage
	h=findobj_local(hf,'type','axes','tag','axes_des_temps');
	hh=findobj_local(hf,'type','uimenu','tag','efface_temps');
	ck = lower(get(hh,'checked'));
	if isempty(h)|strcmp(ck,'off')
		h=subplot(2,1,1);
		hold off
		set(h,'tag','axes_des_temps');
	else
		axes(h);
	end
	if length(time) >2
		hp = plot(time,data,'-');
	else
		hp = plot(time,data,'o');
	end
	hold on
	xlabel(sprintf('%s (%s)',time_desc,time_units));
	ylabel(sprintf('%s (%s)',description,data_units));
	title(sprintf('Shot JET #%d',ref.JET_SHOT_DATA))
	co = get(h,'colororder');
	cog=co(1,:);
	set(h,'colororder',co(cat(2,2:size(co,1),1),:));
	set(h,'tag','axes_des_temps');

	if size(data,2) > 1
		if isempty(x)
			x = 1:size(data,2);
            if isempty(xloc)
                xloc = x;
            end
		end

		% les profils
		h=findobj_local(hf,'type','axes','tag','axes_des_profils');
		hh=findobj_local(hf,'type','uimenu','tag','efface_profil');
		ck = lower(get(hh,'checked'));
		if isempty(h)|strcmp(ck,'off')
			h=subplot(2,3,4);
			hold off
			plot(NaN,NaN);
			set(h,'tag','axes_des_profils');
		else
			axes(h);
		end
		if length(time) > 2
			zplotprof(h,time,xloc,data,'color',cog);
		else
			plot(xloc,data,'linestyle','-','marker','none','color',cog);
		end
		hold on
        xlabel(sprintf('%s (%s)',x_desc_loc,x_units_loc));
        ylabel(sprintf('%s (%s)',description,data_units));
		title('Profiles')
		set(h,'tag','axes_des_profils');

		hm=findobj_local(hf,'type','uimenu','tag','log_image');
		if ~isempty(hm)
			chk = strcmp(get(hm,'checked'),'on');
		else
			chk =0;
		end
		if chk
			datas =log(abs(data))./log(10);
			tnoml = sprintf('log_1_0(%s)',sprintf('%s (%s)',description,data_units));
		else
			datas =data;
			tnoml =sprintf('%s (%s)',description,data_units);
		end

		h = findobj_local(hf,'type','axes','tag','axes_des_images');
		if isempty(h)
		  h  = subplot(2,3,5)
		else
		  axes(h);        
		  if ~verLessThan('matlab','8.5')
		      colorbar(h,'delete');
		  end
		  cla;
		end
		zimagesc(x,time,datas)
		colormap('default')
        if verLessThan('matlab','8.5')
                colorbar;
        else
                colorbar('peer',h);
        end
		set(gca,'ydir','normal')
        xlabel(sprintf('%s (%s)',x_desc,x_units));
        ylabel(sprintf('%s (%s)',time_desc,time_units));
		title(tnoml)
		set(h,'tag','axes_des_images');
	end

end


if strcmp(cmd,'imas_simple')
    % lecture des informations
    fprintf('Access to IMAS description \n');
    if isappdata(0,'DESCRIPTION_IMAS')
        description = getappdata(0,'DESCRIPTION_IMAS');
    elseif exist('DESCRIPTION_IMAS_DATA4ZDATAPLOT.mat','file')
        info = load('DESCRIPTION_IMAS_DATA4ZDATAPLOT');
        if ~isempty(strmatch(info.IMAS_VERSION,getenv('IMAS_VERSION'),'exact'))
            description = info.description;
        else
            [output,idss_list,description,mixed] = litidss;
            setappdata(0,'DESCRIPTION_IMAS',description);
            IMAS_VERSION = getenv('IMAS_VERSION');
            save('DESCRIPTION_IMAS_DATA4ZDATAPLOT','description','IMAS_VERSION');
        end
    else
	try
	  [output,idss_list,description,mixed] = litidss;
	  setappdata(0,'DESCRIPTION_IMAS',description);
	  IMAS_VERSION = getenv('IMAS_VERSION');
	  save('DESCRIPTION_IMAS_DATA4ZDATAPLOT','description','IMAS_VERSION');
	catch
	  % for used with imas installed
	  info = load('DATA4IMAS_WITHOUT_INFRASTRUCTURE');
          description  = info.description;
          IMAS_VERSION = NaN;
	end
    end
    description_simple.core_profiles = description.core_profiles;
    description_simple.core_sources = description.core_sources;
    description_simple.core_transport = description.core_transport;
    description_simple.dataset_description = description.dataset_description;
    description_simple.equilibrium = description.equilibrium;
    description_simple.pulse_schedule = description.pulse_schedule;
    description_simple.summary = description.summary;
    description_simple.transport_solver_numerics = description.transport_solver_numerics;
    description_simple.edge_profiles = description.edge_profiles;
    description_simple.edge_transport = description.edge_transport;
   
    fprintf('Creation of menus , please be patient ...\n');

	% creation du uicontrol num choc
	umem =get(hf,'units');
	set(hf,'units','pixels');
	pos =get(hf,'position');
	set(hf,'units',umem);

	hnum = uicontrol(hf,'style','edit','units','pixels','position',[76,pos(4)-21,480,20],'tag','numshot', ...
	'tooltip','<shot number>.<run number> or matfile name containing exported imas data (without extension) or directory for HDF5 IMAS files');
	set(hnum,'units','normalized');
		
        hracine  = uimenu(hf,'label','IMAS (simplified)','tag','IMAS_SIMPLIFIED');
        creemenu_imas_noerror(hracine,description_simple)
	
	fprintf('\n');
	delete(hcb);
	setappdata(0,'ZDATAPLOT_IMAS_MENU',1);
	
end

if strcmp(cmd,'load_cronos')
    evalin('base','zuiload');
end
if strcmp(cmd,'cronos_data_browser')
    if ~isempty(which('fastzdataplot.fig'))
        hgload('fastzdataplot.fig');
    else
        evalin('base','zdataplot');
    end
end

if strcmp(cmd,'clear_imas')
  noms = fieldnames(getappdata(0));
  for k=1:length(noms)
      if ~isempty(findstr(noms{k},'IMAS_DATA_')) && isappdata(0,noms{k})
            rmappdata(0,noms{k})
      end
  end
end

if strcmp(cmd,'export_imas')
    % acces au info
    nom      = get(hcb,'tag');
    parents  = get(hcb,'userdata');
    aide     = get(hcb,'label');
    cb       = get(hcb,'callback');
    % numreo du choc
    occurrence = '';
    
    hnum    = findobj_local(hf,'type','uicontrol','tag','numshot');
    if ishandle(hnum)
        if isempty(deblank(get(hnum,'string')))
            errordlg('shot number is not defined');
            return
        end
        cache_name = strrep(strrep(get(hnum,'string'),'.','_dot_'),'@','_at_');
        numchoc_str = strrep(get(hnum,'string'),'@0','');
	[numchoc_str,occ] = strtok(numchoc_str,'@');
	if ~isempty(occ)
	    occurrence = deblank(occ(2:end));
	    if strcmp(occurrence,'0')
	      occurrence = '';
	    end
	end
       
        numchoc = str2num(numchoc_str);
        if isempty(numchoc)
            errordlg('shot number is not valid');
            return            
        else
	    shot  = fix(numchoc);
	    [void,run] = strtok(numchoc_str,'.');
	    if isempty(run)
		run = 0;
	    else
		run  = str2num(run(2:end));
	    end
        end
        
        % recover tokamak name and IMAS version 
        if ~isappdata(0,'UAL_TOKAMAK') || ~isappdata(0,'UAL_DATAVERSION') || ~isappdata(0,'UAL_USER')
	  imasdb;
	end
        tokamak = getappdata(0,'UAL_TOKAMAK');
        imas_version = getappdata(0,'UAL_DATAVERSION');
        user = getappdata(0,'UAL_USER');

        [filename, pathname] = uiputfile('*.mat', 'Name of the matfile for exported IMAS data');
	if isequal(filename,0) || isequal(pathname,0)
	    return
	else
        [pseudo_path,filename,ext] = fileparts(filename);
        if ~isempty(occurrence)
            file = fullfile(pathname,sprintf('IMAS_V%s_%s%d_%d_occ_%s_from_%s_%s',imas_version,tokamak,shot,run,occurrence,user,filename));
            [exported_data,idss_list] = litidss('',shot,run,user,tokamak,imas_version,occurrence);
        else
            file = fullfile(pathname,sprintf('IMAS_V%s_%s%d_%d_from_%s_%s',imas_version,tokamak,shot,run,user,filename));
            [exported_data,idss_list] = litidss('',shot,run,user,tokamak,imas_version);
        end
        fprintf('Exporting IMAS data in %s.mat\n',file);
        save(file,'exported_data','imas_version','tokamak','shot','run','filename','occurrence','user');
	end

    else
        errordlg('IMAS menu is not activated');
        return
    end
    

end

if strcmp(cmd,'dataimas')
    
    % acces au info
    nom      = get(hcb,'tag');
    parents  = get(hcb,'userdata');
    aide     = get(hcb,'label');
    cb       = get(hcb,'callback');
    mode_west  = 0;
    other_mode = 0;
    occurrence = '';
    % numreo du choc
    hnum    = findobj_local(hf,'type','uicontrol','tag','numshot');
    if ishandle(hnum)
        if isempty(deblank(get(hnum,'string')))
            errordlg('shot number is not defined');
            return
        end
        cache_name = strrep(strrep(get(hnum,'string'),'.','_dot_'),'@','_at_');
        numchoc_str = strrep(get(hnum,'string'),'@0','');
        [numchoc_str,occ] = strtok(numchoc_str,'@');
        if ~isempty(occ)
            occurrence = deblank(occ(2:end));
            if strcmp(occurrence,'0')
                occurrence = '';
            end
        end
        
        numchoc = str2num(numchoc_str);
        if isempty(numchoc)
            if isappdata(0,sprintf('IMAS_DATA_%s',deblank(get(hnum,'string'))))
                other_mode = 1;
            elseif isdir(strtrim(get(hnum,'string')))
                other_mode = 5;
            elseif exist(sprintf('%s.mat',deblank(get(hnum,'string'))),'file')
                other_mode = 2;
            else
                errordlg('shot number is not defined');
                return
            end
            shot = NaN;
            run  = NaN;
        elseif numchoc < 0
            
            other_mode = 0;
            shot  = fix(numchoc);
            [void,run] = strtok(numchoc_str,'.');
            if isempty(run)
                run = 0;
            else
                run  = str2num(run(2:end));
            end
            cache_name = strrep(strrep(strrep(strrep(get(hnum,'string'),'-','WEST_PUBLIC_'),'.','RUN_'),'.','_dot_'),'@','_at_');
            if isappdata(0,sprintf('IMAS_DATA_%s_%s',parents{1},cache_name))
                other_mode = 4;
            end
            shot = abs(shot);
            mode_west  = 1;
            
        else
            other_mode = 0;
            shot  = fix(numchoc);
            [void,run] = strtok(numchoc_str,'.');
            if isempty(run)
                run = 0;
            else
                run  = str2num(run(2:end));
            end
            if isappdata(0,sprintf('IMAS_DATA_%s_%s',parents{1},cache_name))
                other_mode = 3;
            end
        end
    else
        return
    end
    
    switch other_mode
        case 0
            if mode_west == 1
                if ~isempty(occurrence)
                    fprintf('accessing to IMAS IDS %s for shot %d.%d/%s from WEST data: ',parents{1},shot,run,occurrence);
                else
                    fprintf('accessing to IMAS IDS %s for shot %d.%d from WEST data: ',parents{1},shot,run);
                end
                user           = 'imas_public';
                machine        = 'west';
                version        = getenv('IMAS_VERSION');
                idx  = imas_open_env('ids',shot,run,user, machine,version);
                if ~isempty(occurrence)
                    data = ids_get(idx,sprintf('%s/%s',parents{1},occurrence));
                else
                    data = ids_get(idx,parents{1});
                end
                imas_close(idx);
                if isempty(data)
                    disp('no data')
                    return
                end
                disp('data loaded');
                setappdata(0,sprintf('IMAS_DATA_%s_%s',parents{1},cache_name),data);
                tokamak = 'WEST public';
                
            else
                % acces a la donnees
                if ~isempty(occurrence)
                    fprintf('accessing to IMAS IDS %s for shot %d.%d/%s: ',parents{1},shot,run,occurrence);
                else
                    fprintf('accessing to IMAS IDS %s for shot %d.%d: ',parents{1},shot,run);
                end
                idx  = imas_open('ids',shot,run);
                if ~isempty(occurrence)
                    data = ids_get(idx,sprintf('%s/%s',parents{1},occurrence));
                else
                    data = ids_get(idx,parents{1});
                end
                imas_close(idx);
                if isempty(data)
                    disp('no data')
                    return
                end
                disp('data loaded');
                setappdata(0,sprintf('IMAS_DATA_%s_%s',parents{1},cache_name),data);
                tokamak = getappdata(0,'UAL_TOKAMAK');
            end
        case 1
            fprintf('retrieving IMAS IDS %s from cache: ',parents{1});
            data = getappdata(0,sprintf('IMAS_DATA_%s',cache_name));
            if isfield(data,'run')
                run = data.run;
            end
            if isfield(data,'shot')
                shot = data.shot;
            end
            if isfield(data,'tokamak')
                tokamak = data.tokamak;
            else
                tokamak = 'unknown';
            end
            data = data.exported_data;
            if isfield(data,parents{1})
                data = data.(parents{1});
            else
                disp('no data')
                return
            end
            disp('data loaded');
        case 2
            fprintf('accessing to IMAS IDS %s stored in file %s: ',parents{1},deblank(get(hnum,'string')));
            data = load(deblank(get(hnum,'string')));
            setappdata(0,sprintf('IMAS_DATA_%s',cache_name),data);
            if isfield(data,'run')
                run = data.run;
            end
            if isfield(data,'shot')
                shot = data.shot;
            end
            if isfield(data,'tokamak')
                tokamak = data.tokamak;
            else
                tokamak = 'unknown';
            end
            data = data.exported_data;
            if isfield(data,parents{1})
                data = data.(parents{1});
            else
                disp('no data')
                return
            end
            disp('data loaded');
            
        case 3
            if ~isempty(occurrence)
                fprintf('retrieving IMAS IDS %s from cache for shot %d.%d/%s: ',parents{1},shot,run,occurrence);
            else
                fprintf('retrieving IMAS IDS %s from cache for shot %d.%d: ',parents{1},shot,run);
            end
            data = getappdata(0,sprintf('IMAS_DATA_%s_%s',parents{1},cache_name));
            disp('data loaded');
            tokamak = getappdata(0,'UAL_TOKAMAK');
        case 4
            if ~isempty(occurrence)
                fprintf('retrieving IMAS IDS %s from cache for shot %d.%d/%s from WEST data: ',parents{1},shot,run,occurrence);
            else
                fprintf('retrieving IMAS IDS %s from cache for shot %d.%d from WEST data: ',parents{1},shot,run);
            end
            data = getappdata(0,sprintf('IMAS_DATA_%s_%s',parents{1},cache_name));
            disp('data loaded');
            tokamak = 'WEST public';
        case 5
            path2imasfile = strtrim(get(hnum,'string'));
            if any(path2imasfile) == '@'
                [path2imasfile,occurrence] = strtok(path2imasfile,'@');
                occurrence = occurrence(2:end);
            else
                occurrence = '';
            end
            try
                shot = NaN;
                run  = NaN;
                info = h5info(fullfile(path2imasfile,'master.h5'));
                for k=1:length(info.Attributes)
                   switch lower(info.Attributes(k).Name)
                       case 'shot'
                           shot = info.Attributes(k).Value;
                       case 'run'
                           run = info.Attributes(k).Value;
                   end
                end
            catch
                shot = NaN;
                run  = NaN;
            end
            try
                path2info = split(path2imasfile,'/');
            catch
                path2info = {};
                tobecut   = path2imasfile;
                while ~isempty(tobecut)
                    [path2info{end+1},tobecut] = strtok(tobecut,'/');
                    if ~isempty(tobecut)
                        tobecut = tobecut(2:end);
                    end
                end
            end
                
            % search for tomakak
            for k=1:length(path2info)
                switch strrep(strrep(upper(strtrim(path2info{k})),'/',''),' ','')
                    case {'TS','T-S','TORESUPRA','TORE_SUPRA','TORE-SUPRA'}
                        tokamak = 'Tore_Supra';
                        break;
                    case {'JET','ITER','COMPASS','WEST','EAST','HL-2A','HL-3M','D-IIID','DIIID','D3D','BEAST','ST-40','ST40','CFEDR','CFEFR','SPARK','ARC','JT-60SA'}
                        tokamak = upper(strtrim(path2info{k}));
                        break;
                    otherwise   
                        tokamak = 'unknown';
                end                       
            end
            imasdb(tokamak,getenv('USER'),'3','13');
    end
    
    % parcours de la structure et extraction des donnees
    switch other_mode
        case 5
            % get data
            imas_data_name = parents{2};
            if length(parents) > 2
                for k=3:length(parents)
                    imas_data_name = sprintf('%s.%s',imas_data_name,parents{k});
                end
            end
            data = directread_h5_imasfile(path2imasfile,strtrim(parents{1}),imas_data_name,occurrence);
            leg = sprintf('%s.%s',parents{2},imas_data_name);
            % try to get time
            next_data_name = '';
            time = [];
            time_best = [];           
            for k=2:length(parents)
                if isempty(next_data_name)
                    time_name = 'time';
                    next_data_name = parents{k};
                else
                    time_name = sprintf('%s.%s',next_data_name,'time');
                    next_data_name = sprintf('%s.%s',next_data_name,parents{k});
                end
                try
                   if ~isempty(time)
                       time_best = time;
                   end
                   time = directread_h5_imasfile(path2imasfile,strtrim(parents{1}),time_name,occurrence); 
                   if ~isempty(time) 
                       time = time(:);
                       if any(size(data,1) == length(time))
                          break; 
                       end
                   end
                catch
                   % not that 
                end
            end
            if isempty(time) && ~isempty(time_best)
                time = time_best;
            end
            time = unique(time);
            if ~isempty(time)
                if all(size(data) > 1) && (size(data,1) ~= length(time))
                    if length(size(data)) == 2
                        data = data.';
                    else
                        inds  = find(size(data) == length(time),1);
                        ss    = size(data);
                        index = 1:length(size(data));
                        index(index == inds) = [];
                        if inds == length(ss)
                            index = cat(2,inds,index(end:-1:1));
                       else
                            index = cat(2,inds,index(1:1:end));
                        end
                        data  = reshape(data,ss(index));
                    end
                end
            end
            % renamed
            temps = time;
            
            %  try to search rho_tor_norm
            % try to get time
            next_data_name = '';
            rho_tor_norm = [];
            rho_tor_grid_norm = [];
            rho_tor_norm_best = [];           
            rho_tor_grid_norm_best = [];           
            for k=2:length(parents)
                if isempty(next_data_name)
                    rho_tor_norm_name = 'rho_tor_norm';
                    rho_tor_grid_norm_name = 'grid.rho_tor_norm';
                    next_data_name = parents{k};
                else
                    rho_tor_norm_name = sprintf('%s.%s',next_data_name,'rho_tor_norm');
                    rho_tor_grid_norm_name = sprintf('%s.grid.%s',next_data_name,'rho_tor_norm');
                    next_data_name = sprintf('%s.%s',next_data_name,parents{k});
                end
                try
                   if ~isempty(rho_tor_norm)
                        rho_tor_norm_best = rho_tor_norm;
                   end
                   rho_tor_norm = directread_h5_imasfile(path2imasfile,strtrim(parents{1}),rho_tor_norm_name,occurrence); 
                catch
                   % not that 
                end
                try
                    if ~isempty(rho_tor_grid_norm)
                        rho_tor_grid_norm_best = rho_tor_grid_norm;
                    end
                    rho_tor_grid_norm = directread_h5_imasfile(path2imasfile,strtrim(parents{1}),rho_tor_grid_norm_name,occurrence);
                catch
                    % not that
                end
                
            end
            if isempty(rho_tor_norm) && ~isempty(rho_tor_grid_norm)
                rho_tor_norm = rho_tor_grid_norm;
            end
            if isempty(rho_tor_norm) && ~isempty(rho_tor_norm_best)
                rho_tor_norm = rho_tor_norm_best;
            end
            if isempty(rho_tor_norm) && ~isempty(rho_tor_grid_norm_best)
                rho_tor_norm = rho_tor_grid_norm_best;
            end
            if ~isempty(time)
                if all(size(rho_tor_norm) > 1) && (size(rho_tor_norm,1) ~= length(time))
                    if length(size(rho_tor_norm)) == 2
                        rho_tor_norm = rho_tor_norm.';
                    else
                        inds  = find(size(rho_tor_norm) == length(time),1);
                        ss    = size(rho_tor_norm);
                        index = 1:length(size(rho_tor_norm));
                        index(index == inds) = [];
                        if inds == length(ss)
                            index = cat(2,inds,index(end:-1:1));
                        else
                            index = cat(2,inds,index(1:1:end));
                        end
                        rho_tor_norm  = reshape(rho_tor_norm,ss(index));
                    end
                end
            end
            % try to search identifier
            next_data_name = '';
            identifier = [];
            for k=2:length(parents)
                if isempty(next_data_name)
                    id_name = 'identifier.name';
                    next_data_name = parents{k};
                else
                    id_name = sprintf('%s.%s',next_data_name,'identifier.name');
                    next_data_name = sprintf('%s.%s',next_data_name,parents{k});
                end
                try
                    identifier = directread_h5_imasfile(path2imasfile,strtrim(parents{1}),id_name,occurrence);
                    if ~isempty(identifier)
                        break;
                    end
                catch
                    % not that
                end
            end
            if isempty(identifier)
                %grid_type
                next_data_name = '';
                identifier = [];
                for k=2:length(parents)
                    if isempty(next_data_name)
                        id_name = 'grid_type.name';
                        next_data_name = parents{k};
                    else
                        id_name = sprintf('%s.%s',next_data_name,'grid_type.name');
                        next_data_name = sprintf('%s.%s',next_data_name,parents{k});
                    end
                    try
                        identifier = directread_h5_imasfile(path2imasfile,strtrim(parents{1}),id_name,occurrence);
                        if ~isempty(identifier)
                            break;
                        end
                    catch
                        % not that
                    end
                end
            end
            if ~isempty(identifier)
                % removed repetition
                if (length(size(identifier)) > 1) && (all(size(identifier) >1))
                    identifier = identifier(:,1);
                end                
                if length(identifier) > 1
                    lk = menu('Choose a element in the liste (cell array index): ',identifier);
                else
                    lk = 1;
                end
                ss = size(data);
                indlk = max(find(length(identifier) == ss));
                if ~isempty(indlk)
                    leg = sprintf('%s@%s',leg,identifier{lk});
                    switch length(ss)
                        case 2
                            data = data(:,lk);
                        case 3
                            switch indlk
                                case 2
                                    data = data(:,lk,:);
                                case 3
                                    data = data(:,:,lk);
                                otherwise
                                    warning('not yet implemented');
                            end
                         case 4
                            switch indlk
                                case 2
                                    data = data(:,lk,:,:);
                                case 3
                                    data = data(:,:,lk,:);
                                case 4
                                    data = data(:,:,:,lk);
                                otherwise
                                    warning('not yet implemented');
                            end
                          case 5
                            switch indlk
                                case 2
                                    data = data(:,lk,:,:,:);
                                case 3
                                    data = data(:,:,lk,:,:);
                                case 4
                                    data = data(:,:,:,lk,:);
                                case 5
                                    data = data(:,:,:,:,lk);
                                otherwise
                                    warning('not yet implemented');
                            end
                      otherwise
                            warning('not yet implemented');
                    end
                end
                if ~isempty(rho_tor_norm)
                    ss = size(rho_tor_norm);
                    indlk = max(find(length(identifier) == ss));
                    if ~isempty(indlk)
                        switch length(ss)
                            case 2
                                rho_tor_norm = rho_tor_norm(:,lk);
                            case 3
                                switch indlk
                                    case 2
                                        rho_tor_norm = rho_tor_norm(:,lk,:);
                                    case 3
                                        rho_tor_norm = rho_tor_norm(:,:,lk);
                                    otherwise
                                        warning('not yet implemented');
                                end
                            case 4
                                switch indlk
                                    case 2
                                        rho_tor_norm = rho_tor_norm(:,lk,:,:);
                                    case 3
                                        rho_tor_norm = rho_tor_norm(:,:,lk,:);
                                    case 4
                                        rho_tor_norm = rho_tor_norm(:,:,:,lk);
                                    otherwise
                                        warning('not yet implemented');
                                end
                            case 5
                                switch indlk
                                    case 2
                                        rho_tor_norm = rho_tor_norm(:,lk,:,:,:);
                                    case 3
                                        rho_tor_norm = rho_tor_norm(:,:,lk,:,:);
                                    case 4
                                        rho_tor_norm = rho_tor_norm(:,:,:,lk,:);
                                    case 5
                                        rho_tor_norm = rho_tor_norm(:,:,:,:,lk);
                                    otherwise
                                        warning('not yet implemented');
                                end
                            otherwise
                                warning('not yet implemented');
                        end
                        
                    end
                end
            end
            data = squeeze(data);
            rho_tor_norm = squeeze(rho_tor_norm);
            % search for second identificateur
            if length(size(data)) > 2
                % try to search identifier
                next_data_name = '';
                identifier = [];
                for k=2:length(parents)
                    identifier_mem = identifier;
                    if isempty(next_data_name)
                        id_name = 'identifier.name';
                        next_data_name = parents{k};
                    else
                        id_name = sprintf('%s.%s',next_data_name,'identifier.name');
                        next_data_name = sprintf('%s.%s',next_data_name,parents{k});
                    end
                    try
                        identifier = directread_h5_imasfile(path2imasfile,strtrim(parents{1}),id_name,occurrence);
                    catch
                       identifier = identifier_mem;
                    end
                end
                % removed repetition
                if (length(size(identifier)) > 1) && (all(size(identifier) >1))
                    identifier = identifier(:,1);
                end                
                if ~isempty(identifier) && (any(size(data)) == length(identifier))
                    if length(identifier) > 1
                        lk = menu('Choose a element in the liste (cell array index): ',identifier);
                    else
                        lk = 1;
                    end
                    ss = size(data);
                    indlk = max(find(length(identifier) == ss));
                    if ~isempty(indlk)
                        leg = sprintf('%s+%s',leg,identifier{lk});
                        switch length(ss)
                            case 2
                                data = data(:,lk);
                            case 3
                                switch indlk
                                    case 2
                                        data = data(:,lk,:);
                                    case 3
                                        data = data(:,:,lk);
                                    otherwise
                                        warning('not yet implemented');
                                end
                            case 4
                                switch indlk
                                    case 2
                                        data = data(:,lk,:,:);
                                    case 3
                                        data = data(:,:,lk,:);
                                    case 4
                                        data = data(:,:,:,lk);
                                    otherwise
                                        warning('not yet implemented');
                                end
                            case 5
                                switch indlk
                                    case 2
                                        data = data(:,lk,:,:,:);
                                    case 3
                                        data = data(:,:,lk,:,:);
                                    case 4
                                        data = data(:,:,:,lk,:);
                                    case 5
                                        data = data(:,:,:,:,lk);
                                    otherwise
                                        warning('not yet implemented');
                                end
                            otherwise
                                warning('not yet implemented');
                        end
                    end
                end    
            end
            data = squeeze(data);
             % search for label
            if length(size(data)) > 2
                % try to search label
                next_data_name = '';
                identifier = [];
                for k=2:length(parents)
                    identifier_mem = identifier;
                    if isempty(next_data_name)
                        id_name = 'label';
                        next_data_name = parents{k};
                    else
                        id_name = sprintf('%s.%s',next_data_name,'label');
                        next_data_name = sprintf('%s.%s',next_data_name,parents{k});
                    end
                    try
                       identifier = directread_h5_imasfile(path2imasfile,strtrim(parents{1}),id_name,occurrence);
                    catch
                       identifier = identifier_mem;
                    end
                end
                % removed repetition
                if (length(size(identifier)) > 1) && (all(size(identifier) >1))
                    identifier = identifier(:,1);
                end                
                if ~isempty(identifier)
                    if length(identifier) > 1
                        lk = menu('Choose a element in the liste (cell array index): ',identifier);
                    else
                        lk = 1;
                    end
                    ss = size(data);
                    indlk = max(find(length(identifier) == ss));
                    if ~isempty(indlk)
                        leg = sprintf('%s&%s',leg,identifier{lk});
                        switch length(ss)
                            case 2
                                data = data(:,lk);
                            case 3
                                switch indlk
                                    case 2
                                        data = data(:,lk,:);
                                    case 3
                                        data = data(:,:,lk);
                                    otherwise
                                        warning('not yet implemented');
                                end
                            case 4
                                switch indlk
                                    case 2
                                        data = data(:,lk,:,:);
                                    case 3
                                        data = data(:,:,lk,:);
                                    case 4
                                        data = data(:,:,:,lk);
                                    otherwise
                                        warning('not yet implemented');
                                end
                            case 5
                                switch indlk
                                    case 2
                                        data = data(:,lk,:,:,:);
                                    case 3
                                        data = data(:,:,lk,:,:);
                                    case 4
                                        data = data(:,:,:,lk,:);
                                    case 5
                                        data = data(:,:,:,:,lk);
                                    otherwise
                                        warning('not yet implemented');
                                end
                            otherwise
                                warning('not yet implemented');
                        end
                    end
                end    
            end
            data = squeeze(data);
        otherwise
            
            if isfield(data,'ids_properties') && isfield(data.ids_properties,'homogeneous_time') && (data.ids_properties.homogeneous_time ~= 0)
                temps = data.time;
            else
                temps =[];
            end
            % cas 2D data recover rho_tor_norm
            rho_tor_norm = [];
            leg = parents{1};
            for k=2:length(parents)
                if isempty(temps) && isfield(data,'time')
                    if iscell(data)
                        temps = NaN * ones(length(data),1);
                        for l=1:length(data)
                            temps(l) = data{l}.time;
                        end
                    else
                        temps = data.time;
                    end
                end
                if iscell(data)
                    if length(data) == 1
                        data = data{1};
                        leg  = fullfile(leg,sprintf('%d',1));
                    elseif (length(data) == length(temps)) || isfield(data{1},'time')
                        if length(parents) >= k
                            try
                                if (length(parents) - k) < 8
                                    for l = 1:length(data)
                                        switch length(parents) - k
                                            case 0
                                                if l == 1
                                                    ss = size(data{1}.(parents{k}));
                                                    ss(ss == 1) = [];
                                                    ss = cat(2,length(data),ss);
                                                    res = NaN * ones(ss);
                                                end
                                                inter = data{l}.(parents{k});
                                            case 1
                                                if l == 1
                                                    ss = size(data{1}.(parents{k}).(parents{k+1}));
                                                    ss(ss == 1) = [];
                                                    ss = cat(2,length(data),ss);
                                                    res = NaN * ones(ss);
                                                end
                                                inter = data{l}.(parents{k}).(parents{k+1});
                                            case 2
                                                if l == 1
                                                    ss = size(data{1}.(parents{k}).(parents{k+1}).(parents{k+2}));
                                                    ss(ss == 1) = [];
                                                    ss = cat(2,length(data),ss);
                                                    res = NaN * ones(ss);
                                                end
                                                inter = data{l}.(parents{k}).(parents{k+1}).(parents{k+2});
                                            case 3
                                                if l == 1
                                                    ss = size(data{1}.(parents{k}).(parents{k+1}).(parents{k+2}).(parents{k+3}));
                                                    ss(ss == 1) = [];
                                                    ss = cat(2,length(data),ss);
                                                    res = NaN * ones(ss);
                                                end
                                                inter = data{l}.(parents{k}).(parents{k+1}).(parents{k+2}).(parents{k+3});
                                            case 4
                                                if l == 1
                                                    ss = size(data{1}.(parents{k}).(parents{k+1}).(parents{k+2}).(parents{k+3}).(parents{k+4}));
                                                    ss(ss == 1) = [];
                                                    ss = cat(2,length(data),ss);
                                                    res = NaN * ones(ss);
                                                end
                                                inter = data{l}.(parents{k}).(parents{k+1}).(parents{k+2}).(parents{k+3}).(parents{k+4});
                                            case 5
                                                if l == 1
                                                    ss = size(data{1}.(parents{k}).(parents{k+1}).(parents{k+2}).(parents{k+3}).(parents{k+4}).(parents{k+5}));
                                                    ss(ss == 1) = [];
                                                    ss = cat(2,length(data),ss);
                                                    res = NaN * ones(ss);
                                                end
                                                inter = data{l}.(parents{k}).(parents{k+1}).(parents{k+2}).(parents{k+3}).(parents{k+4}).(parents{k+5});
                                            case 6
                                                if l == 1
                                                    ss = size(data{1}.(parents{k}).(parents{k+1}).(parents{k+2}).(parents{k+3}).(parents{k+4}).(parents{k+6}));
                                                    ss(ss == 1) = [];
                                                    ss = cat(2,length(data),ss);
                                                    res = NaN * ones(ss);
                                                end
                                                inter = data{l}.(parents{k}).(parents{k+1}).(parents{k+2}).(parents{k+3}).(parents{k+4}).(parents{k+6});
                                            case 7
                                                if l == 1
                                                    ss = size(data{1}.(parents{k}).(parents{k+1}).(parents{k+2}).(parents{k+3}).(parents{k+4}).(parents{k+7}));
                                                    ss(ss == 1) = [];
                                                    ss = cat(2,length(data),ss);
                                                    res = NaN * ones(ss);
                                                end
                                                inter = data{l}.(parents{k}).(parents{k+1}).(parents{k+2}).(parents{k+3}).(parents{k+4}).(parents{k+7});
                                            otherwise
                                                error('How do you arrive here ?');
                                        end
                                        switch length(ss)
                                            case {1,2}
                                                res(l,:) = inter;
                                            case 3
                                                res(l,:,:) = inter;
                                            case 4
                                                res(l,:,:,:) = inter;
                                            case 5
                                                res(l,:,:,:,:) = inter;
                                            case 6
                                                res(l,:,:,:,:,:) = inter;
                                            case 7
                                                res(l,:,:,:,:,:,:) = inter;
                                            otherwise
                                                error('Not yet implemented');
                                        end
                                    end
                                    rho_tor_norm = getgridrhotornorm(data,parents);
                                    for ll=k:length(parents)
                                        leg  = fullfile(leg,parents{ll});
                                    end
                                    data = res;
                                    break;
                                    
                                else
                                    disp('to short')
                                    if ~isfield(data{1},parents{k})
                                        fprintf('field "%s" not available in current IMAS dataset\n',parents{k});
                                        return
                                    end
                                    locdat = cell(length(data),1);
                                    for m=1:length(data)
                                        locdat{m}.(parents{k}) = data{m}.(parents{k});
                                        if isfield(data{1},'grid')
                                            locdat{m}.grid = data{m}.grid;
                                        elseif isfield(data{1},'grid_d')
                                            locdat{m}.grid_d = data{m}.grid_d;
                                        elseif isfield(data{1},'grid_d')
                                            locdat{m}.grid_v = data{m}.grid_v;
                                        elseif isfield(data{1},'grid_flux')
                                            locdat{m}.grid_flux = data{m}.grid_flux;
                                        elseif isempty(strcmp(parents{k},'profiles_1d')) && isfield(data{1},'profiles_1d')
                                            locdat{m}.profiles_1d = data{m}.profiles_1d;
                                        end
                                    end
                                    data = swaptime2cronos(locdat);
                                    if isempty(data)
                                        fprintf('field "%s" not available in current IMAS dataset\n',parents{k});
                                        return
                                    end
                                    
                                end
                            catch
                                if ~isfield(data{1},parents{k})
                                    fprintf('field "%s" not available in current IMAS dataset\n',parents{k});
                                    return
                                elseif isstruct(data{1}.(parents{k}))
                                    modesck = 0;
                                else
                                    modesck = 1;
                                    tempvar = data{1}.(parents{k});
                                    liste = cell(length(tempvar),1);
                                    for lk=1:length(tempvar)
                                        % case source
                                        if isfield(tempvar{1},'identifier') && isfield(tempvar{1}.identifier,'name')
                                            if k > 1
                                                if isempty(tempvar{lk}.identifier.name)
                                                    if isfield(tempvar{1}.identifier,'index') && ~isempty(tempvar{1}.identifier.index)
                                                        liste{lk} =  sprintf('%s/%d',parents{k},tempvar{lk}.identifier.index);
                                                    else
                                                        liste{lk} =  sprintf('%s/%d',parents{k},lk);
                                                    end
                                                else
                                                    liste{lk} =  sprintf('%s/%s',parents{k},tempvar{lk}.identifier.name);
                                                end
                                            elseif isempty(tempvar{lk}.identifier.name)
                                                if isfield(tempvar{1}.identifier,'index') && ~isempty(tempvar{lk}.identifier.index)
                                                    liste{lk} =  sprintf('%d',tempvar{lk}.identifier.index);
                                                else
                                                    liste{lk} =  sprintf('%d',lk);
                                                end
                                            else
                                                liste{lk} =  tempvar.identifier.name;
                                            end
                                        elseif isfield(tempvar{1},'identifier') && isfield(tempvar{1}.identifier,'index')
                                            if k > 1
                                                if ~isempty(tempvar{lk}.identifier.index)
                                                    liste{lk} =  sprintf('%s/%d',parents{k},tempvar{lk}.identifier.index);
                                                else
                                                    liste{lk} =  sprintf('%s/%d',parents{k},lk);
                                                end
                                            else
                                                if ~isempty(tempvar{lk}.identifier.index)
                                                    liste{lk} =  sprintf('%d',tempvar{lk}.identifier.index);
                                                else
                                                    liste{lk} =  sprintf('%d',lk);
                                                end
                                            end
                                        else
                                            liste{lk} = sprintf('%s{%d}',parents{k},lk);
                                        end
                                    end
                                    if length(liste) > 1
                                        lk = menu('Choose a element in the liste (cell array index): ',liste);
                                    else
                                        lk = 1;
                                    end
                                end
                                locdat = cell(length(data),1);
                                for m=1:length(data)
                                    if modesck == 1
                                        if length(parents) > k
                                            locdat{m}.(parents{k}){1}.(parents{k+1}) = data{m}.(parents{k}){lk}.(parents{k+1});
                                        else
                                            locdat{m}.(parents{k}){1}.(parents{k+1}) = data{m}.(parents{k}){lk}.(parents{k+1});
                                        end
                                    else
                                        if length(parents) > k
                                            locdat{m}.(parents{k}).(parents{k+1}) = data{m}.(parents{k}).(parents{k+1});
                                        else
                                            locdat{m}.(parents{k}).(parents{k+1}) = data{m}.(parents{k}).(parents{k+1});
                                        end
                                    end
                                    if isfield(data{1},'grid')
                                        locdat{m}.grid = data{m}.grid;
                                    elseif isfield(data{1},'grid_d')
                                        locdat{m}.grid_d = data{m}.grid_d;
                                    elseif isfield(data{1},'grid_d')
                                        locdat{m}.grid_v = data{m}.grid_v;
                                    elseif isfield(data{1},'grid_flux')
                                        locdat{m}.grid_flux = data{m}.grid_flux;
                                    elseif isempty(strcmp(parents{k},'profiles_1d')) && isfield(data{1},'profiles_1d')
                                        locdat{m}.profiles_1d = data{m}.profiles_1d;
                                    end
                                end
                                fprintf('extracting data ');
                                data = swaptime2cronos(locdat);
                                fprintf('-> done\n');
                                if isempty(data)
                                    if length(parents) > k
                                        fprintf('field "%s/%s" not available in current IMAS dataset\n',parents{k},parents{k+1});
                                    else
                                        fprintf('field "%s" not available in current IMAS dataset\n',parents{k});
                                    end
                                    return
                                end
                            end
                        else
                            disp('Why here ? it is not plan !')
                            data = swaptime2cronos(data);
                            if isempty(data)
                                fprintf('field "%s" not available in current IMAS dataset\n',parents{k});
                                return
                            end
                        end
                    else
                        liste = cell(length(data),1);
                        for l=1:length(data)
                            % case source
                            if isfield(data{1},'identifier') && isfield(data{1}.identifier,'name')
                                if k > 1
                                    if isempty(data{l}.identifier.name)
                                        if isfield(data{1}.identifier,'index') && ~isempty(data{l}.identifier.index)
                                            liste{l} =  sprintf('%s/%d',parents{k-1},data{l}.identifier.index);
                                        else
                                            liste{l} =  sprintf('%s/%d',parents{k-1},l);
                                        end
                                    else
                                        liste{l} =  sprintf('%s/%s',parents{k-1},data{l}.identifier.name);
                                    end
                                elseif isempty(data{l}.identifier.name)
                                    if isfield(data{1}.identifier,'index') && ~isempty(data{l}.identifier.index)
                                        liste{l} =  sprintf('%d',data{l}.identifier.index);
                                    else
                                        liste{l} =  sprintf('%d',l);
                                    end
                                else
                                    liste{l} =  data{l}.identifier.name;
                                end
                            elseif isfield(data{1},'identifier') && isfield(data{1}.identifier,'index')
                                if k > 1
                                    if ~isempty(data{l}.identifier.index)
                                        liste{l} =  sprintf('%s/%d',parents{k-1},data{l}.identifier.index);
                                    else
                                        liste{l} =  sprintf('%s/%d',parents{k-1},l);
                                    end
                                else
                                    if ~isempty(data{l}.identifier.index)
                                        liste{l} =  sprintf('%d',data{l}.identifier.index);
                                    else
                                        liste{l} =  sprintf('%d',l);
                                    end
                                end
                            else
                                liste{l} = sprintf('%s{%d}',parents{k},l);
                            end
                        end
                        l = menu('Choose a element in the liste (cell array index): ',liste);
                        data = data{l};
                        leg  = fullfile(leg,sprintf('%d',l));
                    end
                end
                if isempty(rho_tor_norm)
                    rho_tor_norm = getgridrhotornorm(data,parents);
                end
                data = data.(parents{k});
                leg  = fullfile(leg,parents{k});
            end
    end
    
    % affichage des information
    if isempty(data)
        disp('-------------------------------------------')
        fprintf('Name of the data : %s \n',leg);
        if ~isempty(occurrence)
            fprintf('Tokamak : %s  @  shot = %d for run = %d @ occurrence = %s\n',tokamak,shot,run,occurrence);
        else
            fprintf('Tokamak : %s  @  shot = %d for run = %d\n',tokamak,shot,run);
        end
        fprintf('Source of the data : %s\n','IMAS');
        fprintf('Data size : empty\n');
        fprintf('Tooltip : %s\n',aide);
        return
    elseif iscell(data)
        disp('-------------------------------------------')
        fprintf('Name of the data : %s \n',leg);
        if ~isempty(occurrence)
            fprintf('Tokamak : %s  @  shot = %d for run = %d @ occurrence = %s\n',tokamak,shot,run,occurrence);
        else
            fprintf('Tokamak : %s  @  shot = %d for run = %d\n',tokamak,shot,run);
        end
        fprintf('Source of the data : %s\n','IMAS');
        fprintf('Data size : %d\n',size(data));
        fprintf('Data type : timed string\n');
        fprintf('Tooltip : %s\n',aide);
        fprintf('Data:\n')
        disp(data)
        return
    elseif ischar(data)
        disp('-------------------------------------------')
        fprintf('Name of the data : %s \n',leg);
        if ~isempty(occurrence)
            fprintf('Tokamak : %s  @  shot = %d for run = %d @ occurrence = %s\n',tokamak,shot,run,occurrence);
        else
            fprintf('Tokamak : %s  @  shot = %d for run = %d\n',tokamak,shot,run);
        end
        fprintf('Source of the data : %s\n','IMAS');
        fprintf('Data size : %d\n',size(data));
        fprintf('Data type : string\n');
        fprintf('Tooltip : %s\n',aide);
        fprintf('Data:\n')
        disp(data)
        return
    elseif all(size(data) == 1)
        disp('-------------------------------------------')
        fprintf('Name of the data : %s \n',leg);
        if ~isempty(occurrence)
            fprintf('Tokamak : %s  @  shot = %d for run = %d @ occurrence = %s\n',tokamak,shot,run,occurrence);
        else
            fprintf('Tokamak : %s  @  shot = %d for run = %d\n',tokamak,shot,run);
        end
        fprintf('Source of the data : %s\n','IMAS');
        fprintf('Data size : %d\n',size(data));
        fprintf('Data type : scalar\n');
        fprintf('Tooltip : %s\n',aide);
        fprintf('Data:\n')
        disp(data)
        if isempty(temps)
            temps = 1;
        elseif length(temps) > 1
            temps = 1;
        end
    elseif length(size(data)) > 2
        disp('-------------------------------------------')
        fprintf('Name of the data : %s \n',leg);
        if ~isempty(occurrence)
            fprintf('Tokamak : %s  @  shot = %d for run = %d @ occurrence = %s\n',tokamak,shot,run,occurrence);
        else
            fprintf('Tokamak : %s  @  shot = %d for run = %d\n',tokamak,shot,run);
        end
        fprintf('Source of the data : %s\n','IMAS');
        fprintf('Data size : %d\n',size(data));
        fprintf('Data type : multidimentional array\n');
        fprintf('Tooltip : %s\n',aide);
        sd = size(data);
        data = reshape(data,[sd(1),prod(sd(2:end))]);
    elseif all(size(data) > 1)
        if size(data,1) ~= length(temps)
            if size(data,2) == length(temps)
                data = data';
                type_mat = 'timed array';
            else
                temps = 1:size(data,1);
                type_mat = 'non timed array';
            end
        else
            type_mat = 'timed array';
        end
        disp('-------------------------------------------')
        fprintf('Name of the data : %s \n',leg);
        if ~isempty(occurrence)
            fprintf('Tokamak : %s  @  shot = %d for run = %d @ occurrence = %s\n',tokamak,shot,run,occurrence);
        else
            fprintf('Tokamak : %s  @  shot = %d for run = %d\n',tokamak,shot,run);
        end
        fprintf('Source of the data : %s\n','IMAS');
        fprintf('Data size : %d\n',size(data));
        fprintf('Data type : %s\n',type_mat);
        fprintf('Tooltip : %s\n',aide);
    elseif length(data) ~= length(temps)
        data = data(:);
        temps = 1:length(data);
        fprintf('Name of the data : %s \n',leg);
        if ~isempty(occurrence)
            fprintf('Tokamak : %s  @  shot = %d for run = %d @ occurrence = %s\n',tokamak,shot,run,occurrence);
        else
            fprintf('Tokamak : %s  @  shot = %d for run = %d\n',tokamak,shot,run);
        end
        fprintf('Source of the data : %s\n','IMAS');
        fprintf('Data size : %d\n',size(data));
        fprintf('Data type : non timed signal\n');
        fprintf('Tooltip : %s\n',aide);
    else
        fprintf('Name of the data : %s \n',leg);
        if ~isempty(occurrence)
            fprintf('Tokamak : %s  @  shot = %d for run = %d @ occurrence = %s\n',tokamak,shot,run,occurrence);
        else
            fprintf('Tokamak : %s  @  shot = %d for run = %d\n',tokamak,shot,run);
        end
        fprintf('Source of the data : %s\n','IMAS');
        fprintf('Data size : %d\n',size(data));
        fprintf('Data type : timed signal\n');
        fprintf('Tooltip : %s\n',aide);
        if size(data,1) == 1
            data =data';
        end
    end
    set(hf,'name',sprintf('Help zdataplot  @ IMAS: %s',aide));
    
    
    % affichage
    h=findobj_local(hf,'type','axes','tag','axes_des_temps');
    hh=findobj_local(hf,'type','uimenu','tag','efface_temps');
    ck = lower(get(hh,'checked'));
    if isempty(h)|strcmp(ck,'off')
        h=subplot(2,1,1);
        hold off
        set(h,'tag','axes_des_temps');
    else
        axes(h);
    end
    if length(temps) >2
        if size(data,2) > 1
            hp = plot(temps,real(data(:,1)) + imag(data(:,1)),'-');
            hold on
            hp = plot(temps,real(data(:,end)) + imag(data(:,end)),':');
        else
            hp = plot(temps,real(data(:,1)) + imag(data(:,1)),'-');
        end
    else
        hp = plot(temps,real(data) + imag(data),'o');
    end
    hold on
    xlabel('Time (s)');
    ylabel(sprintf('%s',formund_imas(leg)));
    if ~isempty(occurrence)
        title(sprintf('Shot %s %d@%d/%s: %s',tokamak,shot,run,occurrence,aide))
    else
        title(sprintf('Shot %s %d@%d: %s',tokamak,shot,run,aide))
    end
    co = get(h,'colororder');
    cog=co(1,:);
    set(h,'colororder',co(cat(2,2:size(co,1),1),:));
    set(h,'tag','axes_des_temps');
    
    if size(data,2) > 1
        if size(data,2) == size(rho_tor_norm,2)
            xli = rho_tor_norm;
        elseif (size(rho_tor_norm,2) == 1) && size(data,2) == length(rho_tor_norm)
            xli = rho_tor_norm';
        else
            xli = 1:size(data,2);
        end
        
        % les profils
        h=findobj_local(hf,'type','axes','tag','axes_des_profils');
        hh=findobj_local(hf,'type','uimenu','tag','efface_profil');
        ck = lower(get(hh,'checked'));
        if isempty(h)|strcmp(ck,'off')
            h=subplot(2,3,4);
            hold off
            plot(NaN,NaN);
            set(h,'tag','axes_des_profils');
        else
            axes(h);
        end
        if length(temps) > 2
            zplotprof(h,temps,xli,real(data),'color',cog);
            if any(imag(data(:)))
                zplotprof(h,temps,xli,imag(data),'color',cog,'linestyle','--');
            end
            
        else
            plot(xli,real(data),'linestyle','-','marker','none','color',cog);
            if any(imag(data(:)))
                plot(xli,imag(data),'linestyle','--','marker','none','color',cog);
            end
        end
        
        hold on
        xlabel(sprintf('x '))
        ylabel(sprintf('%s',formund_imas(leg)));
        title('Profiles')
        set(h,'tag','axes_des_profils');
        
        hm=findobj_local(hf,'type','uimenu','tag','log_image');
        if ~isempty(hm)
            chk = strcmp(get(hm,'checked'),'on');
        else
            chk =0;
        end
        data = real(data) + imag(data);
        if chk
            datas =log(abs(data))./log(10);
            tnoml = sprintf('log_1_0(%s)',formund_imas(leg));
        else
            datas =data;
            tnoml =formund_imas(leg);
        end
        
        h = findobj_local(hf,'type','axes','tag','axes_des_images');
        if isempty(h)
            h  = subplot(2,3,5)
        else
            axes(h);
            if ~verLessThan('matlab','8.5')
                colorbar(h,'delete');
            end
            cla;
        end
        zimagesc(mean(xli,1),temps,datas)
        colormap('default')
        if verLessThan('matlab','8.5')
            colorbar;
        else
            colorbar('peer',h);
        end
        set(gca,'ydir','normal')
        xlabel(sprintf('x'))
        ylabel(sprintf('time (s)'));
        title(tnoml)
        set(h,'tag','axes_des_images');
    end
end


    
    
% affiche la structure des parametres des modules
function affparametre(info)


% liste des champs la structures
parametre = info.info;
if isempty(parametre)
	disp('Module without any parameter');
	return
end

champ = sort(fieldnames(parametre));
for k=1:length(champ)
	% nom complet pour acceder a la variable
	champ{k}=strcat('parametre.',champ{k});
end

% jusqu'a ce qu'il n'y ait plus de champ
test=0;
while (~isempty(champ))
  	%premier champ de la liste
	champc=champ{1};
	champ(1)=[];
	eval(strcat('test=isstruct(',champc,');'));
	if test
	  	% cas d'une sous structure -> ajout de champs
	   	eval(strcat('champnew=sort(fieldnames(',champc,'));'));
	   	for k=1:length(champnew)
	   		% nom complet pour acceder a la variable
	   		champnew{k}=strcat(champc,'.',champnew{k});
	   	end
	   	% ajout a la liste des champs
	   	if isempty(champ)
	   		champ =champnew;
	   	else
	   		champ=cat(1,champ,champnew);
	   	end
	else
	   	% affichage
		desc = eval(champc,'');
		valeur  = eval(strrep(champc,'parametre.','info.valeur.'),[]);
		borne   = eval(strrep(champc,'parametre.','info.borne.'),[]);
		defaut  = eval(strrep(champc,'parametre.','info.defaut.'),[]);

		fprintf('%s  -> %s :\t',strrep(champc,'parametre.',''),desc);
		disp(valeur);
	end
end


function creemenu(hracine,stinfo,mono)

if nargin <3
  mono =0;
end
cc = sprintf(' ');

if mono == 0
	% liste des producteurs
	noms   = sort(fieldnames(stinfo));
	typep  = get(hracine,'tag');
	% gestion de la longueur du menu
	cliste = 1;
	cprod  = 0;
	nbmax  = 20;
	if length(noms) > nbmax
		hliste = uimenu(hracine,'label',sprintf('list %d',cliste),'tag',sprintf('%s@%d',typep,cliste));
	else
		hliste = hracine;
	end
	% boucle sur les producteur
	for k =1:length(noms)
		nomprod = noms{k};
		stprod  = getfield(stinfo,nomprod);

		% creation du menu producteur
		hprod   = uimenu(hliste,'label',nomprod(nomprod > cc),'tag',nomprod,'userdata',getfield(stprod,'comprod'));
		% commentaire du  procuteur
		h = uimenu(hprod,'label','->','tag','comprod');
		uimenu(h,'label',getfield(stprod,'comprod'));

		% getsion des listes de producteurs
		cprod = cprod + 1;
		if (cprod > nbmax) & (k < (length(noms)-3))
			cliste = cliste + 1;
			hliste = uimenu(hracine,'label',sprintf('list %d',cliste),'tag',sprintf('%s@%d',typep,cliste));
			cprod = 0;
		end


		% menu pour les donnees
		nomdata = sort(fieldnames(stprod));
		% gestion de la longueur du menu
		cliste_data = 1;
		cdata  = 0;
		if length(nomdata) > nbmax
			hld = uimenu(hprod,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
		else
			hld = hprod;
		end

		% bocule sur les donnees
		for l=1:length(nomdata)
			nomc = nomdata{l};
			if ~strcmp(nomc,'comprod')& ~isempty(find(nomc>cc))
				com  = getfield(stprod,nomc);
				h    = uimenu(hld,'label',nomc(nomc >cc),'tag',sprintf('%s@',nomc));
				uimenu(h,'label',com,'tag',nomc,'userdata',com,'callback','zdataplot_metis(''dataTS'')','HandleVisibility','off');
				fprintf('.');
				% gestion des listes
				cdata = cdata + 1;
				if (cdata > nbmax) & (l < (length(nomdata)-3))
					cliste_data = cliste_data + 1;
					hld = uimenu(hprod,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
					cdata  = 0;
				end
			end
		end
	end
else
 	% cas productuer isole
	% menu pour les donnees
	nomdata = sort(fieldnames(stinfo));
	nomprod = get(hracine,'label');
	% gestion de la longueur du menu
	cliste_data = 1;
	cliste = 1;
	cprod  = 0;
	nbmax  = 20;
	cdata  = 0;
	if length(nomdata) > nbmax
		hld = uimenu(hracine,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
	else
		hld = hracine;
	end

	% bocule sur les donnees
	for l=1:length(nomdata)
		nomc = nomdata{l};
		if ~strcmp(nomc,'comprod') & ~isempty(find(nomc>cc))
			com  = getfield(stinfo,nomc);
			h    = uimenu(hld,'label',nomc(nomc>cc),'tag',sprintf('%s@',nomc));
			uimenu(h,'label',com,'tag',nomc,'userdata',com,'callback','zdataplot_metis(''dataTS'')','HandleVisibility','off');
			fprintf('.');
			% gestion des listes
			cdata = cdata + 1;
			if (cdata > nbmax) & (l < (length(nomdata)-3))
				cliste_data = cliste_data + 1;
				hld = uimenu(hracine,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
				cdata  = 0;
			end
		end
	end
end


function cree0d(hracine,stinfo,mono)

if nargin <3
  mono =0;
end
cc = sprintf(' ');

% liste des producteurs
noms   = sort(fieldnames(stinfo));
typep  = get(hracine,'tag');
% gestion de la longueur du menu
cliste = 1;
cprod  = 0;
nbmax  = 20;
% boucle sur les signaux
nomdata = sort(fieldnames(stinfo));
nomprod = get(hracine,'label');
% gestion de la longueur du menu
cliste_data = 1;
cliste = 1;
cprod  = 0;
nbmax  = 20;
cdata  = 0;
if length(nomdata) > nbmax
	hld = uimenu(hracine,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
else
	hld = hracine;
end

% boucle sur les donnees
for l=1:length(nomdata)
	nomc = nomdata{l};
	if ~strcmp(nomc,'comprod') & ~strcmp(nomc,'temps') & ~isempty(find(nomc>cc))
		com  = getfield(stinfo,nomc);
		h    = uimenu(hld,'label',nomc(nomc>cc),'tag',sprintf('%s@',nomc));
		uimenu(h,'label',com,'tag',nomc,'userdata',com,'callback','zdataplot_metis(''data0d'')','HandleVisibility','off');
		fprintf('.');
		% gestion des listes
		cdata = cdata + 1;
		if (cdata > nbmax) & (l < (length(nomdata)-3))
			cliste_data = cliste_data + 1;
			hld = uimenu(hracine,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
			cdata  = 0;
		end
	end
end

function creeexp0d(hracine,stinfo,mono)

if nargin <3
  mono =0;
end
cc = sprintf(' ');

% liste des producteurs
noms   = sort(fieldnames(stinfo));
typep  = get(hracine,'tag');
% gestion de la longueur du menu
cliste = 1;
cprod  = 0;
nbmax  = 20;
% boucle sur les signaux
nomdata = sort(fieldnames(stinfo));
nomprod = get(hracine,'label');
% gestion de la longueur du menu
cliste_data = 1;
cliste = 1;
cprod  = 0;
nbmax  = 20;
cdata  = 0;
if length(nomdata) > nbmax
	hld = uimenu(hracine,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
else
	hld = hracine;
end

% boucle sur les donnees
for l=1:length(nomdata)
	nomc = nomdata{l};
	if ~strcmp(nomc,'comprod') & ~strcmp(nomc,'temps') & ~isempty(find(nomc>cc))
		com  = getfield(stinfo,nomc);
		h    = uimenu(hld,'label',nomc(nomc>cc),'tag',sprintf('%s@',nomc));
		uimenu(h,'label',com,'tag',nomc,'userdata',com,'callback','zdataplot_metis(''exp0d'')','HandleVisibility','off');
		fprintf('.');
		% gestion des listes
		cdata = cdata + 1;
		if (cdata > nbmax) & (l < (length(nomdata)-3))
			cliste_data = cliste_data + 1;
			hld = uimenu(hracine,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
			cdata  = 0;
		end
	end
end


function cree0d1d(hracine,stinfo,mono)

if nargin <3
  mono =0;
end
cc = sprintf(' ');

% liste des producteurs
noms   = sort(fieldnames(stinfo));
typep  = get(hracine,'tag');
% gestion de la longueur du menu
cliste = 1;
cprod  = 0;
nbmax  = 20;
% boucle sur les signaux
nomdata = sort(fieldnames(stinfo));
nomprod = get(hracine,'label');
% gestion de la longueur du menu
cliste_data = 1;
cliste = 1;
cprod  = 0;
nbmax  = 20;
cdata  = 0;
if length(nomdata) > nbmax
	hld = uimenu(hracine,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
else
	hld = hracine;
end

% boucle sur les donnees
for l=1:length(nomdata)
	nomc = nomdata{l};
	if ~strcmp(nomc,'comprod') & ~strcmp(nomc,'xli')  & ~strcmp(nomc,'temps') & ~isempty(find(nomc>cc))
		com  = getfield(stinfo,nomc);
		h    = uimenu(hld,'label',nomc(nomc>cc),'tag',sprintf('%s@',nomc));
		uimenu(h,'label',com,'tag',nomc,'userdata',com,'callback','zdataplot_metis(''data0d1d'')','HandleVisibility','off');
		fprintf('.');
		% gestion des listes
		cdata = cdata + 1;
		if (cdata > nbmax) & (l < (length(nomdata)-3))
			cliste_data = cliste_data + 1;
			hld = uimenu(hracine,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
			cdata  = 0;
		end
	end
end


% fonction pour une image correcte
function zimagesc(x,t,v)

xx = linspace(min(x),max(x),length(x));
tt = linspace(min(t),max(t),length(t))';
try
    vv = interp1(x',v',xx','nearest')';
catch
    disp('non monotonic data')
    return
end
indif = find(diff(t) <=0);
nbmax =100;
if ~isempty(indif)
	while( ~isempty(indif) && (nbmax >0))
		t(indif,:)  = [];
		vv(indif,:) = [];
		indif = find(diff(t) <=0);
		nbmax = nbmax -  1;
	end
end
vv = interp1(t,vv,tt,'nearest');
imagesc(xx,tt,vv);



function aide=litaide(nom,hcb)

aide = get(hcb,'userdata');
if ~isempty(aide)
	if aide(1) ~='?'
		return
	end
end

info  = zinfo;
reste = nom;
while (~isempty(reste))
	[champ,reste] = strtok(reste,'.');
	info = getfield(info,champ);
end
aide   = info;
if ~isempty(aide)
	set(hcb,'userdata',aide);
else
	aide = '???';
end


% traduction des noms francais en anglais
function s = ztr(e)

s ='';
if ~isempty(strfind(e,'idn'))
	s = 'NBI';
end
if ~isempty(strfind(e,'fci'))
	s = 'ICRH';
end
if ~isempty(strfind(e,'fce'))
	s = 'ECRH';
end
if ~isempty(strfind(e,'hyb'))
	s = 'LH';
end
if ~isempty(strfind(e,'fus'))
	s = 'alpha';
end
if ~isempty(strfind(e,'fonction'))
	s = 'module';
end
if ~isempty(strfind(e,'dds'))
	s = 'sawteeth';
end


if ~isempty(s)
	s = sprintf('%s (%s)',e,s);
else
	s = e;
end

function [tor,pol] = decodefce(cons)

angle_mul  = angle(cons);
tor        = fix(abs(angle_mul) .* 1e7) .* 1e-4 - 360;
pol        = (abs(angle_mul) .* 1e7 - fix(abs(angle_mul) .* 1e7)) .* 1e3 - 360;



function nom = formund(nom)

indund = find(nom =='_',1);
other = find(nom =='_');
if ~isempty(indund)
    if length(other) > 1
        other = setxor(other,indund);
        if ~isempty(other)
            nom(other) =[];
        end
    end
    nom = strcat(nom(1:indund-1),'_{',nom(indund+1:end),'}');
end

function nom = formund_imas(nom)

nom = strrep(nom,'_',' ');




function creemenu_imas(hracine,stinfo,parent)

if nargin < 3
	parent ={};
	indice_parent = 1;
else
	indice_parent = length(parent) +1;
end
cc = sprintf(' ');

% liste des producteurs
noms   = sort(fieldnames(stinfo));
typep  = get(hracine,'tag');
% gestion de la longueur du menu
cliste = 1;
cprod  = 0;
nbmax  = 20;
if length(noms) > (nbmax + 3)
	hliste = uimenu(hracine,'label',sprintf('list %d',cliste),'tag',sprintf('%s@%d',typep,cliste),'HandleVisibility','off');
else
	hliste = hracine;
end
% boucle sur les producteur
if indice_parent == 1
	fprintf('IMAS:');
end
for k =1:length(noms)
	nomprod = noms{k};

	% creation du menu producteur
	hprod   = uimenu(hliste,'label',nomprod,'tag',nomprod,'userdata',nomprod,'HandleVisibility','off');
	parent{indice_parent} = nomprod;
	if isstruct(stinfo.(nomprod))
		creemenu_imas(hprod,stinfo.(nomprod),parent)
	else		
		s=char(parent);
		s(s<=32) =[];
		hloc = uimenu(hprod,'label',stinfo.(nomprod),'tag',s(:)','userdata',parent,'callback','zdataplot_metis(''dataimas'')','HandleVisibility','off');
		%get(hloc)
	end
	% getsion des listes de producteurs
	cprod = cprod + 1;
	if (cprod > nbmax) & (k < (length(noms)-3))
		cliste = cliste + 1;
		hliste = uimenu(hracine,'label',sprintf('list %d',cliste),'tag',sprintf('%s@%d',typep,cliste),'HandleVisibility','off');
		cprod = 0;
	end

        if indice_parent == 1
		fprintf('.');
	end
end
if indice_parent == 1
	fprintf('\n');
end



function creemenu_imas_noerror(hracine,stinfo,parent)

if nargin < 3
	parent ={};
	indice_parent = 1;
else
	indice_parent = length(parent) +1;
end
cc = sprintf(' ');

% liste des producteurs
noms   = sort(fieldnames(stinfo)); 
l_noms = 0;
for k=1:length(noms)
  if isempty(strfind(noms{k},'_error_'))
      l_noms = l_noms + 1; 
  end
end
typep  = get(hracine,'tag');
% gestion de la longueur du menu
cliste = 1;
cprod  = 0;
nbmax  = 20;
if l_noms > (nbmax + 3)
	hliste = uimenu(hracine,'label',sprintf('list %d',cliste),'tag',sprintf('%s@%d',typep,cliste),'HandleVisibility','off');
else
	hliste = hracine;
end
% boucle sur les producteur
if indice_parent == 1
	fprintf('IMAS:');
end
for k =1:length(noms)
     nomprod = noms{k};
     if isempty(strfind(nomprod,'_error_'))
	% creation du menu producteur
	hprod   = uimenu(hliste,'label',nomprod,'tag',nomprod,'userdata',nomprod,'HandleVisibility','off');
	parent{indice_parent} = nomprod;
	if isstruct(stinfo.(nomprod))
		creemenu_imas_noerror(hprod,stinfo.(nomprod),parent)
	else
		s=char(parent);
		s(s<=32) =[];
		hloc = uimenu(hprod,'label',stinfo.(nomprod),'tag',s(:)','userdata',parent,'callback','zdataplot_metis(''dataimas'')','HandleVisibility','off');
		%get(hloc)
	end
	% getsion des listes de producteurs
	cprod = cprod + 1;
	if (cprod > nbmax) & (k < (length(noms)-3))
		cliste = cliste + 1;
		hliste = uimenu(hracine,'label',sprintf('list %d',cliste),'tag',sprintf('%s@%d',typep,cliste),'HandleVisibility','off');
		cprod = 0;
	end

        if indice_parent == 1
		fprintf('.');
	end
    end
end
if indice_parent == 1
	fprintf('\n');
end

function rho_tor_norm = getgridrhotornorm(data,parents)

rho_tor_norm =[];

if (isstruct(data) && isfield(data,'grid') && isfield(data.grid,'rho_tor_norm')) || (iscell(data) && isfield(data{1},'grid') && isfield(data{1}.grid,'rho_tor_norm'))
    if iscell(data)
        rho_tor_norm = NaN * ones(length(data),length(data{1}.grid.rho_tor_norm));
        for k=1:length(data)
            rho_tor_norm(k,:) = data{k}.grid.rho_tor_norm;
        end
    else
        rho_tor_norm = data.grid.rho_tor_norm;
    end
elseif (isstruct(data) && isfield(data,'profiles_1d') && isfield(data.profiles_1d,'rho_tor_norm')) || (iscell(data) && isfield(data{1},'profiles_1d') && isfield(data{1}.profiles_1d,'rho_tor_norm'))
    if iscell(data)
        rho_tor_norm = NaN * ones(length(data),length(data{1}.profiles_1d.rho_tor_norm));
        for k=1:length(data)
            rho_tor_norm(k,:) = data{k}.profiles_1d.rho_tor_norm;
        end
    else
        rho_tor_norm = data.profiles_1d.rho_tor_norm;
    end
else
    switch parents{end}
        case 'd'
            if (isstruct(data) && isfield(data,'grid_d') && isfield(data.grid_d,'rho_tor_norm')) || (iscell(data) && isfield(data{1},'grid_d') && isfield(data{1}.grid_d,'rho_tor_norm'))
                if iscell(data)
                    rho_tor_norm = NaN * ones(length(data),length(data{1}.grid_d.rho_tor_norm));
                    for k=1:length(data)
                        if ~isempty(data{k}.grid_d.rho_tor_norm)
                            rho_tor_norm(k,:) = data{k}.grid_d.rho_tor_norm;
                        end
                    end
                else
                    rho_tor_norm = data.grid_d.rho_tor_norm;
                end
            end
        case 'v'
            if (isstruct(data) && isfield(data,'grid_v') && isfield(data.grid_v,'rho_tor_norm')) || (iscell(data) && isfield(data{1},'grid_v') && isfield(data{1}.grid_v,'rho_tor_norm'))
                if iscell(data)
                    rho_tor_norm = NaN * ones(length(data),length(data{1}.grid_v.rho_tor_norm));
                    for k=1:length(data)
                        if  ~isempty(data{k}.grid_v.rho_tor_norm)
                            rho_tor_norm(k,:) = data{k}.grid_v.rho_tor_norm;
                        end
                    end
                else
                    rho_tor_norm = data.grid_v.rho_tor_norm;
                end
            end
        case 'flux'
            if (isstruct(data) && isfield(data,'grid_flux') && isfield(data.grid_flux,'rho_tor_norm')) || (iscell(data) && isfield(data{1},'grid_flux') && isfield(data{1}.grid_flux,'rho_tor_norm'))
                if iscell(data)
                    rho_tor_norm = NaN * ones(length(data),length(data{1}.grid_flux.rho_tor_norm));
                    for k=1:length(data)
                        if  ~isempty(data{k}.grid_flux.rho_tor_norm)
                            rho_tor_norm(k,:) = data{k}.grid_flux.rho_tor_norm;
                        end
                    end
                else
                    rho_tor_norm = data.grid_flux.rho_tor_norm;
                end
            end
    end
    if isempty(rho_tor_norm)
        if (isstruct(data) && isfield(data,'grid_d') && isfield(data.grid_d,'rho_tor_norm')) || (iscell(data) && isfield(data{1},'grid_d') && isfield(data{1}.grid_d,'rho_tor_norm'))
             if iscell(data)
                rho_tor_norm = NaN * ones(length(data),length(data{1}.grid_d.rho_tor_norm));
                for k=1:length(data)
                        if ~isempty(data{k}.grid_d.rho_tor_norm)
                            rho_tor_norm(k,:) = data{k}.grid_d.rho_tor_norm;
                        end
                end
            else
                rho_tor_norm = data.grid_d.rho_tor_norm;
            end
            
        elseif (isstruct(data) && isfield(data,'grid_v') && isfield(data.grid_v,'rho_tor_norm')) || (iscell(data) && isfield(data{1},'grid_v') && isfield(data{1}.grid_v,'rho_tor_norm'))
            if iscell(data)
                rho_tor_norm = NaN * ones(length(data),length(data{1}.grid_v.rho_tor_norm));
                for k=1:length(data)
                        if ~isempty(data{k}.grid_v.rho_tor_norm)
                                rho_tor_norm(k,:) = data{k}.grid_v.rho_tor_norm;
                        end
                end
            else
                rho_tor_norm = data.grid_v.rho_tor_norm;
            end
        elseif (isstruct(data) && isfield(data,'grid_flux') && isfield(data.grid_flux,'rho_tor_norm')) || (iscell(data) && isfield(data{1},'grid_flux') && isfield(data{1}.grid_flux,'rho_tor_norm'))
            if iscell(data)
                rho_tor_norm = NaN * ones(length(data),length(data{1}.grid_flux.rho_tor_norm));
                for k=1:length(data)
                        if ~isempty(data{k}.grid_flux.rho_tor_norm)
                            rho_tor_norm(k,:) = data{k}.grid_flux.rho_tor_norm;
                        end
                end
            else
                rho_tor_norm = data.grid_flux.rho_tor_norm;
            end
        end
    end

end


function creemenu_jet_sal(hracine,data)

%retreive data
shot      = getappdata(hracine,'JET_SHOT_DATA');
data_path = getappdata(hracine,'JET_DATA_PATH');

% create new level of menu
if isfield(data.object,'children') && ~isempty(data.object.children)
    if isfield(data.object.children,'branches') && ~isempty(data.object.children.branches)
        
        nomprod = get(hracine,'label');
        cliste_data = 1;
        cliste = 1;
        cprod  = 0;
        nbmax  = 20;
        if (length(data.object.children.branches)/nbmax) > 40
            nbmax = ceil(sqrt(length(data.object.children.branches)));
        end
        cdata  = 0;
        if length(data.object.children.branches) > nbmax
            hld = uimenu(hracine,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
        else
            hld = hracine;
        end
        
        for k=1:length(data.object.children.branches)
            fpath = sprintf('%s/%s',data_path,data.object.children.branches{k});
            hloc = uimenu(hld,'label',data.object.children.branches{k},'tag',fpath, ...
                'MenuSelectedFcn','zdataplot_metis(''jet_sal_list'')', ...
                'Callback','zdataplot_metis(''jet_sal_list'')', ...
                'ButtonDownFcn','zdataplot_metis(''jet_sal_list'')', ...
                'HandleVisibility','off','BusyAction','cancel');
            setappdata(hloc,'PARENT',hracine);
            setappdata(hloc,'JET_SHOT_DATA',shot);
            setappdata(hloc,'JET_DATA_PATH',fpath);
            
            % gestion des listes
            cdata = cdata + 1;
            if (cdata > nbmax) && (k < (length(data.object.children.branches)-3))
                cliste_data = cliste_data + 1;
                hld = uimenu(hracine,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
                cdata  = 0;
            end
            
        end
    end
    if isfield(data.object.children,'leaves') && ~isempty(data.object.children.leaves)
                
        nomprod = get(hracine,'label');
        cliste_data = 1;
        cliste = 1;
        cprod  = 0;
        nbmax  = 20;
        if (length(data.object.children.leaves)/nbmax) > 40
            nbmax = ceil(sqrt(length(data.object.children.leaves)));
        end

        cdata  = 0;
        if length(data.object.children.leaves) > nbmax
            hld = uimenu(hracine,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
        else
            hld = hracine;
        end
        
        for k=1:length(data.object.children.leaves)
            fpath = sprintf('%s/%s',data_path,data.object.children.leaves(k).name);
            hloc = uimenu(hld,'label',data.object.children.leaves(k).name,'tag',fpath, ...
                'MenuSelectedFcn','zdataplot_metis(''jet_sal_plot'')', ...
                'Callback','zdataplot_metis(''jet_sal_plot'')', ...
                'ButtonDownFcn','zdataplot_metis(''jet_sal_plot'')', ...
                'HandleVisibility','off','BusyAction','cancel');
            setappdata(hloc,'PARENT',hracine);
            setappdata(hloc,'JET_SHOT_DATA',shot);
            setappdata(hloc,'JET_DATA_PATH',fpath);
            
%             rep = sal_list(fpath);
%             if isfield(rep,'object') && isfield(rep.object,'description')
%                 fprintf('%s: %s\n',fpath,rep.object.description);
%             end
%             
            % gestion des listes
            cdata = cdata + 1;
            if (cdata > nbmax) && (k < (length(data.object.children.leaves)-3))
                cliste_data = cliste_data + 1;
                hld = uimenu(hracine,'label',sprintf('list %d',cliste_data),'tag',sprintf('%s@%d',nomprod,cliste_data));
                cdata  = 0;
            end
            
        end
    end
    set(hracine,'MenuSelectedFcn','','callback','','ButtonDownFcn','');

else
   set(hracine,'MenuSelectedFcn','zdataplot_metis(''jet_sal_plot'')', ...
                'Callback','zdataplot_metis(''jet_sal_plot'')', ...
                'ButtonDownFcn','zdataplot_metis(''jet_sal_plot'')');
            
            
   zdataplot_metis('jet_sal_plot');
            
end
drawnow


