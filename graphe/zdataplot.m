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
function zdataplot(cmd,hcb)


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
   		if strmatch(status_debug.cond,'error')

	dbclear error
   		end
   	end
   end

	% ouverture des clipboard
	%zgclipp;
	%zgclipt;

	% creation de la figure
	hf=figure('tag','zdataplot','toolbar','figure','name','zdataplot');
	set(hf,'defaultaxesfontsize',12,'units','normalized','position',[0.1 0.1 0.7 0.7], ...
	    'color',[1 1 1]);
	% liste des donnees
	try
	    data  = evalin('base','data');
	    param = evalin('base','param');
	catch
	    data=[];
	    param=[];
	end
	if isempty(data) | isempty(param)
		[cr,data,param]=zinit('',linspace(0,10,11)',21,'',3,7,2,2,2,2,2,12);
	elseif ~isstruct(data) | ~isstruct(param)
		[cr,data,param]=zinit('',linspace(0,10,11)',21,'',3,7,2,2,2,2,2,12);
	end
	% suppression des donnees du profiler
	param.profile.data=[];
	% information
	try
	   info=zinfo;
	catch
	   disp('tooltip not available');
	   info =[];
	end



	% creation du menu pour param
	hmparam=uimenu(hf,'label','Param','tag','param');
	% % 1 - la structure param
        % liste des champs la structures
	champ = sort(fieldnames(param));
	for k=1:length(champ)
		% nom complet pour acceder a la variable
		champ{k}=strcat('param.',champ{k});
	end


	% jusqu'a ce qu'il n'y ait plus de champ
	test=0;
	while (~isempty(champ))
	   %premier champ de la liste
	   champc=champ{1};
	   champ(1)=[];
      try
	        eval(strcat('test=isstruct(',champc,');'));
      catch
           test = 0;
      end
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
		champ = sort(champ);
	   	% juste pour les tests
		%if ~isempty(findstr(champc,'cons'))
	   	%    disp(champc);
		%    keyboard
		%end
	   else
	   	% juste pour les tests
		%if ~isempty(findstr(champc,'memoire'))
	   	%    disp(champc);
		%    keyboard
		%end

	   	fprintf('.');

	   	% creation du menu
	   	champf = strcat(champc,'.');
	   	ind = find(champf == '.');
	   	hmem = hf;
		hfle=[];
	   	for k = 1:length(ind)
	   		nom = champf(1:(ind(k)-1));
	   		h   = findobj(hf,'type','uimenu','tag',nom);
	   		if isempty(h)
	   			indl = find(nom =='.');
	   			if isempty(indl)
	   				noml = nom;
	   			else
	   				noml = nom((max(indl) +1):end);
	   			end
	   			hmem = uimenu(hmem,'tag',nom,'Label',noml);
				% modification du 08/06/2001 -> affichage strucure complete
				try
				   if eval(strcat('isstruct(',nom,')'))
				       hfle=uimenu(hmem,'Label','->','tag',strcat('@',nom));
				       set(hfle,'callback','zdataplot(''struct'')','userdata','');

				   end
				catch
					%keyboard
				end
	   		else
	   			% test sur le nombres d'enfant
	   			ch=get(h,'children');
	   			n=1;
	   			while (length(ch) >= 25) & (k>1)
	   				n = n+1;
	   				h = findobj(hf,'type','uimenu','tag',strcat(nom,int2str(n)));
	   				if isempty(h)
	   					indl = find(nom =='.');
	   					if isempty(indl)
	   						noml = nom;
	   					else
	   						noml = nom((max(indl) +1):end);
	   					end
						h = uimenu(hmem,'tag',strcat(nom,int2str(n)),'Label',noml);
	   					ch=[];
	   				else
	   					ch=get(h,'children');
	   				end
	   			end
	   			hmem = h;
	   		end
	   	end
	   	txt=get(hmem,'tag');
	   	try
	   	      hlp = getfield(info,txt);
	   	catch
	   	      hlp = '???';
	        end
	   	set(hmem,'callback','zdataplot(''param'')','userdata',hlp);
		if ~isempty(hfle)
	   	   set(hfle,'userdata',hlp);
		end

	   end

	end

	% 2 - la structure data
	% creation du menu pour param
	hmparam=uimenu(hf,'label','Data','tag','data');
	% liste des champs la structures
	champ = sort(fieldnames(data));
	for k=1:length(champ)
		% nom complet pour acceder a la variable
		champ{k}=strcat('data.',champ{k});
	end

	% jusqu'a ce qu'il n'y ait plus de champ
	test=0;
	while(~isempty(champ))
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
		champ = sort(champ);
	   else
	   	% juste pour les tests
	   	%disp(champc);
	   	fprintf('.');

	   	% creation du menu
	   	champf = strcat(champc,'.');
	   	ind = find(champf == '.');
	   	hmem = hf;
		hfle =[];
	   	for k = 1:length(ind)
	   		nom = champf(1:(ind(k)-1));
	   		h   = findobj(hf,'type','uimenu','tag',nom);
	   		if isempty(h)
	   			indl = find(nom =='.');
	   			if isempty(indl)
	   				noml = nom;
	   			else
	   				noml = nom((max(indl) +1):end);
	   			end

	   			hmem = uimenu(hmem,'tag',nom,'Label',noml);
				% modification du 08/06/2001 -> affichage strucure complete
				try
				   if eval(strcat('isstruct(',nom,')'))
				       hfle=uimenu(hmem,'Label','->','tag',strcat('@',nom));
	   	   		       set(hfle,'callback','zdataplot(''struct'')','userdata','');
				   end
				catch
					%keyboard
				end
	   		else

	   			% test sur le nombres d'enfant
	   			ch=get(h,'children');
	   			n=1;
	   			while (length(ch) >= 25) & (k > 1)
	   				n = n+1;
	   				h = findobj(hf,'type','uimenu','tag',strcat(nom,int2str(n)));
	   				if isempty(h)
	   					indl = find(nom =='.');
	   					if isempty(indl)
	   						noml = nom;
	   					else
	   						noml = nom((max(indl) +1):end);
	   					end
	   					h = uimenu(hmem,'tag',strcat(nom,int2str(n)),'Label',noml);
	   					ch=[];
	   				else
	   					ch=get(h,'children');
	   				end
	   			end
	   			hmem = h;
	   		end

	   	end
	   	txt=get(hmem,'tag');
	   	try
	   	      hlp = getfield(info,txt);
	   	catch
	   	      hlp = '???';
	   	end
	   	set(hmem,'callback','zdataplot(''data'')','userdata',hlp);
		if ~isempty(hfle)
	   	   set(hfle,'userdata',hlp);
		end
       end
    end
    fprintf('\n');


    % substituion des labels
    hso = findobj(hf,'type','uimenu');
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
		 'callback','zdataplot(''pick_t'');', ...
		 'tooltip','Select time slice with the mouse', ...
		 'tag','pick_t');
    % suite des axes
    h=subplot(4,1,2);
    set(h,'tag','axes_des_phases','xticklabel',[]);
    %xlabel('temps (s)')
    %title('angle(f(t))')
    h=subplot(2,3,4);
    set(h,'tag','axes_des_profils');
    xlabel('r','Fontname','Symbol')
    title('Profiles')
    % le bouton pour selectionner le temps
    poscc = get(h,'position');
    uicontrol(hf,'style','checkbox', ...
                 'units','normalized', ...
                 'position',[0,poscc(2)+0.5*poscc(4),0.02,0.02], ...
		 'callback','zdataplot(''supperpose'');', ...
		 'tooltip','Hold on time slice', ...
		 'tag','supperpose');
    % suite des axes
    h=subplot(2,3,5);
    set(h,'tag','axes_des_images');
    xlabel('r','Fontname','Symbol')
    title('2D Profiles')
    h=subplot(2,3,6);
    set(h,'tag','axes_des_surfaces');
    xlabel('R (m)')
    ylabel('Z (m)')
    title('Equilibrium')

    % le menu des commandes
    h=uimenu(hf,'label','Commands','tag','commande');
    uimenu(h,'label','Plot also reference', ...
              'Callback','zdataplot(''plot_ref'')', ...
              'tag','plot_ref');
    uimenu(h,'label','Hold on time traces', ...
              'Callback','zdataplot(''efface_temps'')', ...
              'tag','efface_temps');
    uimenu(h,'label','Hold on profiles', ...
              'Callback','zdataplot(''efface_profils'')', ...
              'tag','efface_profil');
    uimenu(h,'label','LogY time traces', ...
              'Callback','zdataplot(''log_temps'')', ...
              'tag','log_temps');
    uimenu(h,'label','LogY profiles', ...
              'Callback','zdataplot(''log_profils'')', ...
              'tag','log_profil');
    uimenu(h,'label','Log10 2D image', ...
              'Callback','zdataplot(''log_image'')', ...
              'tag','log_image');

    uimenu(h,'label','Plot equilibrium', ...
              'Callback','zdataplot(''plot_equi'')', ...
              'tag','plot_equi');

    uimenu(h,'label','Plot profiler', ...
              'Callback','zdataplot(''plot_profile'')', ...
              'tag','plot_profile');
    uimenu(h,'label','Report profiler', ...
              'Callback','zdataplot(''rapport_profile'')', ...
              'tag','rapport_profile');

    uimenu(h,'label','Look for NaN or Inf', ...
              'Callback','zdataplot(''NaNInf'')', ...
              'tag','NaNInf');

    uimenu(h,'label','Extract axes', ...
              'Callback','zdataplot(''extrait'')', ...
              'tag','extrait');

    uimenu(h,'label','Look for variable', ...
              'Callback','zdataplot(''cherche'')', ...
              'tag','cherche');

    uimenu(h,'label','Access to TS data', ...
              'Callback','zdataplot(''creeTS'')', ...
              'tag','creeTS');

    uimenu(h,'label','Access to METIS data', ...
              'Callback','zdataplot(''cree0D'')', ...
              'tag','cree0D');

    uimenu(h,'label','Create fast visualisation window', ...
              'Callback','zdataplot(''fast'')', ...
              'tag','make pre-built , non dynamic, zdataplot');

    uimenu(h,'label','Open Clipboard', ...
              'Callback','zdataplot(''clipboard'')', ...
              'tag','open clipboards for axes');

   % restauration du debgug sur error
   %status_debug = dbstatus;
   if ~isempty(status_debug)
   	if ~isempty({status_debug.cond})
   		if strmatch(status_debug.cond,'error')
   			dbstop if error
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

compare = strcmp(get(findobj(hf,'type','uimenu','tag','plot_ref'),'checked'),'on');

if isempty(ho)&isempty(hcb)
	return
end


if strcmp(cmd,'param')
	nom  = get(hcb,'tag');
	aide = litaide(nom,hcb);
	nomc = strrep(nom,'param.','');
	infofun =[];
	try
	   nomv = nom;
	   data=evalin('base',nom);
	   if strmatch('param.fonction.',nom)
	      infofun = feval(data);
	   end
	catch
	   try
	     nomv = nomc;
	     data=evalin('caller',nomc);
	     if strmatch('param.fonction.',nomc)
	        infofun = feval(data);
	     end
	  catch
	      nomv = nom;
	      data =[];
	   end
	end
	if isempty(data)
		data='<Empty>';
	end

	disp('-------------------------------------------')
	fprintf('Parameter name : %s\n',nom);
	fprintf('Variable  name : %s\n',nomv);
	if ~isempty(infofun)
	    aide =strrep(aide,'{','\n{');
	end
	fprintf('Help : %s\n',aide);
	set(hf,'name',sprintf('Help zdataplot : %s',aide));
	disp('data :')
	disp(' ')
	disp(data)
	disp(' ')
	indnan=find(~isfinite(data));
	if ~isempty(indnan)
	    fprintf('Contains NaN or Inf (#%d)\n\n',length(indnan));
	end
	indimag=find(imag(data));
	if ~isempty(indimag)
	    fprintf('Contains complex numbers (#%d)\n\n',length(indimag));
	end
	disp(' ')
	if ~isempty(infofun)
	   disp('function description :')
	   ind = find(infofun.description =='{')-1;
	   if isempty(ind)
	      ind = length(infofun.description);
	   end
	   disp(infofun.description(1:ind))
	   disp('parameters :')
	   affparametre(infofun)
	end

end

if strcmp(cmd,'data')
	nom  = get(hcb,'tag');
	aide=litaide(nom,hcb);
	nomk = strrep(nom,'data.','datak.');
	nomkp1 = strrep(nom,'data.','datakp1.');

	% gestion des consignes
	modeimag = 0;
	modeylabel  = 0;
	modefce   =0;
	if ~isempty(findstr(nom,'cons.idn'))
	   	modeimag =1;
		modeylabel =0;
	end
	if ~isempty(findstr(nom,'cons.fce'))
		modefce =1;
	end
	if ~isempty(findstr(nom,'mhd.gamma'))
	   	modeimag =1;
		modeylabel =1;
	end
	if ~isempty(findstr(nom,'pfcur'))
	   	modeimag =1;
		modeylabel =1;
	end

	% cas des donnees cronos
	% info choc
	desctitre ='';
	try
	   machine = evalin('base','param.from.machine');
	   choc = evalin('base','param.from.shot.num');
	   tmin = evalin('base','param.gene.tdeb');
	   tmax = evalin('base','param.gene.tfin');
	   desctitre = sprintf('Shot %s #%d, from %g to %g s',machine,fix(choc),tmin,tmax);
	end
	  %
	  ok   = 1;
	  mode = 0;
	  jeux1 =[];
	  try
	     data  = double(evalin('base',nom));
	     temps = evalin('base','data.gene.temps');
	     x     = evalin('base','param.gene.x');
	     kind  = evalin('base','param.gene.k');
	     rho   = evalin('base','data.equi.rhoRZ');
	   catch
	     ok=0;
	   end
	   if ok == 0
	   	ok   = 1;
	   	mode =1;
	   	try
	   	   datak     = double(evalin('caller',nomk));
	   	   datakp1   = double(evalin('caller',nomkp1));
		      rhok      = evalin('caller','datak.equi.rhoRZ');
		      rhokp1    = evalin('caller','datakp1.equi.rhoRZ');
		      if size(datak,2)>1
		   	   x    = linspace(0,1,size(datak,2));
		      else
		   	   x    = linspace(0,1,21);
		      end
		      tempsk    = evalin('caller','datak.gene.temps');
		      tempskp1  = evalin('caller','datakp1.gene.temps');
		      kind      = NaN;
		   catch
		      ok=0;
		   end
		if ok == 1
			data   = cat(1,datak,datakp1);
			temps  = cat(1,tempsk,tempskp1);
			rho    = cat(1,rhok,rhokp1);
		end
	else
		if compare
			try
                    jeux1 = double(evalin('base',strcat('jeux1.',nom)));
                 end
            end
	    if ~isempty(jeux1)
                temps1 = evalin('base','jeux1.data.gene.temps');
                x1     = evalin('base','jeux1.param.gene.x');
            end

	end

	if ok == 0
		ok   = 1;
		mode =1;
		try
		   datak     = double(evalin('caller',nomk));
	           rhok      = evalin('caller','datak.equi.rhoRZ');
	           if size(datak,2)>1
	           	x    = linspace(0,1,size(datak,2));
	           else
	           	x    = linspace(0,1,21);
	           end
		   tempsk    = evalin('caller','datak.gene.temps');
		   kind      = NaN;
		catch
		   ok=0;
		end
		if ok == 1
			data   = cat(1,datak,datak);
			temps  = cat(1,tempsk,tempsk);
			rho    = cat(1,rhok,rhok);
		end
	end

	if ok == 0
		ok   = 1;
		mode =1;
		try
		   datak     = double(evalin('caller',nomkp1));
	           rhok      = evalin('caller','datakp1.equi.rhoRZ');
	           if size(datak,2)>1
	           	x    = linspace(0,1,size(datak,2));
	           else
	           	x    = linspace(0,1,21);
	           end
		   tempsk    = evalin('caller','datakp1.gene.temps');
		   kind      = NaN;
		catch
		   ok=0;
		end
		if ok == 1
			data   = cat(1,datak,datak);
			temps  = cat(1,tempsk,tempsk);
			rho    = cat(1,rhok,rhok);
		end
	end

	if isempty(data)
		data=[];
		temps=[];
		kind=[];
		x=[];
	end

	disp('-------------------------------------------')
	fprintf('Variable name : %s\n',ztr(nom));
	fprintf('Help : %s\n',aide);
	set(hf,'name',sprintf('Help : %s',aide));
	if mode == 0
		disp('time dependent')
	else
		disp('for 2 time slices')
		fprintf('k = %d\n',kind);
		fprintf('time = %g et %g\n',temps(1),temps(2));
	end

	if isstruct(data)
		disp(' ')
		disp(data)
		disp( ' ')
	else
		fprintf('data size :')
		disp(size(data))
		disp(' ')
	end

	indnan=find(~isfinite(data));
	if ~isempty(indnan)
	    fprintf('Contains NaN or Inf (#%d)\n\n',length(indnan));
	end
	indimag=find(imag(data));
	if ~isempty(indimag)
	    fprintf('Contains complex numbers (#%d)\n\n',length(indimag));
	end

	sld =size(data);

	indl = find(nom =='.');
	if isempty(indl)
		noml = nom;
	else
		noml = nom((min(indl) +1):end);
	end

	% plot sur le temps si length(sld) =2
	if length(sld) == 2
		if any(imag(data(:)))
			h=findobj(hf,'type','axes','tag','axes_des_phases');
			hh=findobj(hf,'type','uimenu','tag','efface_temps');
			ck = lower(get(hh,'checked'));
			if isempty(h)|strcmp(ck,'off')
				h=subplot(4,1,2);
				set(h,'tag','axes_des_phases');
				hold off
			else
				axes(h);
			end
			if length(temps) >2
				if size(data,2)>1
					if size(data,2)~=length(x)
					        if modeimag == 1
						    hp = plot(temps,imag(data),'-');
						elseif modefce == 1
						    [tor,pol] = decodefce(data);
						    hp = plot(temps,tor,'-');
						    hold on
						    for kp = 1:length(hp)
						    	plot(temps,pol(:,kp),'linestyle','--','color',get(hp(kp),'color'));
						    end
						else
						    hp = plot(temps,angle(data),'-');
						end
						if ~isempty(jeux1)
							hold on
							if size(data,2) == size(jeux1,2)
							   cmt = get(hp,'color');
							   for k = 1:size(jeux1,2)
							       	if modeimag == 1
									plot(temps1,imag(jeux1(:,k)),'linestyle',':','color',cmt{k});
								elseif modefce == 1
						    			[tor,pol] = decodefce(jeux1(:,k));
									plot(temps1,tor,'linestyle',':','color',cmt{k});
									plot(temps1,pol,'linestyle','-.','color',cmt{k});

								else
									plot(temps1,angle(jeux1(:,k)),'linestyle',':','color',cmt{k});
								end
							   end
							else
					                     if modeimag == 1
							          plot(temps1,imag(jeux1),'linestyle',':');

							    elseif modefce == 1
						    			[tor,pol] = decodefce(jeux1);
									plot(temps1,tor,'linestyle',':','color',get(hp,'color'));
									plot(temps1,pol,'linestyle','-.','color',get(hp,'color'));
							     else
							          plot(temps1,angle(jeux1),'linestyle',':','color',get(hp,'color'));
							     end
							end
						end
					elseif modeimag == 1
						hp1=plot(temps,imag(data(:,1)),'-');
						hold on
						hp2=plot(temps,imag(data(:,end)),':');
						if ~isempty(jeux1)
							plot(temps1,imag(jeux1(:,1)),'linestyle','none','marker','o','color',get(hp1,'color'));
							plot(temps1,imag(jeux1(:,end)),'linestyle','none','marker','d','color',get(hp2,'color'));
						end
					elseif modefce == 1
						[tor,pol] = decodefce(data(:,1));
						hp1=plot(temps,tor,'-');
						hold on
						hp1=plot(temps,pol,'linestyle','--','color',get(hp1,'color'));
						[tor,pol] = decodefce(data(:,end));
						hp2=plot(temps,tor,':');
						hp2=plot(temps,pol,'linestyle','_.','color',get(hp2,'color'));
						if ~isempty(jeux1)
							[tor,pol] = decodefce(jeux1(:,1));
							plot(temps1,tor,'linestyle','none','marker','o','color',get(hp1,'color'));
							plot(temps1,pol,'linestyle','none','marker','s','color',get(hp1,'color'));
							[tor,pol] = decodefce(jeux1(:,end));
							plot(temps1,tor,'linestyle','none','marker','d','color',get(hp2,'color'));
							plot(temps1,pol,'linestyle','none','marker','h','color',get(hp2,'color'));

						end
					else
						hp1=plot(temps,angle(data(:,1)),'-');
						hold on
						hp2=plot(temps,angle(data(:,end)),':');
						if ~isempty(jeux1)
							plot(temps1,angle(jeux1(:,1)),'linestyle','none','marker','o','color',get(hp1,'color'));
							plot(temps1,angle(jeux1(:,end)),'linestyle','none','marker','d','color',get(hp2,'color'));
						end
					end
				elseif modeimag == 1
					hp=plot(temps,imag(data));
					if ~isempty(jeux1)
						hold on
						plot(temps1,imag(jeux1),'color',get(hp,'color'));
					end
				elseif modefce == 1
					[tor,pol] = decodefce(data);
					hp=plot(temps,tor,'-');
					hold on
					hp=plot(temps,pol,'linestyle','--','color',get(hp,'color'));

					if ~isempty(jeux1)
						hold on
						[tor,pol] = decodefce(jeux1);
						hp = plot(temps1,tor,'-.');
						plot(temps1,pol,'linestyle',':','color',get(hp,'color'));
					end
				else
					hp=plot(temps,angle(data));
					if ~isempty(jeux1)
						hold on
						plot(temps1,angle(jeux1),'color',get(hp,'color'));
					end
				end
			elseif modeimag == 1
				if size(data,2)>1
					if (size(data,2)~=length(x))&(size(data,2)<=7)
						plot(temps,imag(data,'o'));
					else
						plot(temps,imag(data(:,1)),'-');
						hold on
						plot(temps,imag(data(:,end)),':');
					end
				else
					plot(temps,angle(data),'o');
				end
			elseif modefce == 1
				if size(data,2)>1
					if (size(data,2)~=length(x))&(size(data,2)<=7)
					        [tor,pol] = decodefce(data);
						hp= plot(temps,tor,'o');
						hold on
						plot(temps,tor,'linestyle','none','marker','d','color',get(hp,'color'));
					else
					        [tor,pol] = decodefce(data(:,1));
						hp= plot(temps,tor,'-');
						hold on
						plot(temps,tor,'linestyle','--','color',get(hp,'color'));
					        [tor,pol] = decodefce(data(:,end));
						hp= plot(temps,tor,':');
						plot(temps,tor,'linestyle','-.','color',get(hp,'color'));
					end
				else
					plot(temps,angle(data),'o');
				end
			else
				if size(data,2)>1
					if (size(data,2)~=length(x))&(size(data,2)<=7)
						plot(temps,angle(data),'o');
					else
						plot(temps,angle(data(:,1)),'-');
						hold on
						plot(temps,angle(data(:,end)),':');
					end
				else
					plot(temps,angle(data),'o');
				end
			end
			hold on
			if  modeimag == 1
			      if modeylabel == 0
			           ylabel(sprintf('d(%s)/dt',noml));
			      else
			           ylabel(sprintf('imag(%s)',noml));
					end
			else
			      ylabel(sprintf('angle(%s)',noml));
			end
			co = get(h,'colororder');
			cog=co(1,:);
			set(h,'colororder',co(cat(2,2:size(co,1),1),:));
			set(h,'tag','axes_des_phases','xticklabel','');
		end
		h=findobj(hf,'type','axes','tag','axes_des_temps');
		hh=findobj(hf,'type','uimenu','tag','efface_temps');
		ck = lower(get(hh,'checked'));
		if isempty(h)|strcmp(ck,'off')
			if any(imag(data(:)))
				h=subplot(4,1,1);
				if  modeimag == 1
					data =real(data);
					jeux1=real(jeux1);
				else
					data =abs(data);
					jeux1=abs(jeux1);
				end
			else
				h=subplot(2,1,1);
			end
			hold off
			set(h,'tag','axes_des_temps');
		else
			axes(h);
			if any(imag(data(:)))
				if  modeimag == 1
					data =real(data);
					jeux1=real(jeux1);
				else
					data =abs(data);
					jeux1=abs(jeux1);
				end
				%data =abs(data);
				%jeux1=abs(jeux1);
			end
		end
		if length(temps) >2
			if size(data,2)>1
				if size(data,2)~=length(x)
					hp=plot(temps,data,'-');
					if ~isempty(jeux1)
						hold on
						if size(data,2) == size(jeux1,2)
						    cmt = get(hp,'color');
						    for k = 1:size(jeux1,2)
							 plot(temps1,jeux1(:,k),'linestyle',':','color',cmt{k});
					            end
						else
						    plot(temps1,jeux1,'linestyle',':');
						end
					end
				else
					hp1=plot(temps,data(:,1),'-');
					hold on
					hp2=plot(temps,data(:,end),':');
					if ~isempty(jeux1)
						plot(temps1,jeux1(:,1),'color',get(hp1,'color'),'linestyle','none','marker','o');
						plot(temps1,jeux1(:,end),'color',get(hp2,'color'),'linestyle','none','marker','d');
					end
				end
			else
				hp=plot(temps,data);
				if ~isempty(jeux1)
					hold on
					plot(temps1,jeux1,'color',get(hp,'color'),'linestyle',':');
				end
			end
		else
			if size(data,2)>1
				if (size(data,2)~=length(x))&(size(data,2)<=7)
					plot(temps,data,'o');
				else
					plot(temps,data(:,1),'-');
					hold on
					plot(temps,data(:,end),':');
				end
			else
				plot(temps,data,'o');
			end
		end
		hold on
		xlabel('Time (s)');
		if modeylabel == 0
			   ylabel(ztr(noml));
		else
			   ylabel(sprintf('real(%s)',ztr(noml)));
		end
      		%ylabel(noml);
		title(['f(t)    ' desctitre])
		co = get(h,'colororder');
		cog=co(1,:);
		set(h,'colororder',co(cat(2,2:size(co,1),1),:));
		set(h,'tag','axes_des_temps');

		if size(data,2)==length(x)
			% les profils
			h=findobj(hf,'type','axes','tag','axes_des_profils');
			hh=findobj(hf,'type','uimenu','tag','efface_profil');
			ck = lower(get(hh,'checked'));
			if isempty(h)|strcmp(ck,'off')
				h=subplot(2,3,4);
				hold off
				plot(NaN,NaN);
				set(h,'tag','axes_des_profils');
			else
				axes(h);
			end
			if length(temps) >2
				zplotprof(h,temps,x,data,'color',cog);
				if ~isempty(jeux1)
					zplotprof(h,temps1,x1,jeux1,'color',cog,'linestyle','none','marker','o');
				end
			else
				plot(x,data(1,:),'linestyle','-','marker','none','color',cog);
				hold on
				plot(x,data(2,:),'linestyle',':','marker','none','color',cog);
			end
			hold on
			xlabel('r','Fontname','Symbol')
			ylabel(ztr(noml));
			title('Profiles')
			set(h,'tag','axes_des_profils');

			hm=findobj(hf,'type','uimenu','tag','log_image');
			if ~isempty(hm)
				chk = strcmp(get(hm,'checked'),'on');
			else
				chk =0;
			end
			if chk
				datas =log(abs(data))./log(10);
				tnoml = sprintf('log_1_0(%s)',noml);
			else
				datas =data;
				tnoml =noml;
			end

			h=subplot(2,3,5);
			zimagesc(x,temps,datas)
			colormap('default')
			colorbar
			set(gca,'ydir','normal')
			xlabel('r','Fontname','Symbol')
			ylabel('Time (s)')
			title(ztr(tnoml))
			set(h,'tag','axes_des_images');
		elseif size(data,2) >= 2
			% les donnees 2D (temps,indice)
			h=findobj(hf,'type','axes','tag','axes_des_profils');
			hh=findobj(hf,'type','uimenu','tag','efface_profil');
			ck = lower(get(hh,'checked'));
			if isempty(h)|strcmp(ck,'off')
				h=subplot(2,3,4);
				hold off
				plot(NaN,NaN);
				set(h,'tag','axes_des_profils');
			else
				axes(h);
			end
			xx =  ones(size(temps,1),1) * (1:size(data,2));
			if length(temps) >2
				zplotprof(h,temps,xx,data,'color',cog,'marker','*','linestyle','none');
				if ~isempty(jeux1)
					xx1 =  ones(size(temps1,1),1) * (1:size(jeux1,2));
					zplotprof(h,temps1,xx1,jeux1,'color',cog,'linestyle','none','marker','o');
				end
			else
				plot(xx,data(1,:),'linestyle','-','marker','none','color',cog);
				hold on
				plot(xx,data(2,:),'linestyle',':','marker','none','color',cog);
			end
			hold on
			xlabel('index')
			ylabel(ztr(noml));
			title('Profiles')
			set(h,'tag','axes_des_profils');

			hm=findobj(hf,'type','uimenu','tag','log_image');
			if ~isempty(hm)
				chk = strcmp(get(hm,'checked'),'on');
			else
				chk =0;
			end
			if chk
				datas =log(abs(data))./log(10);
				tnoml = sprintf('log_1_0(%s)',noml);
			else
				datas =data;
				tnoml =noml;
			end

			h=subplot(2,3,5);
			zimagesc(xx(1,:),temps,datas)
			colormap('default')
			colorbar
			set(gca,'ydir','normal')
			xlabel('index')
			ylabel('Time (s)')
			title(ztr(tnoml))
			set(h,'tag','axes_des_images');

		end
	else
		% debut plot moyenne
		if any(imag(data(:)))
			h=findobj(hf,'type','axes','tag','axes_des_phases');
			hh=findobj(hf,'type','uimenu','tag','efface_temps');
			ck = lower(get(hh,'checked'));
			if isempty(h)|strcmp(ck,'off')
				h=subplot(4,1,2);
				set(h,'tag','axes_des_phases');
				hold off
			else
				axes(h);
			end
			if length(temps) >2
				hp = plot(temps,angle(squeeze(mean(data,2))),'-');
				if ~isempty(jeux1)
					hold on
					plot(temps1,angle(squeeze(mean(jeux1,2))),'linestyle',':','color',get(hp,'color'));
				end
			end
			hold on
			ylabel(sprintf('angle(<%s>)',ztr(noml)));
			co = get(h,'colororder');
			cog=co(1,:);
			set(h,'colororder',co(cat(2,2:size(co,1),1),:));
			set(h,'tag','axes_des_phases','xticklabel','');
		end
		h=findobj(hf,'type','axes','tag','axes_des_temps');
		hh=findobj(hf,'type','uimenu','tag','efface_temps');
		ck = lower(get(hh,'checked'));
		if isempty(h)|strcmp(ck,'off')
			if any(imag(data(:)))
				h=subplot(4,1,1);
				data =abs(data);
				jeux1=abs(jeux1);
			else
				h=subplot(2,1,1);
			end
			hold off
			set(h,'tag','axes_des_temps');
		else
			axes(h);
			if any(imag(data(:)))
				data =abs(data);
				jeux1=abs(jeux1);
			end
		end
		if length(temps) >2
			hp=plot(temps,squeeze(mean(data,2)),'-');
			if ~isempty(jeux1)
				hold on
				for k = 1:size(jeux1,3)
				    plot(temps1,squeeze(mean(jeux1(:,:,k),2)),'color',get(hp(k),'color'),'linestyle',':');
				end
			end
		end
		hold on
		xlabel('Time (s)');
		ylabel(sprintf('<%s>',ztr(noml)));
		title(['f(t)    ' desctitre])
		co = get(h,'colororder');
		cog=co(1,:);
		set(h,'colororder',co(cat(2,2:size(co,1),1),:));
		set(h,'tag','axes_des_temps');


		% fin plot moyene
		% cas a trois dimensions
		if size(data,2) == length(x)
			yy = ones(size(temps,1),1)*x;
			xlab = 'x (su)';
		elseif size(data,2) == length(rho)
			yy = double(rho);
			xlab ='rho (m)';
		else
			yy = ones(size(temps,1),1)*(1:size(data,2));
			xlab = 'index';
		end
		tt = temps(:) * ones(1,size(yy,2));
		if ~isempty(jeux1)
			if size(jeux1,2) == length(x1)
				yy1 = ones(size(temps1,1),1)*x1;
			elseif size(jeux1,2) == length(rho)
				yy1 = double(rho);
			else
				yy1 = ones(size(temps1,1),1)*(1:size(jeux1,2));
			end
			tt1 = temps1(:) *ones(1,size(yy1,2));
		end
		% plot des profils
		h=subplot(2,2,3);
		hold off
		plot(NaN,NaN);
		set(h,'tag','axes_des_profils');
		for l =1:size(data,3);
			co = get(h,'colororder');
			cog=co(1,:);
			set(h,'colororder',co(cat(2,2:size(co,1),1),:));
			zplotprof(h,tt,yy,double(data(:,:,l)),'color',cog);
			if ~isempty(jeux1)
			     zplotprof(h,tt1,yy1,double(jeux1(:,:,l)),'color',cog,'linestyle',':');
			end
			fprintf('.');
		end
		fprintf('\n');
		xlabel(xlab);
		set(h,'tag','axes_des_profils');
		ylabel(ztr(noml));
		title('Profiles')
	end
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
	h=findobj(hf,'type','axes','tag','axes_des_temps');
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
	h=findobj(hf,'type','axes','tag','axes_des_temps');
 	hm=findobj(hf,'type','uimenu','tag','log_temps');
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
	h=findobj(hf,'type','axes','tag','axes_des_profils');
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
	h=findobj(hf,'type','axes','tag','axes_des_profils');
 	hm=findobj(hf,'type','uimenu','tag','log_profil');
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
	    R     = double(evalin('base','data.equi.R'));
	    Z     = double(evalin('base','data.equi.Z'));
	    Rext  = double(evalin('base','data.geo.R'));
	    Zext  = double(evalin('base','data.geo.Z'));
	    temps = evalin('base','data.gene.temps');
	catch
	    ok=0;
	end
	if ok == 0
		ok   = 1;
		mode =1;
		try
		    Rk    = double(evalin('caller','datak.equi.R'));
		    Zk    = double(evalin('caller','datak.equi.Z'));
		    Rkp1  = double(evalin('caller','datakp1.equi.R'));
		    Zkp1  = double(evalin('caller','datakp1.equi.Z'));
		    tempsk    = evalin('caller','datak.gene.temps');
		    tempskp1  = evalin('caller','datakp1.gene.temps');
		    Rextk  = double(evalin('caller','datak.geo.R'));
		    Zextk  = double(evalin('caller','datak.geo.Z'));
		    Rextkp1  = double(evalin('caller','datakp1.geo.R'));
		    Zextkp1  = double(evalin('caller','datakp1.geo.Z'));
		catch
		   ok=0;
		end
		if ok == 1
			R   = cat(1,Rk,Rkp1);
			Z   = cat(1,Zk,Zkp1);
			Rext   = cat(1,Rextk,Rextkp1);
			Zext   = cat(1,Zextk,Zextkp1);
			temps  = cat(1,tempsk,tempskp1);
			rho    = cat(1,rhok,rhokp1);
		end
	else
	  	if compare
	   		try
	   		      R1     = double(evalin('base','jeux1.data.equi.R'));
	   		      Z1     = double(evalin('base','jeux1.data.equi.Z'));
	   		      Rext1  = double(evalin('base','jeux1.data.geo.R'));
	   		      Zext1  = double(evalin('base','jeux1.data.geo.Z'));
	   		      temps1 = evalin('base','jeux1.data.gene.temps');
	   		end
   	  end
	 end
	if ok == 0
		ok   = 1;
		mode =1;
		try
		    Rk    = double(evalin('caller','datak.equi.R'));
		    Zk    = double(evalin('caller','datak.equi.Z'));
		    Rkp1  = double(evalin('caller','datak.equi.R'));
		    Zkp1  = double(evalin('caller','datak.equi.Z'));
		    tempsk    = evalin('caller','datak.gene.temps');
		    tempskp1  = evalin('caller','datak.gene.temps');
	            Rextk  = double(evalin('caller','datak.geo.R'));
	            Zextk  = double(evalin('caller','datak.geo.Z'));
	            Rextkp1  = double(evalin('caller','datak.geo.R'));
	            Zextkp1  = double(evalin('caller','datak.geo.Z'));
		catch
		   ok=0;
		end
		if ok == 1
			R   = cat(1,Rk,Rkp1);
			Z   = cat(1,Zk,Zkp1);
			Rext   = cat(1,Rextk,Rextkp1);
			Zext   = cat(1,Zextk,Zextkp1);
			temps  = cat(1,tempsk,tempskp1);
			rho    = cat(1,rhok,rhokp1);
		end

	end
	if ok == 0
		ok   = 1;
		mode =1;
		try
		    Rk    = double(evalin('caller','datakp1.equi.R'));
		    Zk    = double(evalin('caller','datakp1.equi.Z'));
		    Rkp1  = double(evalin('caller','datakp1.equi.R'));
		    Zkp1  = double(evalin('caller','datakp1.equi.Z'));
		    tempsk    = evalin('caller','datakp1.gene.temps');
		    tempskp1  = evalin('caller','datakp1.gene.temps');
	            Rextk  = double(evalin('caller','datakp1.geo.R'));
	            Zextk  = double(evalin('caller','datakp1.geo.Z'));
	            Rextkp1  = double(evalin('caller','datakp1.geo.R'));
	            Zextkp1  = double(evalin('caller','datakp1.geo.Z'));
		catch
		   ok=0;
		end
		if ok == 1
			R   = cat(1,Rk,Rkp1);
			Z   = cat(1,Zk,Zkp1);
			Rext   = cat(1,Rextk,Rextkp1);
			Zext   = cat(1,Zextk,Zextkp1);
			temps  = cat(1,tempsk,tempskp1);
			rho    = cat(1,rhok,rhokp1);
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
	if mode == 0
		disp('time dependent')
	else
		disp('for 2 time slices')
		fprintf('k = %d\n',kind);
		fprintf('time = %g et %g\n',temps(1),temps(2));
	end

	fprintf('data size :')
	disp(size(R))
	disp(' ')

	% les axes
	h=subplot(2,3,6);
	hold off
	plot(NaN,NaN);
	set(h,'tag','axes_des_surfaces');

	% limitationdu nombre de surfaces
	pas =fix(size(R,2)./21);
	if pas ==0
		pas=1;
	end
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
	   paroi = evalin('base','param.from.paroi');
	   hold on
	   if ~isempty(paroi)
	   	plot(paroi.R,paroi.Z)
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

% plot profile
if strcmp(cmd,'plot_profile')
   info = evalin('base','param.profile.data');
   zplotprofile(info);
end

% rapport profile
if strcmp(cmd,'rapport_profile')
   info = evalin('base','param.profile.data');
   profreport(info);
end



if strcmp(cmd,'NaNInf')
	%
	ok   = 1;
	mode = 0;
	param  = evalin('base','param');
	nbt    = param.gene.nbt;
	try
	    data  = evalin('base','data');
            indk  = param.gene.kmin:param.gene.k;
	catch
	    ok=0;
	end
	if ok == 0
		ok   = 1;
		mode =1;
		try
		   datak     = evalin('caller','datak');
		   datakp1   = evalin('caller','datakp1');
                   indk = 1:2;
		catch
		   ok=0;
		end
		if ok == 1
			data   = cat(1,datak,datakp1);
		end
	end

	if ok == 0
		ok   = 1;
		mode =1;
		try
		   datak     = evalin('caller','datak');
                   indk = 1;
		catch
		   ok=0;
		end
		if ok == 1
			data   = cat(1,datak,datak);
		end
	end

	if ok == 0
		ok   = 1;
		mode =1;
		try
		   datak     = evalin('caller','datakp1');
                   indk = 1;
		catch
		   ok=0;
		end
		if ok == 1
			data   = cat(1,datak,datak);
		end
	end
	if ~isempty(data)
	   % ouverture du fichier
	   file = tempname;
	   fid =fopen(file,'w');

	   % 2 - la structure data
	   fprintf(fid,'---------------------------------\n');
	   fprintf(fid,'Search of NaN and Inf :\n');
	   % liste des champs la structures
	   champ = sort(fieldnames(data));
	   for k=1:length(champ)
		   % nom complet pour acceder a la variable
		   champ{k}=strcat('data.',champ{k});
	   end

	   % jusqu'a ce qu'il n'y ait plus de champ
	   test=0;
	   while(~isempty(champ))
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
	   	   champ=cat(1,champ,champnew);
	      else
	   	   % juste pour les tests
	   	   %disp(champc);
	   	   %fprintf('.');
		   var = double(eval(champc,[]));
                   ss  = size(var);
                   if (ss(1) == 1) & (ss(2) > 1)
                     var = var';
                     ss  = size(var);
                  end
		   if  isempty(var)
                         indnan =[];
                         nb =0;
                   elseif length(ss) == 2
                         indnan = find(~isfinite(var(indk,:)));
		         nb     = prod(size(var(indk,:)));
                   elseif length(ss) == 3
                         indnan = find(~isfinite(var(indk,:,:)));
		         nb     = prod(size(var(indk,:,:)));
                   else
                         indnan = find(~isfinite(var(indk,:,:,:)));
		         nb     = prod(size(var(indk,:,:,:)));
                   end
		   if ~isempty(indnan)
		      ll = length(indnan);
		      if ll == nb
		          fprintf(fid,'found in %s some NaN or Inf in\n',champc);
		      else
		          fprintf(fid,'found in %s  #%d NaN or Inf\n',champc,ll);
		      end
		   end
		   if size(var,1) ~= nbt
		          fprintf(fid,'%s wrong time size \n',champc);
		   end
	      end
           end
        end
	fclose(fid);
	[s,t]=unix(['ne ',file,' &']);
	fprintf('file %s created\n',file);
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
	hnum    = findobj(hf,'type','uicontrol','tag','numchoc');
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
	htop10 = findobj(hf,'type','uimenu','tag','top10ts');
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
	h=findobj(hf,'type','axes','tag','axes_des_temps');
	hh=findobj(hf,'type','uimenu','tag','efface_temps');
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
		h=findobj(hf,'type','axes','tag','axes_des_profils');
		hh=findobj(hf,'type','uimenu','tag','efface_profil');
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

		hm=findobj(hf,'type','uimenu','tag','log_image');
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

		h=subplot(2,3,5);
		zimagesc(rhofit,temps,datas)
		colormap('default')
		colorbar
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

if strcmp(cmd,'cree0D')
	if ~isempty(findobj(hf,'type','uimenu','tag','0Dmenu'))
		return
	end
	% lecture des informations
        try
            zs  = evalin('base','post');
        catch
            zs = [];
        end
	if isempty(zs)
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

end

if strcmp(cmd,'fast')
	hgsave(gcf,fullfile(getenv('HOME'),'/zineb/','fastzdataplot'));
        fprintf('%s cree ...\n',fullfile(getenv('HOME'),'/zineb/','fastzdataplot'));
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
	if isempty(nomvar)
		% liste des mots reservee
		s =load('zineb_mots_reserves.mat','reserve');
		reserve = s.reserve;
		clef ={};
		for k = 3:length(reserve)
			if isempty(findstr(reserve{k},'.'))
				if isempty(strmatch(reserve{k},clef,'exact'))
					clef{end+1} = reserve{k};
				end
			end
		end
		clef =sort(clef);
		[s,v] = listdlg('PromptString','Choose a key :',...
	              		'SelectionMode','single',...
		            	'ListString',clef);
		if v == 0 | isempty(s)
			return
		else
			nomvar =clef{s};
		end
	end

	% essai de lecture
	s =load('zineb_mots_reserves.mat','dico');
	dico = s.dico;
	if isfield(dico,nomvar)
		info = getfield(dico,nomvar);
		if length(info) >1
			[s,v] = listdlg('PromptString','chosse a variable name :',...
					'SelectionMode','single',...
					'ListString',info);
			if v == 0 | isempty(s)
				return
			else
				nomvar =info{s};
			end

		elseif iscell(info)
			nomvar =info{1};
		else
			nomvar =info;
		end
		if ~isempty(nomvar)
	    		hcb =findobj(hf,'type','uimenu','tag',nomvar);
	    		if ~isempty(hcb)
	    			zdataplot(strtok(nomvar,'.'),hcb);
	    		end
		end
	else
		warndlg('this variable is undocumented ...','unlucky !');
		return
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
    hot   = findobj(hf,'type','uicontrol','style','edit','tag','temps');
    if ~isempty(ho)
        set(hot,'string',sprintf('%g',x));
	zplotprof('temps');
    end
    set(ho,'value',0);
end


% supperposition
if strcmp(cmd,'supperpose')
      hax  = findobj(hf,'type','axes','tag','axes_des_profils');
      if isempty(hax)
         set(ho,'value',0);
         return
      end
      axes(hax);
      dy   = abs(diff(get(hax,'ylim')));
      dx   = abs(diff(get(hax,'xlim')));
      hll  = findobj(hax,'type','line','tag','zplotprof');
      if isempty(hll)
          set(ho,'value',0);
          return
      end
      proper = set(hll(1));
      mlist  = proper.Marker(6:(end-1));
      hlm    = findobj(hax,'type','line');
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
            hot   = findobj(hf,'type','uicontrol','style','edit','tag','temps');
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
	fprintf('Nmae of the data : %s [%s@%g]\n',nom,machine,numchoc);
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
	h=findobj(hf,'type','axes','tag','axes_des_temps');
	hh=findobj(hf,'type','uimenu','tag','efface_temps');
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
	fprintf('Nmae of the data : %s [%s@%g]\n',nom,machine,numchoc);
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
	h=findobj(hf,'type','axes','tag','axes_des_temps');
	hh=findobj(hf,'type','uimenu','tag','efface_temps');
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
	h=findobj(hf,'type','axes','tag','axes_des_temps');
	hh=findobj(hf,'type','uimenu','tag','efface_temps');
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
		h=findobj(hf,'type','axes','tag','axes_des_profils');
		hh=findobj(hf,'type','uimenu','tag','efface_profil');
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

		hm=findobj(hf,'type','uimenu','tag','log_image');
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

		h=subplot(2,3,5);
		zimagesc(mean(xli,1),temps,datas)
		colormap('default')
		colorbar
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
				uimenu(h,'label',com,'tag',nomc,'userdata',com,'callback','zdataplot(''dataTS'')');
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
			uimenu(h,'label',com,'tag',nomc,'userdata',com,'callback','zdataplot(''dataTS'')');
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
		uimenu(h,'label',com,'tag',nomc,'userdata',com,'callback','zdataplot(''data0d'')');
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
		uimenu(h,'label',com,'tag',nomc,'userdata',com,'callback','zdataplot(''exp0d'')');
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
		uimenu(h,'label',com,'tag',nomc,'userdata',com,'callback','zdataplot(''data0d1d'')');
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

vv = interp1(x',v',xx','nearest')';
indif = find(diff(t) <=0);
nbmax =100;
if ~isempty(indif)
	while( ~isempty(indif) & (nbmax >0))
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
if ~isempty(indund)
    nom = strcat(nom(1:indund-1),'_{',nom(indund+1:end),'}');
end
