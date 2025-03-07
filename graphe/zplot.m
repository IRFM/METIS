% ZPLOT plot interactif de Cronos-U (zineb)
%-------------------------------------------
% fichier zplot.m ->  zplot
%
%
% fonction Matlab 5 :
%
% Cette fonction calcule la puissance de chauffage du aux alpha et
% les sources de matiere du aux produits de fusion. 
%  
% syntaxe  :
%  
%     [cr,info]=zplot(data,gene,from,info,cons,phys,dk);
%    
% entrees :
%
%     data    =  structure des donnees temporelles
%     gene    =  parametres generique (param.gene)
%     from    =  description des sources de donnees (param.from)
%     info    =  information sur le plot (param.plot)
%     cons    =  parametres de la fonction (param.cons.plot)
%     phys    =  constantes physiques (param.phys)
%     dk      =  nombre de pas de temps a tracer
%
% sorties :
% 
%     cr      =  compte rendu d'execution
%     info    =  informations sur le plot modifiees (param.plot)
% 
% parametres : aucun
% 
%
% fonction ecrite par J-F Artaud , poste 46-78
% version 1.0, du 02/03/2001
% 
% 
% liste des modifications : 
%
%--------------------------------------------------------------
%
function [cr,info]=zplot(data,gene,from,info,cons,phys,dk)

% initialisation de cr
cr =0;

% mode initialisation 
% fonction auto declarante                             
if nargin <=1 

	valeur     = [];            % pas de parametres
	type       = [];
	borne      = [];
	defaut     = [];
	info       = [];
	
	interface.ts        = '';   % pas d'interface  a definir (nom de la fonction d'interfacage avec les donnees TS)
	interface.jet       = '';   % pas d'interface  a definir (nom de la fonction d'interfacage avec les donnees Jet)
	
	sortie.valeur=valeur;
	sortie.type=type;
	sortie.borne=borne;
	sortie.defaut=defaut;
	sortie.info=info;
	sortie.interface=interface;

	sortie.description = 'Plot interratif pour Cronos-Zineb';   % description (une ligne) de la fonction
	sortie.help = '';                            % nom du fichier d'aide s'il existe, sinon aide de la fonction
	sortie.gui  ='';                             % nom de l'interface graphique specifique si elle existe
	sortie.controle = '';                        % nom de la fonction de controle des valeurs si elle existe
	
	cr=sortie;
	
	return

end


% test des handles et creation des figures si necessaire
% figure des controles
if ~ishandle(info.figure1)|isempty(info.figure1)
	cree = 1;
elseif ~strcmp(get(info.figure1,'tag'),'cz_controle') 
	cree = 1;
else
	cree = 0;
end

if cree == 1
	
	h = findobj(0,'type','figure','tag','cz_controle');
	if ~isempty(h)
		delete(h);
	end
	
	
	info.figure1 = figure('name','Cronos-Zineb controle','tag','cz_controle','number','off', ...
	                      'units','normalized','menubar','none','interruptible','on', ...
	                      'busyaction','queue');
	                      
	drawnow;                      
	set(info.figure1,'position',[0.01,0,0.5,0.15]);                      
	% les boutons
	info.h.run  = uicontrol('style','radio','tag','cz_run','string','Run', ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0,0,0.15,0.15], ...
	                      'value',info.run,'tooltipstring','Execution de la simulation', ...
	                      'callback','zrun');                      
	                      
	info.h.pause = uicontrol('style','radio','tag','cz_pause','string','Pause', ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0.17,0,0.15,0.15], ...
	                      'tooltipstring','Suspend l''execution de la simulation', ...
	                      'value',~info.run,'callback','zpause');                     
	                      
	cmd = ['ho=gco;', ...
	       'zverbose(''Keyboard => '');', ...
	       'disp(''Tapez ''''return'''' pour reprendre l''''execution de Zineb'');', ...
	       'set(ho,''value'',0);', ...
	       'keyboard'];
               
	info.h.keyboard = uicontrol('style','radio','tag','cz_keyboard','string','Keyboard', ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0.34,0,0.15,0.15], ...
	                      'tooltipstring','Donne la main dans la fenetre de commande matlab', ...
	                      'value',0,'callback',cmd);  
	                      
	info.h.step = uicontrol('style','radio','tag','cz_step','string','Step', ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0.51,0,0.15,0.15], ...
	                      'tooltipstring','Execution pas a pas (en mode pause)', ...
	                      'value',0);  
	cmd = ['ho=gco;set(gco,''value'',1);', ...
	       'b=questdlg(''Confirmez-vous l''''arret definitif de la simulation'',''Attention -> arret definitif'',''Oui'',''Non'',''Non'');', ...
	       'if ~strcmp(b,''Oui''),set(ho,''value'',0);end'];
	info.h.fin = uicontrol('style','radio','tag','cz_fin','string','Fin', ...
	                      'units','normalized','interruptible','off', ...
	                      'busyaction','queue','position',[0.8,0,0.2,0.15], ...
	                      'backgroundcolor',[1 0 0],'tooltipstring','Arret definitif de la simulation', ...
	                      'value',info.fin,'callback',cmd);  
	                      
	txf = 'Tout ce que je te dis, il faut le croire !';
	
	info.h.info = uicontrol('style','text','tag','cz_info','string',txf, ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0,0.2,1,0.15], ...
	                      'backgroundcolor',[1 1 1],'tooltipstring','Zone d''information');                       
	                      
	 uicontrol('style','text','tag','cz_input_text','string','Fichier d''entree : ', ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0,0.85,0.2,0.15]); 
	                      
	 uicontrol('style','text','tag','cz_input','string',gene.origine, ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0.2,0.85,0.8,0.15]);                       
	                      
	 uicontrol('style','text','tag','cz_output_text','string','Fichier de sortie : ', ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0,0.69,0.2,0.15]); 
	                      
	 uicontrol('style','text','tag','cz_output','string',gene.file, ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0.2,0.69,0.8,0.15]);                       
	                      
	 uicontrol('style','text','tag','cz_machine_text','string','Machine', ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0,0.5,0.245,0.12]); 
	                      
	 uicontrol('style','text','tag','cz_num_text','string','Choc #', ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0.25,0.5,0.205,0.12]);                       
	                      
	 uicontrol('style','text','tag','cz_num_text','string','Occurrence', ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0.46,0.5,0.14,0.12]);                       
	                      
	 uicontrol('style','text','tag','cz_date_text','string','Date', ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0.605,0.5,0.395,0.12]); 
	                      
	 uicontrol('style','text','tag','cz_machine_text','string',from.machine, ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0,0.37,0.245,0.12]); 
	                      
	 uicontrol('style','text','tag','cz_num_text','string',int2str(from.shot.num), ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0.25,0.37,0.205,0.12]);                       
	                      
	 uicontrol('style','text','tag','cz_num_text','string',int2str(fix(10*rem(from.shot.num,1))), ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0.46,0.37,0.14,0.12]);                       
	                      
	 d =from.shot.date;
	 if isempty(d)
	 	d=zeros(1,6);
	 end
	 dds=datestr(datenum(d(1),d(2),d(3),d(4),d(5),d(6)),0);                     
	 uicontrol('style','text','tag','cz_date_text','string',dds, ...
	                      'units','normalized','interruptible','on', ...
	                      'busyaction','queue','position',[0.605,0.37,0.395,0.12]); 
	                      
end
	                     
% figure plot en fonction du temps
if ~ishandle(info.figure2)|isempty(info.figure2)
	cree = 1;
elseif ~strcmp(get(info.figure2,'tag'),'cz_plot_temps') 
	cree = 1;
else
	cree = 0;
end

if cree == 1	
	h = findobj(0,'type','figure','tag','cz_plot_temps');
	if ~isempty(h)
		delete(h);
	end
	
	info.figure2 = figure('name','Cronos-Zineb plot temps','tag','cz_plot_temps','number','off', ...
	                      'units','normalized','menubar','figure','interruptible','on', ...
	                      'busyaction','queue');
	                      
	drawnow;                      
	set(info.figure2,'position',[0.01,0.2,0.5,0.7]);                      
	
	tms = 3;
	y0 = 0.1;
	dy = 0.2;
	x0 = 0.1;
	dx = 0.85;
	xlim = [min(data.gene.temps),max(data.gene.temps)];
	t=data.gene.temps;
	if isempty(xlim)
		xlim=[0,1];
	elseif xlim(1) == xlim(2)
		xlim(2)=xlim(1) +1;
	end
	axes('units','normalized','position',[x0    y0    dx    dy],'tag','qbetali','xlim',xlim,'fontsize',12);
	xlabel('temps (s)'); 
	line(t,rand(size(t)),'color',[1 0 0],'linestyle','-','tag','q0');
	line(t,rand(size(t)),'color',[1 0 1],'linestyle','-','tag','qa');
	line(t,rand(size(t)),'color',[0 1 0],'linestyle','-','tag','ip');
	line(t,rand(size(t)),'color',[0 1 1],'linestyle','-','tag','j0sip');
	line(t,rand(size(t)),'color',[0 0.75 0.5],'linestyle','-','tag','li');
	line(t,rand(size(t)),'color',[0.75 0 0.5],'linestyle','-','tag','betadia');
	line(t,rand(size(t)),'color',[0 0 0],'linestyle','-','tag','vloop');
        if strcmp(from.machine,'TS')
	       legend('q_0','q_a','Ip/10^6','j_0/Ip','li','beta','Vloop',-1);
        else
	       legend('q_0','q_{95}','Ip/10^6','j_0/Ip','li','beta','Vloop',-1);
        end
	line(t,rand(size(t)),'color',[1 0 0],'linestyle','none','marker','o','markersize',tms,'tag','q0cons');
	line(t,rand(size(t)),'color',[0 1 0],'linestyle','none','marker','o','markersize',tms,'tag','ipcons');
	line(t,rand(size(t)),'color',[0 0.75 0.5],'linestyle','none','marker','o','markersize',tms,'tag','licons');
	line(t,rand(size(t)),'color',[0 0 0],'linestyle','none','marker','o','markersize',tms,'tag','vloopcons');
	line(t,rand(size(t)),'color',[1 0 0],'linestyle',':','tag','cq0');
	line(t,rand(size(t)),'color',[1 0 1],'linestyle',':','tag','cqa');
	line(t,rand(size(t)),'color',[0 1 0],'linestyle',':','tag','cip');
	line(t,rand(size(t)),'color',[0 1 1],'linestyle',':','tag','cj0sip');
	line(t,rand(size(t)),'color',[0 0.75 0.5],'linestyle',':','tag','cli');
	line(t,rand(size(t)),'color',[0.75 0 0.5],'linestyle',':','tag','cbetadia');
	line(t,rand(size(t)),'color',[0 0 0],'linestyle',':','tag','cvloop');
	
	y0 = y0 + dy ;
	axes('units','normalized','position',[x0    y0    dx    dy],'tag','neni','xlim',xlim,'xticklabel',[],'fontsize',12);
	line(t,rand(size(t)),'color',[1 0 0],'linestyle','-','tag','ne0');
	line(t,rand(size(t)),'color',[1 0 1],'linestyle','-','tag','nea');
	line(t,rand(size(t)),'color',[0 1 0],'linestyle','-','tag','ni0');
	line(t,rand(size(t)),'color',[0 1 1],'linestyle','-','tag','nia');
	legend('Ne_0   ','Ne_a','Ni_0','Ni_a',-1);
	line(t,rand(size(t)),'color',[1 0 1],'linestyle','none','marker','o','markersize',tms,'tag','neacons');
	line(t,rand(size(t)),'color',[1 0 0],'linestyle',':','tag','cne0');
	line(t,rand(size(t)),'color',[1 0 1],'linestyle',':','tag','cnea');
	line(t,rand(size(t)),'color',[0 1 0],'linestyle',':','tag','cni0');
	line(t,rand(size(t)),'color',[0 1 1],'linestyle',':','tag','cnia');
	ylabel('10^{19} m^{-3}');
	
	y0 = y0 + dy ;
	axes('units','normalized','position',[x0    y0    dx    dy],'tag','teti','xlim',xlim,'xticklabel',[],'fontsize',12);
	line(t,rand(size(t)),'color',[1 0 0],'linestyle','-','tag','te0');
	line(t,rand(size(t)),'color',[1 0 1],'linestyle','-','tag','tea');
	line(t,rand(size(t)),'color',[0 1 0],'linestyle','-','tag','ti0');
	line(t,rand(size(t)),'color',[0 1 1],'linestyle','-','tag','tia');
	legend('Te_0  ','10 Te_a','Ti_0','10 Ti_a',-1);
	line(t,rand(size(t)),'color',[1 0 0],'linestyle','none','marker','o','markersize',tms,'tag','teacons');
	line(t,rand(size(t)),'color',[0 1 0],'linestyle','none','marker','o','markersize',tms,'tag','tiacons');
	line(t,rand(size(t)),'color',[1 0 0],'linestyle',':','tag','cte0');
	line(t,rand(size(t)),'color',[1 0 1],'linestyle',':','tag','ctea');
	line(t,rand(size(t)),'color',[0 1 0],'linestyle',':','tag','cti0');
	line(t,rand(size(t)),'color',[0 1 1],'linestyle',':','tag','ctia');
	ylabel('keV');
	
	
	y0 = y0 + dy ;
	axes('units','normalized','position',[x0    y0    dx    dy],'tag','pepi', ...
	     'xlim',xlim,'xticklabel',[],'tag','axes_avec_titre','fontsize',12);
   title(' ''-'' = simulation, ''o'' = consigne, ''..'' = exp');
	line(t,rand(size(t)),'color',[1 0 0],'linestyle','-','tag','pe0');
	line(t,rand(size(t)),'color',[1 0 1],'linestyle','-','tag','pea');
	line(t,rand(size(t)),'color',[0 1 0],'linestyle','-','tag','pi0');
	line(t,rand(size(t)),'color',[0 1 1],'linestyle','-','tag','pia');
	legend('Pe_0   ','Pe_a','Pion_0','Pion_a',-1);
	line(t,rand(size(t)),'color',[1 0 0],'linestyle',':','tag','cpe0');
	line(t,rand(size(t)),'color',[1 0 1],'linestyle',':','tag','cpea');
	line(t,rand(size(t)),'color',[0 1 0],'linestyle',':','tag','cpi0');
	line(t,rand(size(t)),'color',[0 1 1],'linestyle',':','tag','cpia');
	ylabel('kPa');
	
end


% figure plot des profils
if ~ishandle(info.figure3)|isempty(info.figure3)
	cree = 1;
elseif ~strcmp(get(info.figure3,'tag'),'cz_plot_profil') 
	cree = 1;
else
	cree = 0;
end

if cree == 1	
	h = findobj(0,'type','figure','tag','cz_plot_profil');
	if ~isempty(h)
		delete(h);
	end
	
	info.figure3 = figure('name','Cronos-Zineb profils','tag','cz_plot_profil','number','off', ...
	                      'units','normalized','menubar','figure','interruptible','on', ...
	                      'busyaction','queue');
	                      
	drawnow;
	                      
	set(info.figure3,'position',[0.53,0.2,0.45,0.7]);                      
	
	h=subplot(3,2,1);
	set(h,'tag','pression','xlim',[0,1]);
	ylabel('Pa');
	title('Pe (r), Pion (c) et Ptot (b)')
	drawnow
	info.order=get(h,'colororder');
	
	h=subplot(3,2,2);
	set(h,'tag','electro','xlim',[0,1]);
	ylabel('su');
	title('Jmoy (MA,r), q (su, c) et E// *2 Pi R0 (V, b) ')
	
	h=subplot(3,2,3);
	set(h,'tag','temperature','xlim',[0,1]);
	ylabel('keV');
	title('Te (r) et Ti (c) ')
	
	h=subplot(3,2,5);
	set(h,'tag','densite','xlim',[0,1]);
	ylabel('10^{19} m^{-3} et su');
	xlabel('Rho/Rho_{max}')
	title('Ne (r), Ni (c) et  ae (b), Zeff (c)')
	
	h=subplot(3,2,4);
	set(h,'tag','equilibre','xlim',[0,1]);
	ylabel('keV');
	title('Psi (r), Phi (c) [- = diff, o =equi]')
	xlabel('Rho/Rho_max')
	
	h=subplot(3,2,6);
	set(h,'tag','isoflux');
	ylabel('m');
	xlabel('m');
        axis('square');	
	
end
  
% ici commence le plot
% 1 - en fonction du temps
indt = gene.kmin:(gene.k+dk);
tt=data.gene.temps(indt);
t=data.gene.temps;
ind95 =  min( iround(gene.x,0.95));
if strcmp(from.machine,'TS')
   ind95 =length(gene.x);
end

figure(info.figure2);

hl= findobj(info.figure2,'type','line','tag','q0');
set(hl,'xdata',tt,'ydata',data.prof.q(indt,1));
hl= findobj(info.figure2,'type','line','tag','cq0');
set(hl,'xdata',t,'ydata',data.exp.q0);
hl= findobj(info.figure2,'type','line','tag','qa');
set(hl,'xdata',tt,'ydata',data.prof.q(indt,ind95));
hl= findobj(info.figure2,'type','line','tag','cqa');
set(hl,'xdata',t,'ydata',data.exp.qa);
hl= findobj(info.figure2,'type','line','tag','ip');
set(hl,'xdata',tt,'ydata',data.gene.ip(indt)./1e6);
hl= findobj(info.figure2,'type','line','tag','cip');
set(hl,'xdata',t,'ydata',data.exp.ip./1e6);
hl= findobj(info.figure2,'type','line','tag','j0sip');
set(hl,'xdata',tt,'ydata',data.prof.jmoy(indt,1)./data.gene.ip(indt));
hl= findobj(info.figure2,'type','line','tag','cj0sip');
set(hl,'xdata',t,'ydata',data.exp.j0./data.exp.ip);
hl= findobj(info.figure2,'type','line','tag','li');
set(hl,'xdata',tt,'ydata',data.gene.li(indt));
hl= findobj(info.figure2,'type','line','tag','cli');
set(hl,'xdata',t,'ydata',data.exp.li);
hl= findobj(info.figure2,'type','line','tag','betadia');
set(hl,'xdata',tt,'ydata',data.gene.betadia(indt));
hl= findobj(info.figure2,'type','line','tag','cbetadia');
set(hl,'xdata',t,'ydata',data.exp.betadia);
hl= findobj(info.figure2,'type','line','tag','vloop');
set(hl,'xdata',tt,'ydata',data.gene.vloop(indt));
hl= findobj(info.figure2,'type','line','tag','cvloop');
set(hl,'xdata',t,'ydata',data.exp.vloop);

hl= findobj(info.figure2,'type','line','tag','ne0');
set(hl,'xdata',tt,'ydata',data.prof.ne(indt,1)./1e19);
hl= findobj(info.figure2,'type','line','tag','cne0');
set(hl,'xdata',t,'ydata',data.exp.ne0./1e19);
hl= findobj(info.figure2,'type','line','tag','nea');
set(hl,'xdata',tt,'ydata',data.prof.ne(indt,end)./1e19);
hl= findobj(info.figure2,'type','line','tag','cnea');
set(hl,'xdata',t,'ydata',data.exp.nea./1e19);
hl= findobj(info.figure2,'type','line','tag','ni0');
set(hl,'xdata',tt,'ydata',data.prof.ni(indt,1)./1e19);
hl= findobj(info.figure2,'type','line','tag','cni0');
set(hl,'xdata',t,'ydata',data.exp.ni0./1e19);
hl= findobj(info.figure2,'type','line','tag','nia');
set(hl,'xdata',tt,'ydata',data.prof.ni(indt,end)./1e19);
hl= findobj(info.figure2,'type','line','tag','cnia');
set(hl,'xdata',t,'ydata',data.exp.nia./1e19);

hl= findobj(info.figure2,'type','line','tag','te0');
set(hl,'xdata',tt,'ydata',data.prof.te(indt,1)./1e3);
hl= findobj(info.figure2,'type','line','tag','cte0');
set(hl,'xdata',t,'ydata',data.exp.te0./1e3);
hl= findobj(info.figure2,'type','line','tag','tea');
set(hl,'xdata',tt,'ydata',data.prof.te(indt,end)./1e2);
hl= findobj(info.figure2,'type','line','tag','ctea');
set(hl,'xdata',t,'ydata',data.exp.tea./1e2);
hl= findobj(info.figure2,'type','line','tag','ti0');
set(hl,'xdata',tt,'ydata',data.prof.ti(indt,1)./1e3);
hl= findobj(info.figure2,'type','line','tag','cti0');
set(hl,'xdata',t,'ydata',data.exp.ti0./1e3);
hl= findobj(info.figure2,'type','line','tag','tia');
set(hl,'xdata',tt,'ydata',data.prof.ti(indt,end)./1e2);
hl= findobj(info.figure2,'type','line','tag','ctia');
set(hl,'xdata',t,'ydata',data.exp.tia./1e2);

hl= findobj(info.figure2,'type','line','tag','pe0');
set(hl,'xdata',tt,'ydata',data.prof.pe(indt,1)./1e3);
hl= findobj(info.figure2,'type','line','tag','cpe0');
set(hl,'xdata',t,'ydata',data.exp.te0.*data.exp.ne0 .* phys.e./1e3);
hl= findobj(info.figure2,'type','line','tag','pea');
set(hl,'xdata',tt,'ydata',data.prof.pe(indt,end)./1e3);
hl= findobj(info.figure2,'type','line','tag','cpea');
set(hl,'xdata',t,'ydata',data.exp.tea .* data.exp.nea .* phys.e./1e3 );
hl= findobj(info.figure2,'type','line','tag','pi0');
set(hl,'xdata',tt,'ydata',data.prof.pion(indt,1)./1e3);
hl= findobj(info.figure2,'type','line','tag','cpi0');
set(hl,'xdata',t,'ydata',data.exp.ti0.*data.exp.ni0 .* phys.e./1e3 );
hl= findobj(info.figure2,'type','line','tag','pia');
set(hl,'xdata',tt,'ydata',data.prof.pion(indt,end)./1e3);
hl= findobj(info.figure2,'type','line','tag','cpia');
set(hl,'xdata',t,'ydata',data.exp.tia.*data.exp.nia .* phys.e./1e3 );

% les consignes generales
indc = find(data.mode.cons.psi == 0);
hl= findobj(info.figure2,'type','line','tag','ipcons');
if ~isempty(indc)
	set(hl,'xdata',t(indc),'ydata',data.cons.ip(indc)./1e6);
else
	set(hl,'xdata',t,'ydata',NaN.*t);
end	
indc = find(data.mode.cons.psi == 1);
hl= findobj(info.figure2,'type','line','tag','vloopcons');
if ~isempty(indc)
	set(hl,'xdata',t(indc),'ydata',data.cons.vloop(indc));
else
	set(hl,'xdata',t,'ydata',NaN.*t);
end	
indc = find(data.mode.cons.ne == 0);
hl= findobj(info.figure2,'type','line','tag','neacons');
if ~isempty(indc)
	set(hl,'xdata',t(indc),'ydata',data.cons.ne1(indc)./1e19);
else
	set(hl,'xdata',t,'ydata',NaN.*t);
end	
indc = find(data.mode.cons.pe == 0);
hl= findobj(info.figure2,'type','line','tag','teacons');
if ~isempty(indc)
	set(hl,'xdata',t(indc),'ydata',data.cons.te1(indc)./1e3);
else
	set(hl,'xdata',t,'ydata',NaN.*t);
end	
indc = find(data.mode.cons.pion == 0);
hl= findobj(info.figure2,'type','line','tag','tiacons');
if ~isempty(indc)
	set(hl,'xdata',t(indc),'ydata',data.cons.ti1(indc)./1e3);
else
	set(hl,'xdata',t,'ydata',NaN.*t);
end	

indc = find(data.mode.asser ~= 0);
hlq= findobj(info.figure2,'type','line','tag','q0cons');
hlli= findobj(info.figure2,'type','line','tag','licons');
if ~isempty(indc)
	set(hlq,'xdata',t(indc),'ydata',data.cons.asser.q0(indc));
	set(hlli,'xdata',t(indc),'ydata',data.cons.asser.li(indc));
else
	set(hlq,'xdata',t,'ydata',NaN.*t);
	set(hlli,'xdata',t,'ydata',NaN.*t);
end	

% le titre
htli = findobj(0,'tag','axes_avec_titre');
if ~isempty(htli)
	   hi = get(htli,'title');
	   strt= sprintf('(temps de %g a %g s)',gene.t,gene.t+gene.dt);  
		set(hi,'string',[' ''-'' = simulation, ''o'' = consigne, ''..'' = exp  ' ,strt]);
end

% les consignes des les asservissements


% 2 - les profils
indk = gene.k:1:(gene.k+dk); 
x    = gene.x;
figure(info.figure3);

h = findobj(info.figure3,'type','axes','tag','pression');
axes(h);
set(h,'colororder',info.order,'nextplot','replace');
plot(x,data.prof.pe(indk,:)./1e3,'r',x,data.prof.pion(indk,:)./1e3,'c',x,data.prof.ptot(indk,:)./1e3,'b', ...
     x,data.equi.ptot(indk,:)./1e3,'ob');
ylabel('kPa');
title('Pe (r), Pion (c) et Ptot (b)')
set(h,'tag','pression');

h = findobj(info.figure3,'type','axes','tag','electro');
axes(h);
set(h,'colororder',info.order,'nextplot','replace');
plot(x,data.prof.jmoy(indk,:)./1e6,'r',x,data.prof.q(indk,:),'c', ...
     x,data.prof.epar(indk,:) .* 2 .* pi .* (data.geo.r0(indk)*ones(1,size(x,2))),'b', ...
     x,data.equi.jmoy(indk,:)./1e6,'ro',x,data.equi.q(indk,:),'co',[0 1],[1 1],'g');
ylabel('su');
title('Jmoy (MA,r), q (su, c) et E// * 2 Pi R0 (V, b) ')
set(h,'tag','electro');

h = findobj(info.figure3,'type','axes','tag','temperature');
axes(h);
set(h,'colororder',info.order,'nextplot','replace');
plot(x,data.prof.te(indk,:)./1e3,'r',x,data.prof.ti(indk,:)./1e3,'c');
ylabel('keV');
title('Te (r) et Ti (c) ')
set(h,'tag','temperature');

h = findobj(info.figure3,'type','axes','tag','densite');
axes(h);
set(h,'colororder',info.order,'nextplot','replace');
plot(x,data.prof.ne(indk,:)./1e19,'r',x,data.prof.ni(indk,:)./1e19,'c', ...
     x,data.prof.ae(indk,:),'b',x,data.prof.zeff(indk,:),'m');
ylabel('10^{19} m^{-3} et su');
xlabel('Rho/Rho_{max}')
title('Ne (r), Ni (c) et  ae (b), Zeff (c)')
set(h,'tag','densite');

h = findobj(info.figure3,'type','axes','tag','equilibre');
axes(h);
set(h,'colororder',info.order,'nextplot','replace');
plot(x,data.prof.psi(indk,:),'r',x,data.equi.phi(indk,:)/2/pi,'co', ...
     x,data.equi.psi(indk,:),'or');
ylabel('V*s');
title('Psi (r), Phi/2/pi (c) [- = diff, o =equi]')
xlabel('Rho/Rho_max')
set(h,'tag','equilibre');

h = findobj(info.figure3,'type','axes','tag','isoflux');
pas =fix(size(data.equi.rhoRZ,2)./21);
if pas ==0
	pas=1;
end
indr = 1:pas:size(data.equi.rhoRZ,2);

ptot  = interp1(data.equi.rhomax(max(indk)).*gene.x,data.prof.ptot(max(indk),:),double(data.equi.rhoRZ(max(indk),:)),'spline')';
jmoy  = interp1(data.equi.rhomax(max(indk)).*gene.x,data.prof.jmoy(max(indk),:),double(data.equi.rhoRZ(max(indk),:)),'spline')';
ne    = interp1(data.equi.rhomax(max(indk)).*gene.x,data.prof.ne(max(indk),:)  ,double(data.equi.rhoRZ(max(indk),:)),'spline')';
ptot  = abs(ptot(indr));
jmoy  = abs(jmoy(indr));
ne    = abs(ne(indr));

ii     =  max(indk);
jmax   =  data.gene.ip(ii)./data.gene.surface(ii);
nemax  =  data.gene.ip(ii) ./ 1e6 ./pi ./ data.equi.rhomax(ii) .^ 2 ./ exp(1) .* 1e20;
pmax   =  nemax .* 1e4 .* phys.e;

rouge = min(ones(size(ptot)),tanh(ptot./pmax));
vert  = 1 - min(ones(size(ne)),tanh(ne./nemax));
bleu  = min(ones(size(jmoy)),tanh(jmoy./jmax));

corder = [rouge,vert,bleu];

axes(h);
set(h,'nextplot','replace');
for k =1:length(indr)
	plot(double(squeeze(data.equi.R(min(indk),indr(k),:)))',double(squeeze(data.equi.Z(min(indk),indr(k),:)))','linestyle',':','color',corder(k,:));
	hold on
	plot(double(squeeze(data.equi.R(max(indk),indr(k),:)))',double(squeeze(data.equi.Z(max(indk),indr(k),:)))','linestyle','-','color',corder(k,:));
end
%plot(double(data.geo.R(max(indk),:)),double(data.geo.Z(max(indk),:))-data.geo.z0(max(indk)),'r');
plot(double(data.geo.R(max(indk),:)),double(data.geo.Z(max(indk),:)),'r');

hold off
ylabel('m');
xlabel('m');
axis('equal');	
set(h,'tag','isoflux');

drawnow;
% fin de la fonction
