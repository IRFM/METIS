function rapfile=z0drapport_png(rapfile,nointer,temps_visu)

if nargin < 2
  nointer = 0;
end
if nointer ~= 0
  ws = 'caller';
else
  ws = 'base';
end

info = evalin(ws,'post.z0dinput');
if nargin ==0
	rapfile = sprintf('z0d_%s_%d',info.machine,fix(info.shot));
elseif isempty(rapfile)
	rapfile = sprintf('z0d_%s_%d',info.machine,fix(info.shot));
end

if nargin < 3
  % demande du temps initial et final
  % choix de la langue
  langue =  lower(getappdata(0,'langue_cronos'));
  temps  = evalin(ws,'post.zerod.temps');
  tdeb = num2str(min(temps));
  tfin = num2str(max(temps));
  switch langue
  case 'francais'
	prompt={'temps de debut (s):','temps de fin (s):'};
	def={tdeb,tfin};
	dlgTitle='Base temps pour le 0D';
  otherwise
	prompt={'begin time (s):','end time (s):'};
	def={tdeb,tfin};
	dlgTitle='0D time slices';
  end
  lineNo=1;
  answer=inputdlg(prompt,dlgTitle,lineNo,def);
  if isempty(answer)
	  return
  end
  tdeb  = str2num(answer{1});
  tfin  = str2num(answer{2});
else
  temps  = evalin(ws,'post.zerod.temps');
  tdeb = min(temps);
  tfin = max(temps);  
end

if nargin < 3
  temps_visu = (tdeb + tfin) ./ 2;
end

warning off

etatmem = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');
liste_of_figure = findobj(0,'type','figure');
set(0,'ShowHiddenHandles',etatmem);

evalin(ws,'z0plot_sepa_evolution');
saveas(gcf,sprintf('%s_sepa_evolution',rapfile),'fig');
saveas(gcf,sprintf('%s_sepa_evolution',rapfile),'png');
close(gcf);

evalin(ws,'z0plotsc');
tlim(gcf,tdeb,tfin);
fullscreen
saveas(gcf,sprintf('%s_scenario',rapfile),'fig');
saveas(gcf,sprintf('%s_scenario',rapfile),'png');
close(gcf);

evalin(ws,'z0plotp');
tlim(gcf,tdeb,tfin);
fullscreen(1)
saveas(gcf,sprintf('%s_power',rapfile),'fig');
saveas(gcf,sprintf('%s_power',rapfile),'png');
close(gcf);

evalin(ws,'z0plote');
tlim(gcf,tdeb,tfin);
fullscreen(1)
saveas(gcf,sprintf('%s_energie',rapfile),'fig');
saveas(gcf,sprintf('%s_energie',rapfile),'png');
close(gcf);

evalin(ws,'z0plotj');
tlim(gcf,tdeb,tfin);
fullscreen
saveas(gcf,sprintf('%s_current',rapfile),'fig');
saveas(gcf,sprintf('%s_current',rapfile),'png');
close(gcf);

evalin(ws,'z0plotn');
tlim(gcf,tdeb,tfin);
fullscreen(1)
saveas(gcf,sprintf('%s_density',rapfile),'fig');
saveas(gcf,sprintf('%s_density',rapfile),'png');
close(gcf);

evalin(ws,'z0plott');
tlim(gcf,tdeb,tfin);
fullscreen
saveas(gcf,sprintf('%s_temperature',rapfile),'fig');
saveas(gcf,sprintf('%s_temperature',rapfile),'png');
close(gcf);

evalin(ws,'z0plotc');
tlim(gcf,tdeb,tfin);
fullscreen(1)
saveas(gcf,sprintf('%s_confinement',rapfile),'fig');
saveas(gcf,sprintf('%s_confinement',rapfile),'png');
close(gcf);

evalin(ws,'z0ploteq');
tlim(gcf,tdeb,tfin);
fullscreen
saveas(gcf,sprintf('%s_equilibrium',rapfile),'fig');
saveas(gcf,sprintf('%s_equilibrium',rapfile),'png');
close(gcf);

evalin(ws,'z0plotlh');
tlim(gcf,tdeb,tfin);
fullscreen
saveas(gcf,sprintf('%s_LHCD',rapfile),'fig');
saveas(gcf,sprintf('%s_LHCD',rapfile),'png');
close(gcf);

evalin(ws,'z0plotgeo');
tlim(gcf,tdeb,tfin);
fullscreen
saveas(gcf,sprintf('%s_geometry',rapfile),'fig');
saveas(gcf,sprintf('%s_geometry',rapfile),'png');
close(gcf);

evalin(ws,'z0plotconv');
tlim(gcf,tdeb,tfin);
fullscreen
saveas(gcf,sprintf('%s_convergence',rapfile),'fig');
saveas(gcf,sprintf('%s_convergence',rapfile),'png');
close(gcf);

evalin(ws,'z0plotshine');
tlim(gcf,tdeb,tfin);
fullscreen
saveas(gcf,sprintf('%s_shinethrought',rapfile),'fig');
saveas(gcf,sprintf('%s_shinethrought',rapfile),'png');
close(gcf);

evalin(ws,'z0plotsoq');
tlim(gcf,tdeb,tfin);
fullscreen
saveas(gcf,sprintf('%s_soq',rapfile),'fig');
saveas(gcf,sprintf('%s_soq',rapfile),'png');
close(gcf);

evalin(ws,'z0plotcoherence_0D_1D');
tlim(gcf,tdeb,tfin);
fullscreen
saveas(gcf,sprintf('%s_consistency_0d_1d',rapfile),'fig');
saveas(gcf,sprintf('%s_consistency_0d_1d',rapfile),'png');
close(gcf);

evalin(ws,'z0plot2points');
warning off
drawnow
if ishandle(findobj(0,'type','figure','tag','z0W'))
  tlim(gcf,tdeb,tfin);
  fullscreen;
  saveas(gcf,sprintf('%s_tungsten_effects',rapfile),'fig');
  saveas(gcf,sprintf('%s_tungsten_effects',rapfile),'png');
  close(gcf);
end
tlim(gcf,tdeb,tfin);
fullscreen(1);
saveas(gcf,sprintf('%s_screening',rapfile),'fig');
saveas(gcf,sprintf('%s_screening',rapfile),'png');
close(gcf);
tlim(gcf,tdeb,tfin);
fullscreen;
saveas(gcf,sprintf('%s_2points',rapfile),'fig');
saveas(gcf,sprintf('%s_2points',rapfile),'png');
close(gcf);

evalin(ws,'z0plotmode');
fullscreen
saveas(gcf,sprintf('%s_KBM',rapfile),'fig');
saveas(gcf,sprintf('%s_KBM',rapfile),'png');
close(gcf);
fullscreen
saveas(gcf,sprintf('%s_critical_gradient',rapfile),'fig');
saveas(gcf,sprintf('%s_critical_gradient',rapfile),'png');
close(gcf);
drawnow
hpf_curs =  findobj(gcf,'type','uicontrol','tag','curseur');
if ~isempty(hpf_curs)
	set(hpf_curs,'value',temps_visu);
	zplotprof('curseur',gcf);
end
drawnow
fullscreen
saveas(gcf,sprintf('%s_TEM_ITG',rapfile),'fig');
saveas(gcf,sprintf('%s_TEM_ITG',rapfile),'png');
close(gcf);

evalin(ws,'z0plotstationnary');
tlim(gcf,tdeb,tfin);
fullscreen(1)
saveas(gcf,sprintf('%s_stationnary',rapfile),'fig');
saveas(gcf,sprintf('%s_stationnary',rapfile),'png');
close(gcf);

evalin(ws,'z0plotdivertor');
drawnow
hpf_curs =  findobj(gcf,'type','uicontrol','tag','curseur');
if ~isempty(hpf_curs)
	set(hpf_curs,'value',temps_visu);
	zplotprof('curseur',gcf);
end
drawnow
fullscreen
saveas(gcf,sprintf('%s_divertor',rapfile),'fig');
saveas(gcf,sprintf('%s_divertor',rapfile),'png');
close(gcf);

       
evalin(ws,'z0plotl2h');
tlim(gcf,tdeb,tfin);
fullscreen
saveas(gcf,sprintf('%s_L2H',rapfile),'fig');
saveas(gcf,sprintf('%s_L2H',rapfile),'png');
close(gcf);

evalin(ws,'z0plotrad');              
tlim(gcf,tdeb,tfin);
fullscreen(1)
saveas(gcf,sprintf('%s_radiation_0d',rapfile),'fig');
saveas(gcf,sprintf('%s_radiation_0d',rapfile),'png');
close(gcf);
hpf_curs =  findobj(gcf,'type','uicontrol','tag','curseur');
if ~isempty(hpf_curs)
	set(hpf_curs,'value',temps_visu);
	zplotprof('curseur',gcf);
end
fullscreen(1)
saveas(gcf,sprintf('%s_radiation_1d',rapfile),'fig');
saveas(gcf,sprintf('%s_radiation_1d',rapfile),'png');
close(gcf);

% z0plotflux.m
evalin(ws,'z0plotflux');
fullscreen;
saveas(gcf,sprintf('%s_poloidal_flux',rapfile),'fig');
saveas(gcf,sprintf('%s_poloidal_flux',rapfile),'png');
close(gcf);

% bilan gas
NOINTER = 1;
evalin(ws,'metis_gaz');
if info.option.gaz == 3
	fullscreen;
	saveas(gcf,sprintf('%s_bilan_tritium',rapfile),'fig');
	saveas(gcf,sprintf('%s_bilan_tritium',rapfile),'png');
	close(gcf);
end
fullscreen;
saveas(gcf,sprintf('%s_bilan_gas',rapfile),'fig');
saveas(gcf,sprintf('%s_bilan_gas',rapfile),'png');
close(gcf);
             
% remove unwanted figures
etatmem = get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');
delete(setdiff(findobj(0,'type','figure'),liste_of_figure));
set(0,'ShowHiddenHandles',etatmem);
warning on
    

function tlim(hf,tdeb,tfin)

ha = findobj(hf,'type','axes');
set(ha,'xlim',[tdeb,tfin],'ylimmode','auto','yscale','linear');
drawnow
for k=1:length(ha)
	hc = ha(k);
	ylim = get(hc,'ylim');
	if ylim(1) > 0
		ylim(1) =0;
	end
	set(hc,'ylim',ylim);
end



function fullscreen(update_legende)

warning off
% delete ui objects
delete(findobj(gcf,'type','uicontrol'));

% update legend
if nargin > 0
  ha = findobj(gcf,'type','axes');
  for l=1:length(ha)
      % test if is legend
      ax = ha(l);
      if length(ax) ~= 1 || ~ishandle(ax)
	  tf=false;
      else
	  tf=isa(handle(ax),'scribe.legend');
      end
      if tf
	axes(ax);
	legend('location','best');
      end
  end
end

% set figure size
set(gcf,'Position',[1,1,1920,1080],'PaperUnits','points','PaperOrientation','landscape','PaperPosition',[1,1,1920,1080], ...
        	'PaperPositionMode','manual','PaperSize',[1922 1082],'PaperType','<custom>');
% redraw figure
edition2
drawnow

