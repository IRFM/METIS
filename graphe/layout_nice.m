function layout_nice(hf)

% 
if nargin == 0
  hf = gcf;
end
figure(hf);
fname = get(hf,'name');
% get list of true axes
hall = findobj(hf,'type','axes');
hl = findobj(hf,'type','axes','tag','legend');
ha = [];
xlabel_def = '';
ylabel_def = '';
for k=1:length(hall)
  if ~any(hall(k) == hl)
      ha(end+1) = hall(k);
      if isempty(xlabel_def)
	  xlabel_def = get(get(hall(k),'xlabel'),'string');
      end
      if isempty(ylabel_def)
	  ylabel_def = get(get(hall(k),'ylabel'),'string');
      end
  end
end


% layout
dpos=-20; pos1=256; pos2=247; pos3=364; pos4=425; 
psize1=1; psize10=1; psize2=1; psize3=1; psize4=1; psize5=1; psize6=1;
nfig=0;lls='-'; posfig=0;

% color order
co = [
         0         0    1.0000
         0    0.5000         0
    1.0000         0         0
         0    0.7500    0.7500
    0.7500         0    0.7500
    0.7500    0.7500         0
    0.2500    0.2500    0.2500
];    



% looop on axes
for k=1:length(ha)
    %axes(ha(k));
    %drawnow
    [LEGH,OBJH,OUTH,OUTM] = legend(ha(k));
    % retrieve line
    xlim = get(ha(k),'xlim');
    xlim(1) = floor(xlim(1));
    xlim(2) = ceil(xlim(2));
    ylim = get(ha(k),'ylim');
    ylim(1) = floor(ylim(1));
    ylim(2) = ceil(ylim(2));
    hl = cat(1,findobj(ha(k),'type','line'),findobj(ha(k),'type','patch'));
    % recompute limit to match axes
    xmin = Inf;
    xmax = -Inf;
    ymin = Inf;
    ymax = -Inf;
    for l=1:length(hl)
	xy = get(hl(l),'xdata');
	xmin = min(xmin,min(xy));
	xmax = max(xmax,max(xy));
	xy = get(hl(l),'ydata');
	ymin = min(ymin,min(xy));
	ymax = max(ymax,max(xy));	
    end
    if xmin < xmax
	xlim = cat(1,floor(xmin),ceil(xmax));
    end
    if ymin < ymax
	ylim = cat(1,floor(ymin),ceil(ymax));
    end  
    hui =findobj(hf,'type','uicontrol','style','edit','tag','temps');
    posfig=posfig+dpos;
    nfig=length(findobj(0,'type','figure')) + 1;
    if ~isempty(hui)
      time = get(hui,'string');
      if isempty(xlabel_def)
	xlabel_def = 'x';
      end
      if length(ha) > 1
	  figuren(nfig); set(gcf,'name',sprintf('%s (%d/%d)',fname,k,length(ha)),'Position',[pos1+posfig pos2+posfig pos3*psize1 pos4*psize1]);
      else
	  figuren(nfig); set(gcf,'name',fname,'Position',[pos1+posfig pos2+posfig pos3*psize1 pos4*psize1]);      
      end
      hloc = gca;
      set(hloc,'FontSize',fix(18*psize10),'defaultlinelinewidth',2,'defaultlinelinestyle',lls,'box','on')
      hold on
      %drawnow
      % copy lines
      hlnew = copyobj(hl,hloc);
      set(findobj(hlnew,'type','line'),'linewidth',2);
      tit=['t = ',time,' s'];
      title(tit);
      if isempty(get(get(ha(k),'xlabel'),'string'))
	  xlabel(xlabel_def);
      else
	  xlabel(get(get(ha(k),'xlabel'),'string'));
      end
      if isempty(get(get(ha(k),'ylabel'),'string'))
         ylabel(ylabel_def);
      else
         ylabel(get(get(ha(k),'ylabel'),'string'));
      end 
      legend(OUTM,'Location','best');
      set(hloc,'xlim',xlim);
      set(hloc,'ylim',ylim);
    else
      time = 'NaN';
      if isempty(xlabel_def)
	xlabel_def = 'time (s)';
      end
      if length(ha) > 1
	  figuren(nfig); set(gcf,'name',sprintf('%s (%d/%d)',fname,k,length(ha)),'Position',[pos1+posfig pos2+posfig pos3*psize3 ceil(pos4*psize3 * 0.44)]);
      else
 	  figuren(nfig); set(gcf,'name',fname,'Position',[pos1+posfig pos2+posfig pos3*psize3 ceil(pos4*psize3 * 0.44)]);     
      end
      hloc = gca;
      set(hloc,'fontsize',fix(14*psize3),'defaultlinelinewidth',2,'defaultlinelinestyle',lls,'box','on');
      % copy lines
      hlnew = copyobj(hl,hloc);
      set(findobj(hlnew,'type','line'),'linewidth',2);
      if isempty(get(get(ha(k),'xlabel'),'string'))
	  xlabel(xlabel_def);
      else
	  xlabel(get(get(ha(k),'xlabel'),'string'));
      end
      if isempty(get(get(ha(k),'ylabel'),'string'))
         ylabel(ylabel_def);
      else
         ylabel(get(get(ha(k),'ylabel'),'string'));
      end 
      legend(OUTM,'Location','best');
      set(hloc,'xlim',xlim);
      set(hloc,'ylim',ylim);
    end
    drawnow

end

% recover figure handle
function h = figuren(n)

tag = sprintf('nice_layout_%d',n);
h = findobj(0,'type','figure','tag',tag);
if isempty(h)
       h=figure('tag',tag,'color',[1 1 1]);
else
       figure(h);
       set(h,'color',[1 1 1]);
end   
clf
setcoloroldcolororder(h);

function setcoloroldcolororder(hfig)

if nargin == 0
  hfig = gcf;
end

co = [
         0         0    1.0000
         0    0.5000         0
    1.0000         0         0
         0    0.7500    0.7500
    0.7500         0    0.7500
    0.7500    0.7500         0
    0.2500    0.2500    0.2500
];    



set(hfig,'defaultAxesColorOrder',co);



function xlabelrho

if verLessThan('matlab','8.6')
  xlabel('r','fontname','symbol');
else
  xlabel('\rho');
end


function liste = switchuid(liste,shot)

switch shot
case 86614
  for k=1:length(liste)
    liste{k} = sprintf('%s?uid=cgiroud+seq=977',liste{k});
  end
end
