function zgclip(action)

% gestion des entrees
if nargin < 1
   action ='init';
elseif isempty(action)
   action ='init';
end

% constantes
kmax  = 5;
h = findobj(0,'type','figure','tag','clip_t');


switch action

case 'init'
    if isempty(h)
       h = figure('name','time viewgraph','tag','clip_t','color',[1 1 1]);
    else
       figure(h);
       clf
    end
    for k = 1:kmax
       ha = subplot(kmax,1,k);
       set(ha,'tag',int2str(k),'units','pixels');
       pos = get(ha,'position');
       hc = uicontrol(h,'style','checkbox', ...
                   'units','pixels', ...
                   'position',[0,pos(2)+pos(4)/2,18,18], ...
		   'callback',sprintf('zgclipt(''%d'');',k), ...
		   'tooltip','Copy of window current axes "dataplot"', ...
                   'tag',int2str(k));
       set(ha,'units','normalized');
       set(hc,'units','normalized');
    end
    uicontrol(h,'style','push','units','normalized','tag','xmin', ...
                'position',[0,0,0.1,0.05],'callback','zgclipt(''xmin'');', ...
		'tooltip','time minvalue choice', ...
		'string','Xmin');
    uicontrol(h,'style','push','units','normalized','tag','xmax', ...
                'position',[0.9,0,0.1,0.05],'callback','zgclipt(''xmax'');', ...
		'tooltip','time maxvalue choice', ...
		'string','Xmax');
    uicontrol(h,'style','push','units','normalized','tag','auto', ...
                'position',[0.45,0,0.1,0.05],'callback','zgclipt(''auto'');', ...
		'tooltip','All the avialable time', ...
		'string','Auto');
case {'1','2','3','4','5'}
   nb  = str2num(action);
   hdp = findobj(0,'type','figure','tag','zdataplot_metis');
   if isempty(hdp)
        hdp = findobj(0,'type','figure','tag','zdataplot_dc');
   end
   if isempty(hdp)
        hdp = findobj(0,'type','figure','tag','zdataplot');
   end
   if isempty(hdp)
       hfl = findobj(0,'type','figure');
       % due to sort of handle, the second one should be the last view
       % figure before click in clipboard figure.
       if length(hfl) > 1
           hdp = hfl(2);
       end
   end
   if isempty(hdp)
       return
   end
   hdep = get(hdp,'CurrentAxes');
   if isempty(hdep)
       return
   end
   hcp    = findobj(h,'type','axes','tag',action);
   posmem = get(hcp,'position');
   umem   =  get(hcp,'units');
   delete(hcp);
   hnew   = copyobj(hdep,h);
   set(hnew,'units',umem,'position',posmem,'tag',action);
   if nb > 1
      tt = get(get(hnew,'title'),'string');
      set(get(hnew,'title'),'string','');
      h1 = findobj(h,'type','axes','tag','1');
      set(get(h1,'title'),'string',tt);
   end
   if nb < kmax
      tt = get(get(hnew,'xlabel'),'string');
      set(get(hnew,'xlabel'),'string','');
      h1 = findobj(h,'type','axes','tag',int2str(kmax));
      set(get(h1,'xlabel'),'string',tt);
   end
   set(findobj(h,'type','uicontrol','tag',action),'value',0);
case 'xmin'
    [x,y] = ginput(1);
    hh    = findobj(h,'type','axes');
    for k =1:length(hh)
        xlim = get(hh(k),'xlim');
	if x < xlim(2)
	    xlim(1) = x;
	end
	set(hh(k),'xlim',xlim);
    end
case 'xmax'
    [x,y] = ginput(1);
    hh    = findobj(h,'type','axes');
    for k =1:length(hh)
        xlim = get(hh(k),'xlim');
	if x > xlim(1)
	    xlim(2) = x;
	end
	set(hh(k),'xlim',xlim);
    end
case 'auto'
   hh    = findobj(h,'type','axes');
   for k =1:length(hh)
      set(hh(k),'xlimmode','auto','ylimmode','auto');
   end

otherwise
      disp('Action not taken into account')
end
