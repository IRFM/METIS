function zgclipp(action)

% gestion des entrees
if nargin < 1
   action ='init';
elseif isempty(action)
   action ='init';
end
langue = getappdata(0,'langue_cronos');
% constantes
kmax  = 3;
h = findobj(0,'type','figure','tag','clip_p');


switch action

case 'init'
    if isempty(h)
       if strcmp(langue,'francais')
         h = figure('name','Fenetre des profils','tag','clip_p','color',[1 1 1]);
       end
       if strcmp(langue,'anglais')
         h = figure('name','profile viewgraph','tag','clip_p','color',[1 1 1]);
       end
    else
       figure(h);
       clf
    end
    for k = 1:(kmax*2)
       ha = subplot(kmax,2,k);
       set(ha,'tag',int2str(k),'units','pixels');
       pos = get(ha,'position');
       if rem(k,2) == 0
          xpos = pos(1) + pos(3) +20;
       else 
          xpos = 0;
       end
       hc = uicontrol(h,'style','checkbox', ...
                   'units','pixels', ...
                   'position',[xpos,pos(2)+pos(4)/2,18,18], ...
		   'callback',sprintf('zgclipp(''%d'');',k), ...
 		   'tooltip','Copy of window current axes "dataplot"', ...
                  'tag',int2str(k));
       set(ha,'units','normalized');
       set(hc,'units','normalized');
    end
    zplotprof(findobj(h,'type','axes','tag','1'),[0,1],[0,1],NaN.*ones(2,2));
case {'1','2','3','4','5','6'}
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
   set(findobj(h,'type','uicontrol','tag',action),'value',0);
   
  
otherwise
   if strcmp(langue,'francais')
      disp('Action non prise en compte')
   end
   if strcmp(langue,'anglais')
      disp('Action not taken into account')
   end
end
