% script de plot des temperatures
tt   =  input('temps ? ');
vece =  fix(input(' voies ECE ?'));


chargedata
indt1 = fix(interp1(t,1:length(t),tt,'nearest'));
indt2 = fix(interp1(bile.times,1:length(bile.times),tt,'nearest'));
indth =  find(sbile.te.espace == 1);
indsh =  find(sbile.te.espace == 2  | sbile.te.espace == 6)
if ~isempty(indsh)
   indsh  = indsh(vece);
end
xte   = tsplinet(ones(length(bile.times),1)*bile.rhofit,xbile,bile.rhote);

h = findobj(0,'type','figure','tag','compare');
if isempty(h)
       h=figure('tag','compare');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',9)
leg = {};
cc  = get(gca,'colororder');
for k =1:length(tt)
  plot(x,data.prof.te(indt1(k),:)/1e3,'color',cc(k,:));
  leg{end+1} = sprintf('Cronos @ %6.3g s',tt(k));
  hold on
  if ~isempty(indth)
    indl = indth(find(bile.te(indt2(k),indth)>0));
    plot(xte(indt2(k),indl),bile.te(indt2(k),indl),'linestyle','+','color',cc(k,:));
    leg{end+1} = sprintf('Thomson @ %6.3g s ',bile.times(indt2(k)));
  end
  if ~isempty(indsh)
    indl = indsh(find(bile.te(indt2(k),indsh)>0));
    plot(xte(indt2(k),indl),bile.te(indt2(k),indl),'linestyle','o','color',cc(k,:));
    leg{end+1} = sprintf('ECE @ %6.3g s',bile.times(indt2(k)));
  end
  title(sprintf('Electronic temperature, shot # %d',fix(param.from.shot.num)));
  ylabel('keV')
  xlabel('x')
  legend(leg{:},-1);
end  
  % compare mesure
  if ~isempty(bile)
    hh = findobj(0,'type','figure','tag','comp@mes_th');
    if isempty(hh)
       hh=figure('tag','comp@mes_th');
    else
       figure(hh);
    end   
    clf
    set(hh,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',9)
    xte = tsplinet(ones(length(bile.times),1)*bile.rhofit,xbile,bile.rhote);
 
    if ~isempty(indth)
      nb = length(find(any(bile.te(:,indth)>0,1)));
      if ~isempty(nb)
         ha = joint_axes(hh,nb);
         iax =1;
         for k =1:length(indth);
            indc = indth(k);
            tem  = bile.te(:,indc);
            if any( tem >0)
               indx = min(find(mean(xte(:,indc)) <= param.gene.x));
               axes(ha(iax));
               iax =iax +1;
               plot(bile.times,bile.te(:,indc),'b',data.gene.temps,data.prof.te(:,indx)./1e3,'r'); 
               if (iax-1) == nb
                  ylabel(sprintf('x = %g',param.gene.x(indx)));
               else
                  ylabel(sprintf('%g',param.gene.x(indx)));
               end
               if k == 1
                  title(sprintf('electronic temperature : Thomson (blue) / simulation (red) ; shot # %d',fix(param.from.shot.num)));
               end
            end
         end
         xlabel('time (s)');
         forme_axes(ha);
         set(ha,'xlim',data.gene.temps([param.gene.kmin,param.gene.k]));
       end
    end
    hh = findobj(0,'type','figure','tag','comp@mes_sh');
    if isempty(hh)
       hh=figure('tag','comp@mes_sh');
    else
       figure(hh);
    end   
    clf
    set(hh,'defaultaxesfontsize',16,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',9);
    if ~isempty(indsh)
      nc = length(find(any(bile.te(:,indsh)>0,1)));
      if ~isempty(nc)
         ha = joint_axes(hh,nc);
         iax = 1;
         for k =1:length(indsh)
            indc = indsh(k);
            tem  = bile.te(:,indc);
            if any( tem >0)
               indx = min(find(mean(xte(tem >0 ,indc)) <= param.gene.x));
               axes(ha(iax));
               iax =iax +1;
               plot(bile.times,bile.te(:,indc),'b',data.gene.temps,data.prof.te(:,indx)./1e3,'r'); 
               if (iax-1) == nc
                  ylabel(sprintf('x = %g',param.gene.x(indx)));
               else
                  ylabel(sprintf('%g',param.gene.x(indx)));
               end
               if k == 1
                  title(sprintf('electronic temperature : ECE (blue) / simulation (red) ; shot # %d',fix(param.from.shot.num)));
               end
            end
         end
         xlabel('time (s)');
         forme_axes(ha);
         set(ha,'xlim',data.gene.temps([param.gene.kmin,param.gene.k]));
      end
   end 
end
 
