% script de plot des temperatures

chargedata
h = findobj(0,'type','figure','tag','compare');
if isempty(h)
       h=figure('tag','compare');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])


switch param.from.machine

case 'TS'


  subplot(2,2,1)
  plotprof(gca,t,x,data.prof.pe,'color','r');
  if ~isempty(t1)
     plotprof(gca,t1,x,jeux1.data.prof.pe,'linestyle','o');
  end
  title(sprintf('Pel shot %d',fix(param.from.shot.num)))
  ylabel('Pa')

  legend('CRONOS','reference');

  subplot(2,2,2)
  plotprof(gca,t,x,data.prof.pion,'color','r');
  if ~isempty(t1)
    plotprof(gca,t1,x,jeux1.data.prof.pion,'linestyle','o');
  end

  title('Pion ')
  ylabel('Pa')
  legend('CRONOS','reference');


  subplot(2,2,3)
  plotprof(gca,t,x,data.prof.te/1e3,'color','r');
  if ~isempty(t1)
    plotprof(gca,t1,x,jeux1.data.prof.te/1e3,'linestyle','o');
  end
  title('Te ')
  ylabel('keV')
  xlabel('rho')

  legend('CRONOS','reference');

  subplot(2,2,4)
  plotprof(gca,t,x,data.prof.ti/1e3,'color','r');
  if ~isempty(t1)
    plotprof(gca,t1,x,jeux1.data.prof.ti/1e3,'linestyle','o');
  end
  title('Ti')
  ylabel('keV')
  xlabel('rho')
  legend('CRONOS','reference');

%
% trace temporelle
%
  if ~isempty(t1)
    hh = findobj(0,'type','figure','tag','compare2');
    if isempty(hh)
       hh=figure('tag','compare2');
    else
       figure(hh);
    end   
    clf
    set(hh,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1])


	 fin1    = iround(x,0.2);
    v3      = min(min(data.prof.te(:,fin1)/1e3),min(jeux1.data.prof.te(:,fin1)/1e3))*0.9;
    v4      = max(max(data.prof.te(:,fin1)/1e3),max(jeux1.data.prof.te(:,fin1)/1e3))*1.1;
    tdeb    = max(param.gene.tdeb,jeux1.param.gene.tdeb);
	 tfin    = min(param.gene.tfin,jeux1.param.gene.tfin);
    vgraph1 = [tdeb tfin v3 v4];    
    h1=subplot(4,1,1);
    plot(t,data.prof.te(:,fin1)/1e3,'b',t1,jeux1.data.prof.te(:,fin1)/1e3,'ro')
    axis(vgraph1)
    ylabel('keV')
    set(h1,'xticklabel','')
  
    fin2    = iround(x,0.4);
    v3      = min(min(data.prof.te(:,fin2)/1e3),...
                  min(jeux1.data.prof.te(:,fin2)/1e3))*0.9;
    v4      = max(max(data.prof.te(:,fin2)/1e3),max(jeux1.data.prof.te(:,fin2)/1e3))*1.1;
    vgraph2 = [tdeb tfin v3 v4];    
    h2=subplot(4,1,2);
    plot(t,data.prof.te(:,fin2)/1e3,'b',t1,jeux1.data.prof.te(:,fin2)/1e3,'ro')
    title(['Te(x = ',num2str(x(fin2),3),')'])
    axis(vgraph2)
    ylabel('keV')
    set(h2,'xticklabel','')

    fin3    = iround(x,0.6);
    v3      = min(min(data.prof.te(:,fin3)/1e3),min(jeux1.data.prof.te(:,fin3)/1e3))*0.9;
    v4      = max(max(data.prof.te(:,fin3)/1e3),max(jeux1.data.prof.te(:,fin3)/1e3))*1.1;
    vgraph3 = [tdeb tfin v3 v4];    
    h3=subplot(4,1,3);
    plot(t,data.prof.te(:,fin3)/1e3,'b',t1,jeux1.data.prof.te(:,fin3)/1e3,'ro')
    title(['Te(x = ',num2str(x(fin3),3),')'])
    axis(vgraph3)
    ylabel('keV')
    set(h3,'xticklabel','')

    fin4    = iround(x,0.8);
	 v3      = min(min(data.prof.te(:,fin4)/1e3),min(jeux1.data.prof.te(:,fin4)/1e3))*0.9;
    v4      = max(max(data.prof.te(:,fin4)/1e3),max(jeux1.data.prof.te(:,fin4)/1e3))*1.1;
    vgraph4 = [tdeb tfin v3 v4];    
    h4=subplot(4,1,4);
    plot(t,data.prof.te(:,fin4)/1e3,'b',t1,jeux1.data.prof.te(:,fin4)/1e3,'ro')
    title(['Te(x = ',num2str(x(fin4),3),')'])
    axis(vgraph4)
    ylabel('keV')
    xlabel('temps') 
  end  
%
% erreur par rapport au fit
%
  if ~isempty(t1)
    hh = findobj(0,'type','figure','tag','compare3');
    if isempty(hh)
       hh=figure('tag','compare3');
    else
       figure(hh);
    end   
    clf
    set(hh,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1])

    tdeb    = max(param.gene.tdeb,jeux1.param.gene.tdeb);
	 tfin    = min(param.gene.tfin,jeux1.param.gene.tfin);
    vgraph1 = [tdeb tfin 0 1];    
    h1      = subplot(4,1,1);
	 fin1    = iround(x,0.2);
	 nte     = tsample(data.prof.te(:,fin1),t,t1);
	 expte   = tsample(bile.tefit(:,iround(bile.rhofit,0.2)),bile.times,t1)*1000;
    plot(t1,abs(expte-jeux1.data.prof.te(:,fin1))./expte,'b',t1,abs(expte-nte)./expte,'b--')    
    axis(vgraph1)
	 grid
	 ylabel('a)')
    set(h1,'xticklabel','')
  
    vgraph2 = [tdeb tfin 0 1];    
    h2=subplot(4,1,2);
	 fin2    = iround(x,0.4);
	 nte     = tsample(data.prof.te(:,fin2),t,t1);
	 expte   = tsample(bile.tefit(:,iround(bile.rhofit,0.4)),bile.times,t1)*1000;
    plot(t1,abs(expte-jeux1.data.prof.te(:,fin2))./expte,'b',t1,abs(expte-nte)./expte,'b--')    
	 ylabel('b)')
	 grid
    axis(vgraph2)

    set(h2,'xticklabel','')
    vgraph3 = [tdeb tfin 0 1];    
    h3      = subplot(4,1,3);
	 fin3    = iround(x,0.6);
	 grid
	 nte     = tsample(data.prof.te(:,fin3),t,t1);
	 expte   = tsample(bile.tefit(:,iround(bile.rhofit,0.6)),bile.times,t1)*1000;
    plot(t1,abs(expte-jeux1.data.prof.te(:,fin3))./expte,'b',t1,abs(expte-nte)./expte,'b--')    
	 ylabel('c)')
    axis(vgraph3)
    set(h3,'xticklabel','')

    vgraph4 = [tdeb tfin 0 1];    
    h4      = subplot(4,1,4);
	 fin4    = iround(x,0.8);
	 nte     = tsample(data.prof.te(:,fin4),t,t1);
	 expte   = tsample(bile.tefit(:,iround(bile.rhofit,0.8)),bile.times,t1)*1000;
    plot(t1,abs(expte-jeux1.data.prof.te(:,fin4))./expte,'b',t1,abs(expte-nte)./expte,'b--')    
	 ylabel('d)')
	 grid
    axis(vgraph4)
    xlabel('temps') 
  end  

  set(gcf,'ResizeFcn','')
  
  % compare mesure
  if ~isempty(bile)
    hh = findobj(0,'type','figure','tag','comp@mes_th');
    if isempty(hh)
       hh=figure('tag','comp@mes_th');
    else
       figure(hh);
    end   
    clf
    set(hh,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1])
    xte = tsplinet(ones(length(bile.times),1)*bile.rhofit,xbile,bile.rhote);
    %plotprof(gca,bile.times,xte(:,1:12),bile.te(:,1:12),'linestyle','+','color','c');
    %plotprof(gca,bile.times,xte(:,13:end),bile.te(:,13:end),'linestyle','x','color','m');
  
    indth =  find(sbile.te.espace == 1);
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
               ylabel(sprintf('%g',param.gene.x(indx)));
               if k == 1
                  title(sprintf('electronic temperature : Thomson (blue) / simulation (red) ; choc # %d',param.from.shot.num));
               end
            end
         end
         xlabel('temps (s)');
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
    set(hh,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1]);
    indsh =  find(sbile.te.espace == 2  | sbile.te.espace == 6);
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
               ylabel(sprintf('%g',param.gene.x(indx)));
               if k == 1
                  title(sprintf('electronic temperature : ECE (blue) / simulation (red) ; choc # %d',param.from.shot.num));
               end
            end
         end
         xlabel('temps (s)');
         forme_axes(ha);
         set(ha,'xlim',data.gene.temps([param.gene.kmin,param.gene.k]));
      end
   end 
end
 
case 'JET'


  subplot(2,2,1)
  plotprof(gca,t,x,data.prof.pe,'color','r');
  if ~isempty(t1)
     plotprof(gca,t1,x,jeux1.data.prof.pe,'linestyle','o');
  end
  title(sprintf('Pel shot %d',fix(param.from.shot.num)))
  ylabel('Pa')
  if ~isempty(t1)
    legend('CRONOS','reference');
  else
    legend('CRONOS');
  end
  subplot(2,2,2)
  plotprof(gca,t,x,data.prof.pion,'color','r');
  title('Pion ')
  ylabel('Pa')
  legend('CRONOS');


  subplot(2,2,3)
  plotprof(gca,t,x,data.prof.te/1e3,'color','r');
  if ~isempty(t1)
    plotprof(gca,t1,x,jeux1.data.prof.te/1e3,'linestyle','o');
  end
  if ~isempty(jetteshthom)
    plotprof(gca,jetteshthom.ttex,x,jetteshthom.tex/1e3,'linestyle','--','color','b');
  else
    if ~isempty(jettethom)  
      plotprof(gca,jettethom.ttex,x,jettethom.tex/1e3,'linestyle','--','color','b');
    end  
  end
  if ~isempty(jetprof)
    if ~isempty(jetprof.terhonth)
      plotprof(gca,jetprof.trhonth,jetprof.rhonth,jetprof.terhonth/1e3,'linestyle','x','color','m');
    end
    if ~isempty(jetprof.terhonshf)
      plotprof(gca,jetprof.trhonshf,jetprof.rhonshf,jetprof.terhonshf/1e3,'linestyle','+','color','c');
    end    
  end  
  title('Te ')
  ylabel('keV')
  xlabel('rho')
  if ~isempty(t1)
    if ~isempty(jetprof)
      if ~isempty(jetprof.terhonth)
        if ~isempty(jetprof.terhonshf)
          legend('CRONOS','reference','fit','Thomson','ECE');
        else
          legend('CRONOS','reference','fit','Thomson');
        end	
      else
        legend('CRONOS','reference');
      end
    else
        legend('CRONOS','reference');        
    end      
  else
    if ~isempty(jetprof)
      if ~isempty(jetprof.terhonth)
        if ~isempty(jetprof.terhonshf)
          legend('CRONOS','fit','Thomson','ECE');
        else
          legend('CRONOS','fit','Thomson');
        end	
      else
        legend('CRONOS');
      end
    else
        legend('CRONOS');        
    end   
  end
  subplot(2,2,4)
  plotprof(gca,t,x,data.prof.ti/1e3,'color','r');
  if ~isempty(t1)
    plotprof(gca,t1,x,jeux1.data.prof.ti/1e3,'linestyle','o');
  end
  tiok = 0;
  if ~isempty(jetprof)
    if ~isempty(jetprof.ticx)
      tiok = 1;
      plotprof(gca,jetprof.tcx,jetprof.rhoncx,jetprof.ticx/1e3,'linestyle','x','color','r');
    end
  end  
  
  if ~isempty(t1) & tiok == 1
    legend('CRONOS ','reference','CX')
  elseif ~isempty(t1)
    legend('CRONOS ','reference')
  elseif tiok == 1
    legend('CRONOS ','CX')
  else
    legend('CRONOS ')
  end
  ylabel('keV')
  xlabel('rho')
  set(gcf,'ResizeFcn','')
  
otherwise
        disp('pas de donnees pour cette machine')
end


