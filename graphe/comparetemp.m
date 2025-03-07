% comparaison des donnees de la chaleur en fonction du temps
times = [];
chargedata

switch param.from.machine

case 'TS'

   h = findobj(0,'type','figure','tag','comparetemp_1');
   if isempty(h)
       h=figure('tag','comparetemp_1');
   else
       figure(h);
   end
   clf
   set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   ipmax    = max(data.exp.ip/1e6)*1.1;
   phybmax  = max(data.gene.paddhyb/1e6)*1.1;
   pfcimax  = max(data.gene.paddfci/1e6)*1.1;
   pfcicon  = -sum(real(data.cons.fci)')/1e6;
   nbarmax  = max(data.gene.nbar/1e19)*1.1;
   wemax1   = max(bile.wetot)*1.1;
   wemax2   = max(data.gene.we/1e6)*1.1;
   rlwmax   = max(bile.rlw)*1.1;
   wdiamax1 = max(bile.wdia)*1.1;
   wdiamax2 = max(data.gene.wdia/1e6)*1.1;
   scalmax  = max(bile.ploss.*bile.scalts)*1.1;
	tmax     = t(param.gene.k);
	tmin     = t(param.gene.kmin);
	if ~isempty(jeux1)
	  tmax1  = t1(jeux1.param.gene.k);
	  tmax   = max(tmax,tmax1);
	end
   vgraph1  = [tmin tmax 0 max([ipmax phybmax pfcicon*1.1 nbarmax])];
   vgraph2  = [tmin tmax 0 max([wemax1 wemax2 rlwmax])];
   vgraph3  = [tmin tmax 0 max([wdiamax1 wdiamax2 scalmax])];
   tickgr   = linspace(round(tmin),round(tmax),5);
%
%
%
   h1=subplot(3,1,1);
   plot(t,data.gene.ip/1e6,t,data.gene.paddhyb/1e6,t,data.gene.nbar/1e19,t,data.gene.ihyb/1e6, ...
        t,data.gene.paddfci/1e6,t,pfcicon);
   set(h1,'xticklabel','')
   taille = get(gca,'fontsize');
   set(gca,'fontsize',12)   
   hleg=legend('Ip (MA)','PaddLH (MW)','nbar (10^1^9 m^-^2)','ILH (MA)','PaddFCI (MW)','PFCI cons',-1);
   set(gca,'fontsize',taille)
   posf = get(gca,'position');
   title(['shot ',int2str(fix(param.from.shot.num)),', B_0=',num2str(mean(data.geo.b0),2), ...
           ' T, R_0=',num2str(mean(data.geo.r0),3), ...
           ' m, a=',num2str(mean(data.geo.a),2),' m'])
   axis(vgraph1)
   set(h1,'xtick',tickgr)
   set(h1,'position',[0.15 0.69 0.5 0.25])
   set(hleg,'position',[0.66 0.71 0.3 0.25])
   grid
%
%
%
   h2=subplot(3,1,2);
   if isempty(t1) & ~isempty(bile.wetot)
      plot(times(1:3:end),bile.wetot(1:3:end),'bo',t,data.gene.we/1e6,'k',times,bile.rlw,'r')
      set(gca,'fontsize',12)   
      hleg=legend('We exp','We CRONOS','We RLW',-1);
      set(gca,'fontsize',taille)
   elseif ~isempty(t1) & ~isempty(bile.wetot)
      plot(times(1:3:end),bile.wetot(1:3:end),'bo',t,data.gene.we/1e6,'k',t1,jeux1.data.gene.we/1e6,'-.m',times,bile.rlw,'r')
      set(gca,'fontsize',12)
      hleg=legend('We exp','We CRONOS','reference','We RLW',-1);
      set(gca,'fontsize',taille)
   elseif ~isempty(t1) & isempty(bile.wetot)
      plot(t,data.gene.we/1e6,'k', ...
           t1,jeux1.data.gene.we/1e6,'-.m')
      set(gca,'fontsize',12)   
      hleg=legend('We CRONOS','reference');
      set(gca,'fontsize',taille)
   elseif isempty(t1) & isempty(bile.wetot)
      plot(t,data.gene.we/1e6,'k')
      set(gca,'fontsize',12)   
      hleg=legend('We CRONOS',-1);
      set(gca,'fontsize',taille)
   end
   set(h2,'xticklabel','')
   axis(vgraph2)  
   ylabel('MJ')
   pos2      = get(gca,'position');
   pos2(3:4) = posf(3:4);
   drawnow
   ymax = max(get(gca,'ylim'));
   set(gca,'ylim',[0,ymax]);
   set(h2,'xtick',tickgr) 
   grid  
   set(gca,'position',pos2);
   set(h2,'xtick',tickgr)
   set(h2,'position',[0.15 0.42 0.5 0.25])
   set(hleg,'position',[0.66 0.44 0.3 0.2])

%
%
%   
   h3=subplot(3,1,3);
   if isempty(t1) & ~isempty(bile.wdia)
      plot(times(1:3:end),bile.wdia(1:3:end),'bo',t,data.gene.wdia/1e6,'k',times,bile.ploss.*bile.scalts,'r')
      hold on
      plot(t,data.gene.wth/1e6,'g.-')
      hold off
      set(gca,'fontsize',12)   
      hleg=legend('Wdia exp.','Wdia CRONOS','TS scaling','Wther CRONOS',-1);
      set(gca,'fontsize',taille)      
   elseif ~isempty(t1) & ~isempty(bile.wdia)
      plot(times(1:3:end),bile.wdia(1:3:end),'bo',t,data.gene.wdia/1e6,'k', ...
           t1,jeux1.data.gene.wdia/1e6,'-.m',times,bile.ploss.*bile.scalts,'r')
      hold on
      plot(t,data.gene.wth/1e6,'g.-')
      hold off
      set(gca,'fontsize',12)   
      hleg=legend('Wdia exp','Wdia CRONOS','Wdia reference','Scaling TS','Wther CRONOS',-1);
      set(gca,'fontsize',taille)   
   elseif ~isempty(t1) & isempty(bile.wdia)
      plot(t,data.gene.wdia/1e6,'k', ...
           t1,jeux1.data.gene.wdia/1e6,'-.m')
      hold on
      plot(t,data.gene.wth/1e6,'g-')
      plot(t1,jeux1.data.gene.wth/1e6,'r.-')
      hold off
      set(gca,'fontsize',12)   
      hleg=legend('Wdia CRONOS','Wdia reference','Wther CRONOS','Wther ref',-1);
      set(gca,'fontsize',taille)   
   elseif isempty(t1) & isempty(bile.wdia)
      plot(t,data.gene.wdia/1e6,'k')
      hold on
      plot(t,data.gene.wth/1e6,'g.-')
      hold off
      set(gca,'fontsize',12)   
      hleg=legend('Wdia CRONOS','Wther CRONOS',-1);
      set(gca,'fontsize',taille)   
   end

   ylabel('MJ')
   xlabel('times (s)')
   set(h2,'xtick',tickgr)   
   axis(vgraph3)  
   pos3      = get(gca,'position');
   pos3(3:4) = posf(3:4);
   set(gca,'position',pos3);
  
   ymax = max(get(gca,'ylim'));
   set(gca,'ylim',[0,ymax]);
   set(h3,'xtick',tickgr)
   grid
   set(gcf,'ResizeFcn','')
   set(h3,'xtick',tickgr)
   set(h3,'position',[0.15 0.15 0.5 0.25])
   set(hleg,'position',[0.66 0.16 0.3 0.2])

%
%
%
   h = findobj(0,'type','figure','tag','comparetemp_2');
   if isempty(h)
       h=figure('tag','comparetemp_2');
   else
       figure(h);
   end
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   vgraph1  = [tmin tmax max([-0.5 min(data.gene.vloop) min(bile.vs)])*0.8 min([2,max([max(data.gene.vloop) max(bile.vs)])*1.1])];
   vgraph2  = [tmin tmax min([min(data.gene.li) min(bile.li)])*0.9 max([max(data.gene.li) max(bile.li)])*1.1];
   vgraph3  = [tmin tmax 0 max([max(data.gene.paddohm/1e6) max(bile.poh)])];
   tickgr   = linspace(round(tmin),round(max(t)),5);
%
   h1=subplot(3,1,1);
   taille = get(gca,'fontsize');
   if ~isempty(t1)& ~isempty(bile.vs)
     plot(times,bile.vs,'r',t,medfilt1(data.gene.vloop,3),'k')
	  hold on
	  plot(t1(1:3:end),jeux1.data.gene.vres(1:3:end),'-.m')
	  hold off
     set(gca,'fontsize',12)
	  hleg=legend('measurement','CRONOS','reference',-1);
	  set(gca,'fontsize',taille)
   elseif isempty(t1) & ~isempty(bile.li)
      plot(times,bile.vs,'r',t,medfilt1(data.gene.vloop,3),'k')
      set(gca,'fontsize',12)
      hleg=legend('measurement','CRONOS',-1);
	   set(gca,'fontsize',taille)
   elseif ~isempty(t1) & isempty(bile.li)
     plot(t,medfilt1(data.gene.vloop,3),'k')
	  hold on
	  plot(t1(1:3:end),jeux1.data.gene.vres(1:3:end),'-.m')
	  hold off
     set(gca,'fontsize',12)
     hleg=legend('CRONOS','reference',-1);
	  set(gca,'fontsize',taille)
  elseif isempty(t1) & isempty(bile.li)
     plot(t,medfilt1(data.gene.vloop,3),'k')
     set(gca,'fontsize',12)
     hleg=legend('CRONOS',-1);
	  set(gca,'fontsize',taille)
   end
   grid
   posf = get(gca,'position');
   axis(vgraph1)
   set(h1,'xtick',tickgr)
   title(['shot ',int2str(fix(param.from.shot.num)),', B_0=',num2str(mean(data.geo.b0),2), ...
           ' T, R_0=',num2str(mean(data.geo.r0),3), ...
           ' m, a=',num2str(mean(data.geo.a),2),' m'])
   ymax = max(get(gca,'ylim'));
   set(gca,'ylim',[0,ymax]);
   ylabel('V_l_o_o_p (V)')
   set(h1,'xticklabel','')
   set(h1,'position',[0.15 0.67 0.5 0.25])
   set(hleg,'position',[0.66 0.71 0.3 0.1],'fontsize',12)

%   
   h2=subplot(3,1,2);
   set(h2,'xticklabel','')
   taille = get(gca,'fontsize');
   if ~isempty(t1)& ~isempty(bile.li)
     plot(times,bile.li,'r',t,data.gene.li,'k')
	  hold on
	  plot(t1(1:3:end),jeux1.data.gene.li(1:3:end),'-.m')
	  hold off
     set(gca,'fontsize',12)
	  hleg=legend('measurement','CRONOS','reference',-1);
	  set(gca,'fontsize',taille)
   elseif isempty(t1)& ~isempty(bile.li)
     plot(times,bile.li,'r',t,data.gene.li,'k')
     set(gca,'fontsize',12)
     hleg=legend('measurement','CRONOS',-1);
	  set(gca,'fontsize',taille)
   elseif ~isempty(t1)& isempty(bile.li)
     plot(t,data.gene.li,'k')
 	  hold on
	  plot(t1(1:3:end),jeux1.data.gene.li(1:3:end),'-.m')
	  hold off
     set(gca,'fontsize',12)
     hleg=legend('CRONOS','reference',-1);
	  set(gca,'fontsize',taille)
   elseif isempty(t1)& isempty(bile.li)
     plot(t,data.gene.li,'k')
     set(gca,'fontsize',12)
     hleg=legend('CRONOS',-1);
	  set(gca,'fontsize',taille)
   end
   grid
   posf = get(gca,'position');
   axis(vgraph2)
   set(h2,'xtick',tickgr)
   ymax = max(get(gca,'ylim'));
   set(gca,'ylim',[0,ymax]);
   ylabel('li ')
   pos2      = get(gca,'position');
   pos2(3:4) = posf(3:4);
   set(gca,'position')
   set(h2,'xticklabel','')
   set(h2,'position',[0.15 0.41 0.5 0.25])
   set(hleg,'position',[0.66 0.45 0.3 0.1],'fontsize',12)

%
%   
   h3=subplot(3,1,3);
   set(h3,'xticklabel','')
   taille = get(gca,'fontsize');   
   if ~isempty(t1) & ~isempty(bile.poh)
       plot(times,bile.poh,'r',t,data.gene.paddohm/1e6,'k')
       hold on
       plot(t1(1:3:end),jeux1.data.gene.paddohm(1:3:end)/1e6,'-.m')
       hold off
       set(gca,'fontsize',12)
       hleg=legend('measurement','CRONOS','reference',-1)
       set(gca,'fontsize',taille)
   elseif isempty(t1) & ~isempty(bile.poh)
       plot(times,bile.poh,'r',t,data.gene.paddohm/1e6,'k')
       set(gca,'fontsize',12)
       hleg=legend('measurement','CRONOS',-1)
       set(gca,'fontsize',taille)
   elseif ~isempty(t1) & isempty(bile.poh)
       plot(t,data.gene.paddohm/1e6,'k')
       hold on
       plot(t1(1:3:end),jeux1.data.gene.paddohm(1:3:end)/1e6,'-.m')
       hold off      
		 set(gca,'fontsize',12)
       hleg=legend('CRONOS','reference',-1)
       set(gca,'fontsize',taille)
   elseif isempty(t1) & isempty(bile.poh)
       plot(t,data.gene.paddohm/1e6,'k')
       set(gca,'fontsize',12)
       hleg=legend('CRONOS',-1);       
       set(gca,'fontsize',taille)
   end
   set(h3,'position',[0.15 0.15 0.5 0.25])
   set(hleg,'position',[0.66 0.2 0.3 0.1],'fontsize',12)
   
   ylabel('P_\Omega (MW)')
   ymax = max(get(gca,'ylim'));
   grid
   set(gca,'ylim',[0,ymax]);
   axis(vgraph3)   
   xlabel('times (s)')
   pos3      = get(gca,'position');
   pos3(3:4) = posf(3:4);
   set(gca,'position')
   set(gcf,'ResizeFcn','')

%
%
%
   h = findobj(0,'type','figure','tag','comparetemp_3');
   if isempty(h)
       h=figure('tag','comparetemp_3');
   else
       figure(h);
   end
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1]);

   if ~isempty(t1) & ~isempty(bile.rmaj)
        plot(t,data.equi.d(:,1),'k',t1,jeux1.data.equi.d(:,1),'-.m',times,bile.rm(:,1)-bile.rmaj,'r')
	     legend('CRONOS','reference','measurement')
   elseif isempty(t1) & ~isempty(bile.rmaj)
        plot(t,data.equi.d(:,1),'k',times,bile.rm(:,1)-bile.rmaj,'r')
	     legend('CRONOS','measurement')
   elseif ~isempty(t1) & isempty(bile.rmaj)
        plot(t,data.equi.d(:,1),'k',t1,jeux1.data.equi.d(:,1),'-.m')
	     legend('CRONOS','reference')
   elseif isempty(t1) & isempty(bile.rmaj)
        plot(t,data.equi.d(:,1))
	     legend('CRONOS')	
   end
   title('Axes magetique -grand rayon')
   xlabel('temps (s)')
   ylabel('d0 (m)')
	axis([tmin tmax 0 0.4])

case 'JET'

   h = findobj(0,'type','figure','tag','comparetemp_1');
   if isempty(h)
       h=figure('tag','comparetemp_1');
   else
       figure(h);
   end
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   if ~isempty(jeux1)
     plot(t,data.gene.ip,'r',t1,jeux1.data.gene.ip,':r',jettemp.tIp,-jettemp.Ip,'or', ...
          t,data.gene.iboot,'b',t1,jeux1.data.gene.iboot,':b',...
          t,data.gene.ihyb,'c',t1,jeux1.data.gene.ihyb,':c',...
          t,data.gene.iidn,'m',t1,jeux1.data.gene.iidn,':m');
   else
     plot(t,data.gene.ip,'r',jettemp.tIp,-jettemp.Ip,'or', ...
          t,data.gene.iboot,'b',...
          t,data.gene.ihyb,'c',...
          t,data.gene.iidn,'m');
   end
   if ~isempty(jetefit) & ~isempty(jettransp)
      hold on
      tidn  = jetefit.tefit;
      jidn  = interp1(jettransp.ttransp,jettransp.jnbtrx,jetefit.tefit);
      spr   = jetefit.dSdrhonx;
      iidn  = trapz(jetefit.xefit',spr'.*jidn');
      plot(tidn,iidn,'mo');
      if ~isempty(jeux1)
        plot(t,data.gene.ini,'g',t1,jeux1.data.gene.ini,':g');
      else
        plot(t,data.gene.ini,'g');
      end
      hold off
      if ~isempty(jeux1)
        legend('Ip cronos','Ip reference','Ip efit','Iboot cronos','Iboot reference','Ihyb cronos','Ihyb reference', ...
               'Iidn cronos','Iidn reference','Iidn transp','Ini cronos','Ini reference');
      else
        legend('Ip cronos','Ip efit','Iboot cronos','Ihyb cronos', ...
               'Iidn cronos','Iidn transp','Ini cronos');      
      end
   else
      hold on
      if ~isempty(jeux1)
          plot(t,data.gene.ini,'g',t1,jeux1.data.gene.ini,':g');
	   else
		    plot(t,data.gene.ini,'g');
		end
      hold off
      if ~isempty(jeux1)
        legend('Ip cronos','Ip reference','Ip efit','Iboot cronos','Iboot reference','Ihyb cronos','Ihyb reference', ...
               'Iidn cronos','Iidn reference','Ini cronos','Ini reference');
      else
        legend('Ip cronos','Ip efit','Iboot cronos','Ihyb cronos', ...
               'Iidn cronos','Ini cronos');
      end
   end
   axis([param.gene.tdeb,param.gene.tfin,0,inf]) 
   title(sprintf('Courant plasma et courants generes (%s@%d)',param.from.machine,fix(param.from.shot.num)))
   xlabel('temps (s)')
   ylabel('I (A)')   
   set(gcf,'ResizeFcn','')

   h = findobj(0,'type','figure','tag','comparetemp_2');
   if isempty(h)
       h=figure('tag','comparetemp_2');
   else
       figure(h);
   end
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   plotprof(gca,t,x,data.prof.jmoy,'color','r');
   if ~isempty(jeux1)
     plotprof(gca,t1,x,jeux1.data.prof.jmoy,'color','r','linestyle',':');
   end
   if ~isempty(jetefit)
     plotprof(gca,jetefit.tefit,jetefit.xefit,-real(jetefit.Jx),'linestyle','none','marker','o','color','r');
   end
   plotprof(gca,t,x,data.neo.jboot,'color','b');
   if ~isempty(jeux1)
     plotprof(gca,t1,x,jeux1.data.neo.jboot,'color','b','linestyle',':');
   end
   plotprof(gca,t,x,data.source.hyb.j,'color','c');
   if ~isempty(jeux1)
     plotprof(gca,t1,x,jeux1.data.source.hyb.j,'color','c','linestyle',':');
   end
   plotprof(gca,t,x,data.source.idn.j,'color','m');
   if ~isempty(jeux1)
     plotprof(gca,t1,x,jeux1.data.source.idn.j,'color','m','linestyle',':');
   end
   if ~isempty(jettransp)
      plotprof(gca,jettransp.ttransp,jetefit.xefit,jettransp.jnbtrx,'linestyle','none','marker','o','color','m');
   end
   plotprof(gca,t,x,data.prof.jmoy-data.source.totale.j,'color','g');
   if ~isempty(jeux1)
     plotprof(gca,t1,x,jeux1.data.prof.jmoy-jeux1.data.source.totale.j,'color','g','linestyle',':');
   end
   if isempty(jetefit) & isempty(jettransp)
       if ~isempty(jeux1)
         legend('Jmoy cronos','Jmoy reference','Jboot cronos','Jboot reference','Jhyb cronos','Jhyb reference', ...
                'Jidn cronos','Jidn reference','Johm cronos','Johm reference');
       else
         legend('Jmoy cronos','Jboot cronos','Jhyb cronos', ...
                'Jidn cronos','Johm cronos');
       end
   elseif isempty(jettransp)
       if ~isempty(jeux1)   
         legend('Jmoy cronos','Jmoy reference','Jmoy EFIT','Jboot cronos','Jboot reference','Jhyb cronos','Jhyb reference', ...
                'Jidn cronos','Jidn reference','Johm cronos','Johm reference');  
       else
         legend('Jmoy cronos','Jmoy EFIT','Jboot cronos','Jhyb cronos', ...
                'Jidn cronos','Johm cronos');  
       end
   else
        if ~isempty(jeux1)
          legend('Jmoy cronos','Jmoy reference','Jmoy EFIT','Jboot cronos','Jboot reference','Jhyb cronos','Jhyb reference', ...
                'Jidn cronos','Jidn reference','Jidn TRANSP','Johm cronos','Johm reference');
        else
          legend('Jmoy cronos','Jmoy EFIT','Jboot cronos','Jhyb cronos', ...
                'Jidn cronos','Jidn TRANSP','Johm cronos');	
	end
   end
   title(sprintf('Densite de courant plasma et de courants generes (%s@%d)',param.from.machine,fix(param.from.shot.num)))
   xlabel('sqrt(Phi/pi/B0) normalise')
   ylabel('J (A*m^-^2)')
   set(gcf,'ResizeFcn','')
   
   
   
   h = findobj(0,'type','figure','tag','comparetemp_3');
   if isempty(h)
       h=figure('tag','comparetemp_3');
   else
       figure(h);
   end
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   if ~isempty(jeux1)
     plot(t,data.gene.li,'r',t1,jeux1.data.gene.li,':r',jettemp.tli,jettemp.li,'or')
     legend('CRONOS','reference','EFIT')
   else
     plot(t,data.gene.li,'r',jettemp.tli,jettemp.li,'or')
     legend('CRONOS','EFIT')
   end
   axis([param.gene.tdeb,param.gene.tfin,0,2]) 
   title(sprintf('Induction (%s@%d)',param.from.machine,fix(param.from.shot.num)))
   xlabel('temps (s)')
   ylabel('li')   

      
   h = findobj(0,'type','figure','tag','comparetemp_4');
   if isempty(h)
       h=figure('tag','comparetemp_4');
   else
       figure(h);
   end
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   if ~isempty(jeux1)
     plot(t,-data.gene.vloop,'r',t1,-jeux1.data.gene.vloop,':r',jettemp.tvsur,jettemp.vsur,'b')
     legend('CRONOS','reference','EFIT') 
   else
     plot(t,-data.gene.vloop,'r',jettemp.tvsur,jettemp.vsur,'b')
     legend('CRONOS','EFIT') 
   end
   axis([param.gene.tdeb,param.gene.tfin,-2,2]) 
   title(sprintf('tension par tour (%s@%d)',param.from.machine,fix(param.from.shot.num)))
   xlabel('temps (s)')
   ylabel('V_l_o_o_p')   
   set(gcf,'ResizeFcn','')

case 'DIIID'

   h = findobj(0,'type','figure','tag','comparetemp_1');
   if isempty(h)
       h=figure('tag','comparetemp_1');
   else
       figure(h);
   end
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   if ~isempty(jeux1)
     plot(t,data.gene.ip,'r',t1,jeux1.data.gene.ip,':r',diiidtemp.tip,diiidtemp.ip,'or', ...
          t,data.gene.iboot,'b',t1,jeux1.data.gene.iboot,':b',...
           t,data.gene.iidn,'m',t1,jeux1.data.gene.iidn,':m');
   else
     plot(t,data.gene.ip,'r',diiidtemp.tip,diiidtemp.ip,'or', ...
          t,data.gene.iboot,'b',...
           t,data.gene.iidn,'m');
   end
   if ~isempty(diiidtemp) 
      hold on
      if ~isempty(jeux1)
        plot(t,data.gene.ini,'g',t1,jeux1.data.gene.ini,':g');
      else
        plot(t,data.gene.ini,'g');
      end
      hold off
      if ~isempty(jeux1)
        legend('Ip cronos','Ip reference','Ip efit','Iboot cronos','Iboot reference', ...
               'Inbi cronos','Inbi reference','Ini cronos','Ini reference');
      else
        legend('Ip cronos','Ip efit','Iboot cronos', ...
               'Inbi cronos','Ini cronos');      
      end
   else
      hold on
      if ~isempty(jeux1)
          plot(t,data.gene.ini,'g',t1,jeux1.data.gene.ini,':g');
	   else
		    plot(t,data.gene.ini,'g');
		end
      hold off
      if ~isempty(jeux1)
        legend('Ip cronos','Ip reference','Ip efit','Iboot cronos','Iboot reference', ...
               'Inbi cronos','Inbi reference','Ini cronos','Ini reference');
      else
        legend('Ip cronos','Ip efit','Iboot cronos', ...
               'Inbi cronos','Ini cronos');
      end
   end
   axis([param.gene.tdeb,param.gene.tfin,0,inf]) 
   title(sprintf('Currents (%s@%d)',param.from.machine,fix(param.from.shot.num)))
   xlabel('time (s)')
   ylabel('I (A)')   
   set(gcf,'ResizeFcn','')

   h = findobj(0,'type','figure','tag','comparetemp_2');
   if isempty(h)
       h=figure('tag','comparetemp_2');
   else
       figure(h);
   end
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   plotprof(gca,t,x,data.prof.jmoy,'color','r');
   if ~isempty(jeux1)
     plotprof(gca,t1,x,jeux1.data.prof.jmoy,'color','r','linestyle',':');
   end
   plotprof(gca,t,x,data.neo.jboot,'color','b');
   if ~isempty(jeux1)
     plotprof(gca,t1,x,jeux1.data.neo.jboot,'color','b','linestyle',':');
   end
   plotprof(gca,t,x,data.source.idn.j,'color','m');
   if ~isempty(jeux1)
     plotprof(gca,t1,x,jeux1.data.source.idn.j,'color','m','linestyle',':');
   end
   plotprof(gca,t,x,data.prof.jmoy-data.source.totale.j,'color','g');
   if ~isempty(jeux1)
     plotprof(gca,t1,x,jeux1.data.prof.jmoy-jeux1.data.source.totale.j,'color','g','linestyle',':');
   end
       if ~isempty(jeux1)
         legend('Jmoy cronos','Jmoy reference','Jboot cronos','Jboot reference', ...
                'Jnbi cronos','Jnbi reference','Johm cronos','Johm reference');
       else
         legend('Jmoy cronos','Jboot cronos', ...
                'Jnbi cronos','Johm cronos');
       end

   title(sprintf('current profiles (%s@%d)',param.from.machine,fix(param.from.shot.num)))
   xlabel('sqrt(Phi/pi/B0) normalized')
   ylabel('J (A*m^-^2)')
   set(gcf,'ResizeFcn','')
   
   
   
   h = findobj(0,'type','figure','tag','comparetemp_3');
   if isempty(h)
       h=figure('tag','comparetemp_3');
   else
       figure(h);
   end
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   if ~isempty(jeux1)
     plot(t,data.gene.li,'r',t1,jeux1.data.gene.li,':r',diiidtemp.tli,diiidtemp.li,'or')
     legend('CRONOS','reference','EFIT')
   else
     plot(t,data.gene.li,'r',diiidtemp.tli,diiidtemp.li,'or')
     legend('CRONOS','EFIT')
   end
   axis([param.gene.tdeb,param.gene.tfin,0,2]) 
   title(sprintf('Induction (%s@%d)',param.from.machine,fix(param.from.shot.num)))
   xlabel('time (s)')
   ylabel('li')   

      
   h = findobj(0,'type','figure','tag','comparetemp_4');
   if isempty(h)
       h=figure('tag','comparetemp_4');
   else
       figure(h);
   end
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   clf
   set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
   if ~isempty(jeux1)
     plot(t,-data.gene.vloop,'r',t1,-jeux1.data.gene.vloop,':r',jettemp.tvsur,jettemp.vsur,'b')
     legend('CRONOS','reference','EFIT') 
   else
     plot(t,-data.gene.vloop,'r',diiidtemp.tvl,diiidtemp.vl,'b')
     legend('CRONOS','EFIT') 
   end
   axis([param.gene.tdeb,param.gene.tfin,-2,2]) 
   title(sprintf('tension par tour (%s@%d)',param.from.machine,fix(param.from.shot.num)))
   xlabel('time (s)')
   ylabel('V_l_o_o_p')   
   set(gcf,'ResizeFcn','')
   
   
otherwise
        disp('no input data for this tokamak')
end

