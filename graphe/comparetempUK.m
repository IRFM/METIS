% comparaison des donnees de la chaleur en fonction du temps
chargedata





   h = findobj(0,'type','figure','tag','comparetemp_1');
   if isempty(h)
       h=figure('tag','comparetemp_1');
   else
       figure(h);
   end
   clf
   set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','normal','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
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
   vgraph0  = [tmin tmax 0 ipmax*1.1];
   vgraph2  = [tmin tmax 0 max([wemax1 wemax2 rlwmax])];
   vgraph3  = [tmin tmax 0 max([wdiamax1 wdiamax2 scalmax])];
   tickgr   = linspace(round(tmin),round(tmax),5);
%
%
%
   h1=subplot(3,1,1);
   plot(times,bile.ip)
   set(h1,'xticklabel','')
   taille = get(gca,'fontsize');
   set(gca,'fontsize',12)   
	text('units','normalized','position',[1.05 0.5],'string','a)');
   axis(vgraph0)
	
	
   h2=subplot(3,1,2);
	
	plot(t,data.gene.paddhyb/1e6,'r-',t,data.gene.nbar/1e19,'b.-', ...
        t,abs(pfcicon),'g--');
   set(h2,'xticklabel','')
   taille = get(gca,'fontsize');
   set(gca,'fontsize',12)   
   posf = get(gca,'position');
   st=['shot ',int2str(fix(param.from.shot.num)),', B_0=',num2str(mean(data.geo.b0),2), ...
           ' T, R_0=',num2str(mean(data.geo.r0),3), ...
           ' m, a=',num2str(mean(data.geo.a),2),' m'];
   axis(vgraph1)
   set(h2,'xtick',tickgr)
   
	text('units','normalized','position',[1.05 0.5],'string','b)');

%
%
%
   h3=subplot(3,1,3);

      plot(times(1:3:end),bile.wdia(1:3:end),'r')
      hold on
   axis(vgraph3)  


   ylabel('MJ')
   xlabel('times (s)')

text('units','normalized','position',[1.05 0.5],'string','c)');
