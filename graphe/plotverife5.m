% script de  la figure de verification de l'execution de zineb en diffusion du courant
%chargedata
x=param.gene.x;
t=data.gene.temps;

if length(t)> 100
	pas =fix(length(t)/100);
	if pas ==0
		pas =1;
	end
else
	pas =1;
end
if param.gene.k > param.gene.kmin
  ind  = param.gene.kmin:min(param.gene.k,length(t));
else
  ind  = param.gene.kmin:length(t);
end
tt   = t(ind);
inde = 1:pas:length(t);
t    = t(inde);
tmin = t(1);

h = findobj(0,'type','figure','tag','compare');
if isempty(h)
       h=figure('tag','compare');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])

switch param.from.machine

case 'TS'
imax = find(isnan(data.gene.ip(ind)) == 1);
if ~isempty(imax)
 imax    = imax(1);
else
 imax = length(tt);
end
if exist('jeux1')
  if ~isempty(jeux1)
    imax1 = find(isnan(jeux1.data.gene.ip(ind1)) == 1);
    if ~isempty(imax1)
      imax1    = imax1(1);
    else
      imax1 = length(tt1);
    end
    tmax = max(tt(imax),tt1(imax1));
  end
else
    tmax = tt(imax);
 end
tmax=100;
if isempty(bile.ip)
  vgraph1 = [tmin tmax 0 max(data.gene.ip(ind)/1e6)*1.1];
elseif max(data.exp.ip(inde))/1e6 > 0.4
  vgraph1 = [tmin tmax 0 max(data.exp.ip(inde)/1e6)*1.1];
else
  vgraph1 = [tmin tmax 0 max(data.cons.ip(inde)/1e6)*1.1];
end
vgraph2 = [tmin 100 0 2];
vgraph3 = [tmin 100 0.5 2];
ind95   =length(param.gene.x);

vgraph4 = [tmin tmax 2 max(data.prof.q(ind,ind95))*1.2];
tickgr  = linspace(0,round(tt(imax)),6);

h1=subplot(4,1,1);
if ~isempty(bile.ip)
  plot(t,data.exp.ip(inde)/1e6,'r',tt,data.gene.ip(ind)/1e6,'r.',tt,data.gene.iboot(ind)/1e6,'.m', ...
       tt,data.gene.ini(ind)/1e6,'k-');
  set(h1,'xticklabel','')
  axis(vgraph1)
  set(h1,'xtick',tickgr)
  taille = get(gca,'fontsize');
  set(gca,'fontsize',12)
  legend(gca,'Ip mesure','Ip CRONOS','I boot','I NI',2)
  ylabel('MA')
else
  plot(tt,data.gene.ip(ind)/1e6,'r',tt,data.gene.iboot(ind)/1e6,'m', ...
       tt,data.gene.ini(ind)/1e6,'k');
  set(h1,'xticklabel','')
  axis(vgraph1)
  set(h1,'xtick',tickgr)
  taille = get(gca,'fontsize');
  set(gca,'fontsize',12)
  if exist('jeux1','var')
    hold on
    plot(tt1,jeux1.data.gene.ip(ind1)/1e6,'c.')
    plot(tt1,jeux1.data.gene.ini(ind1)/1e6,'b.')
    hold off
%    legend(gca,'Ip CRONOS','I boot','I NI','Ip ref','I NI ref',2)
  else
%    legend(gca,'Ip CRONOS','I boot','I NI',2)
end
end
set(gca,'fontsize',taille)

if ~isempty(bile.ip)
  st =sprintf('choc %s #%d, t = %4.2g:%4.2g:%4.2g, nbrho = %d',param.from.machine, ...   
             fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
else
  st =sprintf('choc %s simulation, t = %4.2g:%4.2g:%4.2g, nbrho = %d',param.from.machine, ...   
             param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);

end
title(st);
grid

h2=subplot(4,1,2);

if ~isempty(bile.ip)
  plot(t,data.exp.vloop(inde),'r',tt,data.gene.vloop(ind),'m', ...
       tt,medfilt1(data.gene.vsurf(ind),5),'k');
  set(h2,'xticklabel','')
  ylabel('Vloop (V)')
  set(gca,'ylim',[0 2])
  if exist('jeux1','var')
    hold on
    plot(tt1,jeux1.data.gene.vloop(ind1),'r.')
    hold off
  end
  if exist('jeux1','var')
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)
    legend(gca,'mesure','boucle','surface','reference',2)
    set(gca,'fontsize',taille)
  else
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)
    legend(gca,'mesure','boucle','surface',2)
    set(gca,'fontsize',taille)
  end
else
  plot(tt,data.gene.vloop(ind),'m', ...
       tt,medfilt1(data.gene.vsurf(ind),5),'k');
  set(h2,'xticklabel','')
  ylabel('Vloop (V)')
  set(gca,'ylim',[0 2])
  if exist('jeux1','var')
    hold on
    plot(tt1,jeux1.data.gene.vloop(ind1),'r.')
    hold off
  end
  if exist('jeux1','var')
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)
    legend(gca,'boucle','surface','reference',2)
    set(gca,'fontsize',taille)
  else
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)
    legend(gca,'boucle','surface',2)
    set(gca,'fontsize',taille)
  end
end
axis(vgraph2)
set(h2,'xtick',tickgr)
grid

h3=subplot(4,1,3);
if ~isempty(bile.ip)
  plot(bile.times,bile.beli,'r',tt,data.gene.betap(ind)+data.gene.li(ind)/2,'b.-')
  if exist('jeux1')
    hold on
    plot(tt1,jeux1.data.gene.betap(ind1)+jeux1.data.gene.li(ind1)/2,'r:')
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)  
    legend(gca,'mesure','CRONOS','reference',2)
    set(gca,'fontsize',taille)
    hold off
  else
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)
    legend(gca,'mesure','CRONOS',2)
    set(gca,'fontsize',taille)
  end
else
  plot(tt,data.gene.betap(ind)+data.gene.li(ind)/2,'b.-')
  if exist('jeux1')
    hold on
    plot(tt1,jeux1.data.gene.betap(ind1)+jeux1.data.gene.li(ind1)/2,'r:')
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)  
    legend(gca,'CRONOS','reference',2)
    set(gca,'fontsize',taille)
    hold off
  else
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)
    set(gca,'fontsize',taille)
  end
end
set(h3,'xticklabel','')
ylabel('beta+li/2')
axis(vgraph3)
set(h3,'xtick',tickgr)
grid

h4=subplot(4,1,4);

disp('q au bord')
ind95 =length(param.gene.x);
if ~isempty(bile.ip)
  plot(t,data.exp.qa(inde),'r',tt,data.prof.q(ind,ind95),'b.-')
  if exist('jeux1')
    hold on
    plot(tt1,jeux1.data.prof.q(ind1,ind95),'r:')
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)  
    legend(gca,'mesure','CRONOS','reference',2)
    set(gca,'fontsize',taille)
    hold off
  else
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)  
    legend(gca,'mesure','CRONOS',2)
    set(gca,'fontsize',taille)    
  end
else
  plot(tt,data.prof.q(ind,ind95),'b.-')
  if exist('jeux1')
    hold on
    plot(tt1,jeux1.data.prof.q(ind1,ind95),'r:')
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)  
    legend(gca,'CRONOS','reference',2)
    set(gca,'fontsize',taille)
    hold off
  else
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)  
    set(gca,'fontsize',taille)    
  end
end
ylabel('q_a')
axis(vgraph4)
set(h4,'xtick',tickgr)

grid
xlabel('temps (s)')
set(gcf,'color',[1 1 1])
set(gcf,'ResizeFcn','')

case 'JET'

imax = find(isnan(data.gene.ip(ind)) == 1);
if ~isempty(imax)
 imax    = imax(1);
 if imax < 10
   imax = length(tt);
 end
else
 imax = length(tt);
end

if tt(imax) < min(jettemp.tli)
  indjetli  = iround(jettemp.tli,tt(imax)+40);
else
  indjetli  = iround(jettemp.tli,tt(imax));
end
elar = (tt(imax)-tt(1))/10;
vgraph1 = [tt(1)-elar tt(imax) 0 max(data.exp.ip(inde)/1e6)*1.3];
vgraph2 = [tt(1)-elar tt(imax) -1 2];
if ~isempty(jettemp.betap)
  mvgraph3 = max(jettemp.betap(1:indjetli)+jettemp.li(1:indjetli)/2);
else
  mvgraph3 = max(jettemp.li(1:indjetli));
end
vgraph3 = [tt(1)-elar tt(imax) 0.5 mvgraph3*1.2];
mvgraph4 = max(jettemp.li(1:indjetli));
vgraph4 = [tt(1)-elar tt(imax) 0.5 mvgraph4*1.2];
ind95   =length(param.gene.x);
vgraph5 = [tt(1)-elar tt(imax) 2 max(data.prof.q(ind,ind95))*1.2];
tickgr  = round(linspace(tt(1),round(tt(imax)),5)*10)/10;



h1=subplot(5,1,1);
plot(t,data.exp.ip(inde)/1e6,'r',tt,data.gene.ip(ind)/1e6,'r.',tt,data.gene.iboot(ind)/1e6,'om', ...
     tt,data.gene.ini(ind)/1e6,'k-')
taille = get(gca,'fontsize');
if exist('jeux1','var')
  hold on
  plot(tt1,jeux1.data.gene.ini(ind1)/1e6,'k:')
  hold off
end
set(gca,'fontsize',12)  
if ~exist('jeux1','var')
  legend(gca,'Ip EFIT','Ip CRONOS','I boot','I NI',2)
else
  legend(gca,'Ip EFIT','Ip CRONOS','I boot','I NI','I NI ref.',2)
end
set(gca,'fontsize',taille)  
pos1 =get(gca,'position');
ylabel('MA')
set(h1,'xticklabel','')
axis(vgraph1)
set(h1,'xtick',tickgr)
st =sprintf('choc %s #%d, t = %4.2g:%4.2g:%4.2g, nbrho = %d',param.from.machine, ...   
             fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
title(st);
grid

h2=subplot(5,1,2);
si = param.gene.signe.ip*param.gene.signe.b0;
plot(t,data.exp.vloop(inde)*si,'r',tt,data.gene.vloop(ind),'b');
if exist('jeux1','var')
  hold on
  plot(tt1,jeux1.data.gene.vloop(ind1),'r.')
  hold off
end
if exist('jeux1','var')
  taille = get(gca,'fontsize');
  set(gca,'fontsize',12)  
  legend(gca,'EFIT','CRONOS','reference',2)
  set(gca,'fontsize',taille) 
else
  taille = get(gca,'fontsize');
  set(gca,'fontsize',12)  
  legend(gca,'EFIT','CRONOS',2)
  set(gca,'fontsize',taille)  
end     
pos2 =get(gca,'position');
pos2(3:4) = pos1(3:4);
set(gca,'position',pos2)
ylabel('Vloop (V)')
set(h2,'xticklabel','')
axis(vgraph2)
set(h2,'xtick',tickgr)
grid
%
%=====================
%
h3=subplot(5,1,3);
if ~isempty(jettemp.betap)
  if max(tt) < min(jettemp.tli)
    plot(jettemp.tli-40,jettemp.betap+jettemp.li/2,'r',tt,data.gene.betap(ind)+data.gene.li(ind)/2,'b')
  else
    plot(jettemp.tli,jettemp.betap+jettemp.li/2,'r',tt,data.gene.betap(ind)+data.gene.li(ind)/2,'b')
  end
else
  if max(tt) < min(jettemp.tli)
    plot(jettemp.tli-40,jettemp.li,'r',tt,data.gene.li(ind),'b')
  else
    plot(jettemp.tli,jettemp.li,'r',tt,data.gene.li(ind),'b')
  end
end
if exist('jeux1','var')
  hold on
  plot(tt1,jeux1.data.gene.betap(ind1)+jeux1.data.gene.li(ind1)/2,'r.')
  hold off
end
if exist('jeux1','var')
  taille = get(gca,'fontsize');
  set(gca,'fontsize',12)
  legend(gca,'EFIT','CRONOS','reference',2)
  set(gca,'fontsize',taille)  
else
  taille = get(gca,'fontsize');
  set(gca,'fontsize',12)
  legend(gca,'EFIT','CRONOS',2)
  set(gca,'fontsize',taille)  
end
if ~isempty(jettemp.betap)
  ylabel('beta+li/2')
else
  ylabel('li')
end
pos2 =get(gca,'position');
pos2(3:4) = pos1(3:4);
set(gca,'position',pos2)
set(h3,'xticklabel','')
axis(vgraph3)
set(h3,'xtick',tickgr)
grid
%
%=====================
%
h4=subplot(5,1,4);
disp('Provisoire : jettemp.tli - 40 s')
if max(tt) < min(jettemp.tli)
  plot(jettemp.tli-40,jettemp.li,'r',tt,data.gene.li(ind),'b')
else
  plot(jettemp.tli,jettemp.li,'r',tt,data.gene.li(ind),'b')
end
if exist('jeux1','var')
  hold on
  plot(tt1,jeux1.data.gene.li(ind1),'r.')
  hold off
end
if exist('jeux1','var')
  taille = get(gca,'fontsize');
  set(gca,'fontsize',12)
  legend(gca,'EFIT','CRONOS','reference',2)
  set(gca,'fontsize',taille)  
else
  taille = get(gca,'fontsize');
  set(gca,'fontsize',12)
  legend(gca,'EFIT','CRONOS',2)
  set(gca,'fontsize',taille)  
end
ylabel('li')
pos2 =get(gca,'position');
pos2(3:4) = pos1(3:4);
set(gca,'position',pos2)
set(h4,'xticklabel','')
axis(vgraph4)
set(h4,'xtick',tickgr)
grid
%
%=====================
%


h5=subplot(5,1,5);
ind95 =  min( iround(param.gene.x,0.95));

plot(t,data.exp.qa(inde),'r',tt,data.prof.q(ind,ind95),'b')
taille = get(gca,'fontsize');
set(gca,'fontsize',12)
legend(gca,'EFIT','CRONOS',2)
set(gca,'fontsize',taille)  
pos2 =get(gca,'position');
pos2(3:4) = pos1(3:4);
set(gca,'position',pos2)
axis(vgraph5)
set(h5,'xtick',tickgr)

ylabel('q_9_5')

grid
xlabel('temps (s)')

set(gcf,'color',[1 1 1])
set(gcf,'ResizeFcn','')
otherwise
        disp('pas de donnees pour cette machine')
end
