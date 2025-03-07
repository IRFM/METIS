% script de  la figure de verification de l'execution de zineb en diffusion du courant
chargedata
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

ind  = param.gene.kmin:min(param.gene.k,length(t));
tt   = t(ind);
%inde = 1:pas:length(t);
%t    = t(inde);

if exist('jeux1','var')
  if ~isempty(jeux1)
    ind1=jeux1.param.gene.kmin:min(jeux1.param.gene.k,length(t1));
    tt1=t1(ind1);
  else
    jeux1 = [];
  end
else
  jeux1 = [];
end

h = findobj(0,'type','figure','tag','compare');
if isempty(h)
       h=figure('tag','compare');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','normal','defaultaxesfontname','times', ...
	'defaultlinelinewidth',2,'color',[1 1 1],'defaultlinemarkersize',4)

imax = find(isnan(data.gene.ip(ind)) == 1);
if ~isempty(imax)
 imax    = imax(1);
else
 imax = length(tt);
end
if ~isempty(jeux1)
    imax1 = find(isnan(jeux1.data.gene.ip(ind1)) == 1);
    if ~isempty(imax1)
      imax1    = imax1(1);
    else
      imax1 = length(tt1);
    end
    tmax = max(tt(imax),tt1(imax1));
	 tmin = min(param.gene.tdeb,jeux1.param.gene.tdeb);
else
    tmax = tt(imax);
	 tmin = param.gene.tdeb;
end

vgraph1 = [tmin tmax 0 max(data.gene.ip(ind)/1e6)*1.1];
vgraph2 = [tmin tmax -0.1 0.7];
vgraph3 = [tmin tmax 0.5 2.0];
ind95   =length(param.gene.x);
qmax    = max(data.prof.q(:));
if ~isempty(jeux1)
  qmax   = max(qmax,max(jeux1.data.prof.q(:)));
end

vgraph4 = [tmin tmax 9 10];
tickgr  = fix(linspace(0,round(tt(imax)),10));
if ~isempty(jeux1)
  tickgr  = fix(linspace(0,round(max(tt(imax),tt1(imax1))),10));
end
phybmax = max(data.gene.paddhyb/1e6)*1.1;
pfcimax = max(data.gene.paddfci/1e6)*1.1;
pfcicon = -sum(real(data.cons.fci)')/1e6;
nbarmax = max(data.gene.nbar/1e19)*1.1;
ipmax   = max(data.gene.ip)*1.1/1e6;
if ~isempty(jeux1)
  ipmax   = max(ipmax,max(jeux1.data.gene.ip)*1.1/1e6);
  phybmax = max(phybmax,max(jeux1.data.gene.paddhyb/1e6)*1.1);
end
	
h1=subplot(5,1,1);
vgraph1  = [tmin tmax 0 max(phybmax,nbarmax)];
%plot(t,data.gene.paddhyb/1e6,t,data.gene.nbar/1e19,'r--',t,data.gene.ihyb/1e6,'b.-', ...
%        t,data.gene.paddfci/1e6);
if ~isempty(jeux1)
  plot(t,data.gene.paddhyb/1e6,'r',t1,jeux1.data.gene.paddhyb/1e6,'ro',...
  t1,jeux1.data.gene.nbar/1e19,'b--')
else
  plot(t,data.gene.paddhyb/1e6,'r',t,data.gene.nbar/1e19,'b--',...
       t,data.gene.paddfci/1e6,'mo' )
end
set(h1,'xticklabel','')
taille = get(gca,'fontsize');
set(gca,'fontsize',12)   
axis(vgraph1)
grid
set(h1,'xtick',tickgr);
text('units','normalized','position',[1.05 0.5],'string','a)');

h1=subplot(5,1,2);
vgraph1 = [tmin tmax 0 max(data.gene.ip(ind)/1e6)*1.1];

plot(tt,data.gene.ip(ind)/1e6,'r', ...
       tt,data.gene.ini(ind)/1e6,'b--');
set(h1,'xticklabel','')
axis(vgraph1)
set(h1,'xtick',tickgr)
taille = get(gca,'fontsize');
set(gca,'fontsize',12)
if ~isempty(jeux1)
    hold on
    plot(tt1,jeux1.data.gene.ip(ind1)/1e6,'ro')
    plot(tt1,jeux1.data.gene.ini(ind1)/1e6,'b+')
    hold off
end

set(gca,'fontsize',taille)


grid
text('units','normalized','position',[1.05 0.5],'string','b)');

h2=subplot(5,1,3);
plot(tt,medfilt1(data.gene.vsurf(ind),5),'b');
  set(h2,'xticklabel','')
%  ylabel('Vloop (V)')
  set(gca,'ylim',[0 2])
  if ~isempty(jeux1)
    hold on
    plot(tt1,jeux1.data.gene.vloop(ind1),'bo')
    hold off
  end
  if exist('jeux1','var')
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)
%    legend(gca,'CRONOS','reference',2)
    set(gca,'fontsize',taille)
  else
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)
%    legend(gca,'CRONOS',2)
    set(gca,'fontsize',taille)
  end

axis(vgraph2)
set(h2,'xtick',tickgr)
grid
text('units','normalized','position',[1.05 0.5],'string','c)');

h3=subplot(5,1,4);
%plot(tt,data.gene.betap(ind)+data.gene.li(ind)/2,'b')
if ~isempty(frsup)
  plot(tt,data.gene.betap(ind).*(1-frsup(ind)'/100),'b')
else
  plot(tt,data.gene.betap(ind),'b')
end
if ~isempty(jeux1)
  hold on
%    plot(tt1,jeux1.data.gene.betap(ind1)+jeux1.data.gene.li(ind1)/2,'bo')
  plot(tt1,jeux1.data.gene.betap(ind1),'bo')
  taille = get(gca,'fontsize');
  set(gca,'fontsize',12)  
%    legend(gca,'CRONOS','reference',2)
  set(gca,'fontsize',taille)
  hold off
else
  taille = get(gca,'fontsize');
  set(gca,'fontsize',12)
  set(gca,'fontsize',taille)
end



set(h3,'xticklabel','')
%ylabel('beta+li/2')
axis(vgraph3)
set(h3,'xtick',tickgr)
grid
text('units','normalized','position',[1.05 0.5],'string','d)');

h4=subplot(5,1,5);

disp('q au bord')
ind95 =length(param.gene.x);
plot(tt,data.prof.q(ind,ind95),'b')
  if ~isempty(jeux1)
    hold on
    plot(tt1,jeux1.data.prof.q(ind1,ind95),'bo')
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)  
%    legend(gca,'CRONOS','reference',2)
    set(gca,'fontsize',taille)
    hold off
  else
    taille = get(gca,'fontsize');
    set(gca,'fontsize',12)  
    set(gca,'fontsize',taille)    
  end

%ylabel('q_a')
axis(vgraph4)
set(h4,'xtick',tickgr)

grid
text('units','normalized','position',[1.05 0.5],'string','e)');

xlabel('time (s)')
set(gcf,'color',[1 1 1])
set(gcf,'ResizeFcn','')

