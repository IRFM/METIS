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

valip = data.gene.ip./data.gene.ini;
indip = find(valip < 1);
ind   = param.gene.kmin:min(param.gene.k,length(t));
corip = iround(ind,indip);
if ~isempty(corip)
  ind(corip) = [];
end
tt   = t(ind);
%inde = 1:pas:length(t);
%t    = t(inde);

if exist('jeux1','var')
  if ~isempty(jeux1)
    ind1=jeux1.param.gene.kmin:min(jeux1.param.gene.k,length(t1));
    tt1=t1(ind1);
  else
    clear jeux1
  end
end

h = findobj(0,'type','figure','tag','compare');
if isempty(h)
       h=figure('tag','compare');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',4)

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
tmin = param.gene.tdeb;
vgraph1 = [tmin tmax 0 max(data.gene.ip(ind)/1e6)*1.1];
vgraph2 = [tmin tmax 0 2];
vgraph3 = [tmin tmax 0.5 2];
ind95   =length(param.gene.x);

vgraph4 = [tmin tmax 2 max(data.prof.q(ind,ind95))*1.2];
tickgr  = linspace(0,round(tt(imax)),6);
   phybmax  = max(data.gene.paddhyb/1e6)*1.1;
   pfcimax  = max(data.gene.paddfci/1e6)*1.1;
   pfcicon  = -sum(real(data.cons.fci)')/1e6;
   nbarmax  = max(data.gene.nbar/1e19)*1.1;
	
h1=subplot(4,1,1);
vgraph1  = [tmin tmax 0 max([phybmax pfcicon*1.1 nbarmax])];
plot(t,data.gene.paddhyb/1e6,t,data.gene.nbar/1e19, ...
        t,abs(data.cons.fci)/1e6);
set(h1,'xticklabel','')
taille = get(gca,'fontsize');
set(gca,'fontsize',12)   
set(h1,'xtick',tickgr)
axis(vgraph1)
grid


h1=subplot(4,1,2);
vgraph1 = [tmin tmax 0 max(data.gene.ip(ind)/1e6)*1.1];

plot(tt,data.gene.ip(ind)/1e6,'r',tt,data.gene.iboot(ind)/1e6,'.m', ...
       tt,data.gene.ini(ind)/1e6,'k.-');
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
end

set(gca,'fontsize',taille)


grid

h2=subplot(4,1,3);
plot(tt,medfilt1(data.gene.vsurf(ind),5),'k-');
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

h3=subplot(4,1,4);
plot(tt,data.gene.betap(ind)+data.gene.li(ind)/2,'b.-')
  if exist('jeux1')
    hold on
    plot(tt1,jeux1.data.gene.betap(ind1)+jeux1.data.gene.li(ind1)/2,'r:')
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

%set(h3,'xticklabel','')
ylabel('beta+li/2')
axis(vgraph3)
set(h3,'xtick',tickgr)
grid

xlabel('time (s)')
set(gcf,'color',[1 1 1])
set(gcf,'ResizeFcn','')

