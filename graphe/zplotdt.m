% plot l'histogramme des dt
function zplotdt(file)

gene = evalin('base','data.gene');
param = evalin('base','param');

if nargin ==0
	file =param.gene.file;
elseif isempty(file)
	file =param.gene.file;
end

[dt,tmin,tmax,res,restot] = zgetdt(file);


h=findobj(0,'type','figure','tag','zplotdt');
if isempty(h)
	h=figure('tag','zplotdt');
else
	figure(h);
end
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	   'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',3,'name',file)
	   
h1=subplot(2,3,1);
plot(gene.temps,gene.conv,'o')
nbconv = sum(gene.conv(isfinite(gene.conv)));
title(sprintf('Convergence (%d/dt)',ceil(nbconv./length(dt))))
ylabel('Nombre de boucles')
xlabel('temps')
set(gca,'xlim',[min(tmin),max(tmax)])
subplot(2,3,4)
hist(gene.conv,20)
ylabel('frequences')
xlabel('Nombre de boucles')
title(sprintf('moy = %4.2g, std = %4.2g',mean(gene.conv(isfinite(gene.conv))),std(gene.conv(isfinite(gene.conv)))))
h2=subplot(2,3,2);
plot(gene.temps,gene.nbsplit,'o')
indok =find(isfinite(gene.nbsplit))
title(sprintf('Decoupe temporelle (%d/temps)',ceil(sum(gene.nbsplit(indok))/length(indok))))
ylabel('Nombre d''intervalles ')
set(gca,'xlim',[min(tmin),max(tmax)])
xlabel('temps')
subplot(2,3,5)
hist(gene.nbsplit,20)
ylabel('frequences')
xlabel('Nombre d''intervalles')
title(sprintf('moy = %4.2g, std = %4.2g',mean(gene.nbsplit(isfinite(gene.nbsplit))),std(gene.nbsplit(isfinite(gene.nbsplit)))))
h3=subplot(2,3,3);
semilogy(tmin,dt,'o')
title(sprintf('pas de temps interne (%s)',file))
ylabel('dt (s)')
xlabel('temps')
subplot(2,3,6)
hist(log(dt)/log(10),20)
ylabel('frequences')
xlabel('log(dt) ')
title(sprintf('moy = %4.2g, std = %4.2g (ms)',mean(dt)*1000,std(dt)*1000))

drawnow
xl=[fix(min(tmin)),ceil(max(tmax))];
set([h1,h2,h3],'xlim',xl);


h=findobj(0,'type','figure','tag','zplotdt2');
if isempty(h)
	h=figure('tag','zplotdt2');
else
	figure(h);
end
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	   'defaultlinelinewidth',3,'color',[1 1 1],'defaultlinemarkersize',3,'name',file)
nbt =length(gene.temps);	   
hp1=subplot(2,2,1);
plot(gene.temps,gene.flops)
ylabel('Flops')
title(sprintf('%d MFlops/temps',ceil((max(gene.flops)-min(gene.flops))/length(indok)/1e6)))
hp4=subplot(2,2,4);

cpu =real(gene.cputime);
ind = find(isfinite(cpu));
dcpudt = zdxdt(cpu(ind),gene.temps(ind));
cpu = cumtrapz(gene.temps(ind),dcpudt .* (dcpudt >=0));

plot(gene.temps(ind),real(cpu))
ylabel('CPU (s)')
title(sprintf('%d s/temps (#%d)',ceil((max(real(cpu))-min(real(cpu)))/length(indok)),length(indok)))
xlabel('temps (s)')
hp3=subplot(2,2,3);
plot(gene.temps,real(gene.memory))
ylabel('Memoire (Mo)')
xlabel('temps (s)')
title(sprintf('%d Mo/temps (#%d)',ceil(max(real(gene.memory))/length(indok)),length(indok)))
hp2=subplot(2,2,2);
dd =real(gene.datation);
ind = find(isfinite(dd));
ddddt = zdxdt(dd(ind),gene.temps(ind));
dd = cumtrapz(gene.temps(ind),ddddt .* (ddddt >=0));
plot(gene.temps(ind),dd)
ylabel('temps machine (s)')
title(file,'interpreter','none');

set([hp1,hp2,hp3,hp4],'xlim',xl);
