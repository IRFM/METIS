% script de plot du scenario
h = findobj(0,'type','figure','tag','tsscenario');
if isempty(h)
       h=figure('tag','tsscenario');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',18,'defaultaxesfontweight','bold','defaultaxesfontname','helvetica', ...
	'defaultlinelinewidth',2,'color',[1 1 1])
x=param.gene.x;
t=data.gene.temps;
if strcmp(param.gene.filetype,'source')
   ip   = data.cons.ip./1e6;
   nbar = trapz(x,data.prof.ne,2)./1e19;  
   idn  = sum(real(data.cons.idn),2)./1e6;
   fci  = sum(abs(data.cons.fci),2)./1e6;
   hyb  = sum(abs(data.cons.hyb),2)./1e6;
   fce  = sum(abs(data.cons.fce),2)./1e6;
else
   ip   = data.gene.ip./1e6;
   nbar = data.gene.nbar./1e19; 
   idn  = data.gene.paddidn./1e6;
   fci  = data.gene.paddfci./1e6;
   hyb  = data.gene.paddhyb./1e6;
   fce  = data.gene.paddfce./1e6;
end
plot(t,ip,'g',t,nbar,'b',t,hyb,'r',t,fci,'m',t,fce,'k',t,idn,'c');
st =sprintf('%s %d',param.from.machine, param.from.shot.num);
title(st);
xlabel('Time (s)')
legend('Ip (MA)','Nbar (1e19 m^-^3)','LH (MW)','ICRH (MW)','ECRH (MW)','NBI (MW)')
hm = findobj(gcf,'type','line','color',[0 1 1]);
set(hm,'color',[0 0.72 0.72]);
hm = findobj(gcf,'type','line','color',[1 1 0]);
set(hm,'color','k');
hm = findobj(gcf,'type','line','color',[0 1 0]);
set(hm,'color',[0 0.72 0]);
