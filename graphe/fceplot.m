% script de plot du depot FCE
h1 =findobj(0,'type','figure','tag','fceplot');
if isempty(h1)
	  h1=figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname', ...
  	          'times','defaultlinelinewidth',3,'tag','fceplot','name','ECRH');
else
	  figure(h1);
end
clf

t= data.gene.temps;
x = param.gene.x;

subplot(3,1,1)
plot(t,data.gene.ifce./1e3,'b')
st =sprintf('choc %s #%d',param.from.machine, param.from.shot.num);
title(st);
	
ylabel('I_E_C_R_H (kA)')

subplot(3,1,2)
plot(t,data.gene.paddfce,'r')
ylabel('P_E_C_R_H (MW)')
xlabel('temps (s)')

subplot(3,2,5)
zplotprof(gca,t,x,data.source.fce.el,'color','r');
ylabel('P_E_C_R_H (W m^-^3)')
xlabel('r/a')
set(gca,'ylim',[0,inf]);
subplot(3,2,6)
zplotprof(gca,t,x,data.source.fce.j,'color','b');
ylabel('Jl_E_C_R_H (A m^-^2)')
xlabel('r/a')
	
h1 =findobj(0,'type','figure','tag','fceplot2');
if isempty(h1)
	  h1=figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname', ...
  	          'times','defaultlinelinewidth',3,'tag','fceplot2','name','Source ECRH');
else
	  figure(h1);
end
clf
zplotprof(gca,t,x,data.source.fce.el,'color','b');
ylabel('P_E_C_R_H (W m^-^3)')
xlabel('r/a')
set(gca,'ylim',[0,inf]);
