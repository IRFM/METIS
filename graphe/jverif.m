% script de plot des courants
x=param.gene.x;
t=data.gene.temps;
figure('color',[1 1 1 ]);
je = data.prof.epar ./ data.coef.eta +data.source.totale.j;
jo = data.prof.ej ./  data.prof.epar;
plotprof(gca,t,x,data.prof.jmoy,'linestyle','-','color',[1 0 0]);
plotprof(gca,t,x,data.prof.jeff,'linestyle','-.','color',[1 0 1]);
plotprof(gca,t,x,data.prof.jphi,'linestyle',':','color',[0 1 1]);
plotprof(gca,t,x,je,'linestyle','o','color',[0 0 0]);
plotprof(gca,t,x,jo,'linestyle','+','color',[0 0 0]);
title('Courants')
xlabel('sqrt(Phi/pi/B0)');
ylabel('J (A/m^2)');
legend('Jmoy','Jeff','Jphi','E/eta+Jni','<E.J>/E');

