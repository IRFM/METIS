% script de test de l'hypothese de l'etat critique dans le plasma

t        = data.gene.temps;
x        = param.gene.x;    
jcrossb_  = data.prof.bpol .* data.prof.jphi - data.prof.bphi .* data.prof.jpol;
gradptot = pdederive(x,data.prof.ptot,0,2,2,1) ./ data.equi.rhomax .* data.equi.grho;

h=figure;
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',3,'color',[1 1 1]);
subplot(2,2,2)
plotprof(gca,t,x,jcrossb_,'color','r','linestyle','-')
plotprof(gca,t,x,gradptot,'color','b','linestyle','o')
legend('JxB','\NablaPtot');
ylabel('Pa/m')
xlabel('x (su)')
