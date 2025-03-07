% script de test de l'hypothese de l'etat critique dans le plasma

t        = data.gene.temps;
x        = param.gene.x;    
jcrossb  = data.prof.bpol .* data.prof.jphi - data.prof.bphi .* data.prof.jpol;
gradptot = pdederive(x,data.prof.ptot,0,2,2,1) ./ (data.equi.rhomax*ones(size(x))) .* data.equi.grho;

h=figure;
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',3,'color',[1 1 1]);
subplot(2,2,1)
plotprof(gca,t,x,jcrossb,'color','r','linestyle','-');
plotprof(gca,t,x,-gradptot,'color','b','linestyle','o');
legend('JxB','grad(Ptot)');
ylabel('Pa/m')
xlabel('x (su)')
subplot(2,2,2)
plotprof(gca,t,x,data.prof.jphi,'color','r','linestyle','-');
plotprof(gca,t,x,data.prof.jpol,'color','b','linestyle','-');
legend('J\phi','Jpol');
ylabel('A*m^-^2')
xlabel('x (su)')
subplot(2,2,3);
plotprof(gca,t,x,data.prof.bphi,'color','r','linestyle','-');
plotprof(gca,t,x,data.prof.bpol,'color','b','linestyle','-');
legend('B\phi','Bpol');
ylabel('A*m^-^2')
xlabel('x (su)')
