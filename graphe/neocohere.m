% plot la coherence entre les flux et les coef de transport
rm      = data.equi.rhomax * ones(1,param.gene.nbrho);
%neoge   = - data.neo.coef.nn .* data.prof.ned1 ./ rm  -  data.neo.coef.vn .* data.prof.ne;
neoge   = - data.neo.coef.nn .* data.prof.gne ./ data.equi.grho  -  data.neo.coef.vn .* data.prof.ne;
neoqe   = - data.neo.coef.ee .* data.prof.gte .*  param.phys.e ./ data.equi.grho    -   ...
            data.neo.coef.ve .* data.prof.pe;
neoqi   = - data.neo.coef.ii  .* data.prof.gti .*  param.phys.e ./ data.equi.grho  -   ...
            data.neo.coef.vi .* data.prof.pion;


t = data.gene.temps;
x = param.gene.x;            
vpr = data.equi.vpr;

%[ge,qe,ge_b,qe_b,ge_p,qe_p]= neocertif(data,param);

h=figure;
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
     'defaultlinelinewidth',3,'color',[1 1 1])
subplot(2,2,1)
plotprof(gca,t,x,neoqe,'linestyle','-','color',[1 0 0]);
plotprof(gca,t,x,data.neo.flux.qe,'linestyle','-','color',[0 1 1]);
%plotprof(gca,t,x,qe,'linestyle','-','color',[0 1 0]);
title('electron neoclassical heat flux' )
legend('from profiles','from coeff')
ylabel('Qe (W*m^-^2)p')
subplot(2,2,2)
plotprof(gca,t,x,neoqi,'linestyle','-','color',[1 0 0]);
plotprof(gca,t,x,data.neo.flux.qion,'linestyle','-','color',[0 1 1]);
title('ion neoclassical heat flux' )
legend('from profiles','from coeff')
ylabel('Qi (W*m^-^2)p')
subplot(2,2,3)
plotprof(gca,t,x,neoge,'linestyle','-','color',[1 0 0]);
plotprof(gca,t,x,data.neo.flux.ne,'linestyle','-','color',[0 1 1]);
%plotprof(gca,t,x,-ge,'linestyle','-','color',[0 1 0]);
title('neoclassical electron flux' )
legend('from profiles','from coeff')
ylabel('Ge (m^-^2*s^-^1)')
xlabel('x')
subplot(2,2,4);
%plotprof(gca,t,x,data.prof.pion,'linestyle','o','color',[1 0 0]);
%pion2 = sum(data.impur.impur,3) .* data.prof.ti .* param.phys.e;
%plotprof(gca,t,x,pion2,'linestyle','-','color',[0 1 1]);
plotprof(gca,t,x,data.neo.flux.ne ./neoge,'linestyle','-','color',[1 0 0])
plotprof(gca,t,x,neoge ./ data.neo.flux.ne ,'linestyle','-','color',[0 1 1])
title('electron density flux -> verification')
xlabel('x')
