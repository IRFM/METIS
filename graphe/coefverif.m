% script de plot des courants
x=param.gene.x;
t=data.gene.temps;

% calcul des coefficients a partir des flux
keeff      =  data.prof.flux.keeff; 
kieff      =  data.prof.flux.kieff; 
chieeff    =  keeff ./ data.prof.ne;
chiieff    =  kieff ./ data.prof.ni;
deeff      =  data.prof.flux.deeff; 
kecal      =  data.prof.flux.kean; 
kical      =  data.prof.flux.kian; 
chiecal    =  kecal ./ data.prof.ne;
chiical    =  kical ./ data.prof.ni;
decal      =  data.prof.flux.dean; 
chiean     =  data.coef.ee ./ data.prof.ne;
chiian     =  data.coef.ii ./ data.prof.ni;
dean       =  data.coef.nn;
chieneo     =  data.neo.coef.ee ./ data.prof.ne;
chiineo     =  data.neo.coef.ii ./ data.prof.ni;
deneo       =  data.neo.coef.nn;


figure('color',[1 1 1 ],'defaultaxesyscale','log','defaultaxesfontsize',12);

subplot(2,2,1)
plotprof(gca,t,x,chieeff,'linestyle','-','color',[1 0 0]);
plotprof(gca,t,x,chiean,'linestyle','-','color',[1 0 1]);
plotprof(gca,t,x,chiecal,'linestyle','-','color',[0 1 1]);
plotprof(gca,t,x,chieneo,'linestyle','-','color',[0 0 0]);
title('Chie')
xlabel('sqrt(Phi/pi/B0)');
ylabel('Chie (m^2*s^-^1');
axis([0 1 1e-4 1e2])
%legend('interpretatif eff','modele','interpretatif an','Nclass');

subplot(2,2,2)
plotprof(gca,t,x,chiieff,'linestyle','-','color',[1 0 0]);
plotprof(gca,t,x,chiian,'linestyle','-','color',[1 0 1]);
plotprof(gca,t,x,chiical,'linestyle','-','color',[0 1 1]);
plotprof(gca,t,x,chiineo,'linestyle','-','color',[0 0 0]);
title('Chii')
xlabel('sqrt(Phi/pi/B0)');
ylabel('Chii (m^2*s^-^1');
axis([0 1 1e-4 1e2])
%legend('interpretatif eff','modele','interpretatif an','Nclass');

subplot(2,2,3)
plotprof(gca,t,x,deeff,'linestyle','-','color',[1 0 0]);
plotprof(gca,t,x,dean,'linestyle','-','color',[1 0 1]);
plotprof(gca,t,x,decal,'linestyle','-','color',[0 1 1]);
plotprof(gca,t,x,deneo,'linestyle','-','color',[0 0 0]);
title('De')
xlabel('sqrt(Phi/pi/B0)');
ylabel('De (m^2*s^-^1');
axis([0 1 1e-4 1e2])
legend('interpretatif effectif','modele','interpretatif anormal','Nclass');

subplot(2,2,4)
plotprof(gca,t,x,abs(data.coef.ve),'linestyle','-','color',[1 0 0]);
plotprof(gca,t,x,abs(data.coef.vi),'linestyle','-','color',[1 0 1]);
plotprof(gca,t,x,abs(data.coef.vn),'linestyle','-','color',[0 1 1]);
plotprof(gca,t,x,abs(data.neo.coef.ve),'linestyle',':','color',[1 0 0]);
plotprof(gca,t,x,abs(data.neo.coef.vi),'linestyle',':','color',[1 0 1]);
plotprof(gca,t,x,abs(data.neo.coef.vn),'linestyle',':','color',[0 1 1]);
title('V convection')
xlabel('sqrt(Phi/pi/B0)');
ylabel('V (m*s^-^1');
legend('ve an','vi an','vn an','ve neo','vi neo','vn neo');


uicontrol('style','check','units','normalized','position',[0 0.95 0.05 0.05],'callback', ...
          'zlinlog','value',1);

uicontrol('style','push','units','normalized','position',[0 0.75 0.05 0.05],'callback', ...
          'axis([0,1,0,10])','value',1);

zoom yon
