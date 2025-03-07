% script de plot des flux
% creation de la figure
h = findobj(0,'type','figure','tag','flux_1');
if isempty(h)
       h=figure('tag','flux_1');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',14,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',3,'color',[1 1 1])
% script de plot des courants
x=param.gene.x;
t=data.gene.temps;

% les flux neo
neoge  = data.neo.flux.ne;
neogi  = data.neo.flux.nion;
neoqe  = data.neo.flux.qe;
neoqi  = data.neo.flux.qion;

% les flux neo reconstitue
rm       = data.equi.rhomax * ones(1,param.gene.nbrho);
rneoge   = - data.neo.coef.nn .* data.prof.gne ./ data.equi.grho  -  data.neo.coef.vn .* data.prof.ne;
rneoqe   = - data.neo.coef.ee .* data.prof.gte .*  param.phys.e ./ data.equi.grho    -   ...
            data.neo.coef.ve .* data.prof.pe;
rneoqi   = - data.neo.coef.ii  .* data.prof.gti .*  param.phys.e ./ data.equi.grho  -   ...
            data.neo.coef.vi .* data.prof.pion;


% les flux totaux
fluxge  = data.prof.flux.ge; 
fluxgi  = data.prof.flux.gi; 
fluxqe  = data.prof.flux.qe; 
fluxqi  = data.prof.flux.qi; 

% les flux totaux reconstitue
rfluxge   = - (data.coef.nn + param.cons.neomulti.nn .* data.neo.coef.nn).* data.prof.gne ./ data.equi.grho -  ...
               (data.neo.coef.vn + param.cons.neomulti.vn .* data.neo.coef.vn) .* data.prof.ne;
rfluxqe   = - (data.coef.ee + param.cons.neomulti.ee .* data.neo.coef.ee ).* data.prof.gte .*   ...
               param.phys.e ./ data.equi.grho -  ...
              (data.neo.coef.ve + param.cons.neomulti.ve .* data.neo.coef.ve) .* data.prof.pe;
rfluxqi   = - (data.coef.ii + param.cons.neomulti.ii .* data.neo.coef.ii ).* data.prof.gti .*   ...
               param.phys.e ./ data.equi.grho -  ...
              (data.neo.coef.vi + param.cons.neomulti.vi .* data.neo.coef.vi) .* data.prof.pion;


% les flux anormaux
ange  = data.prof.flux.ge - data.neo.flux.ne; 
angi  = data.prof.flux.gi - data.neo.flux.nion; 
anqe  = data.prof.flux.qe - data.neo.flux.qe; 
anqi  = data.prof.flux.qi - data.neo.flux.qion; 

% flux anormaux reconstitue
range   = - data.coef.nn .* data.prof.gne ./ data.equi.grho  -  data.coef.vn .* data.prof.ne;
ranqe   = - data.coef.ee .* data.prof.gte .*  param.phys.e ./ data.equi.grho    -   ...
            data.coef.ve .* data.prof.pe;
ranqi   = - data.coef.ii .* data.prof.gti .*  param.phys.e ./ data.equi.grho  -   ...
            data.coef.vi .* data.prof.pion;


subplot(2,2,1)
plotprof(gca,t,x,neoge,'linestyle','-','color',[1 0 0]);
plotprof(gca,t,x,rneoge,'linestyle','o','color',[1 0 0]);
plotprof(gca,t,x,fluxge,'linestyle','-','color',[1 0 1]);
plotprof(gca,t,x,rfluxge,'linestyle','o','color',[1 0 1]);
plotprof(gca,t,x,ange,'linestyle','-','color',[0 1 1]);
plotprof(gca,t,x,range,'linestyle','o','color',[0 1 1]);
title('flux d''electrons')
%xlabel('sqrt(Phi/pi/B0)');
ylabel('(m^-^2*s^-^1)');
legend('Neo','Neo recalcule','Total','Total recalcule','Anormal','Anormal recalcule');

subplot(2,2,2)
plotprof(gca,t,x,neogi,'linestyle','-','color',[1 0 0]);
plotprof(gca,t,x,fluxgi,'linestyle','-','color',[1 0 1]);
plotprof(gca,t,x,angi,'linestyle','-','color',[0 1 1]);
title('flux d''ions (somme des especes)')
%xlabel('sqrt(Phi/pi/B0)');
ylabel('(m^-^2*s^-^1)');
%legend('Neo','Total','Anormal');

subplot(2,2,3)
plotprof(gca,t,x,neoqe,'linestyle','-','color',[1 0 0]);
plotprof(gca,t,x,rneoqe,'linestyle','o','color',[1 0 0]);
plotprof(gca,t,x,fluxqe,'linestyle','-','color',[1 0 1]);
plotprof(gca,t,x,rfluxqe,'linestyle','o','color',[1 0 1]);
plotprof(gca,t,x,anqe,'linestyle','-','color',[0 1 1]);
plotprof(gca,t,x,ranqe,'linestyle','o','color',[0 1 1]);
title('flux de chaleur electronique')
xlabel('sqrt(Phi/pi/B0)');
ylabel('(W*m^-^2*s^-^1)');
%legend('Neo','Neo recalcule','Total','Total recalcule','Anormal','Anormal recalcule');

subplot(2,2,4)
plotprof(gca,t,x,neoqi,'linestyle','-','color',[1 0 0]);
plotprof(gca,t,x,rneoqi,'linestyle','o','color',[1 0 0]);
plotprof(gca,t,x,fluxqi,'linestyle','-','color',[1 0 1]);
plotprof(gca,t,x,rfluxqi,'linestyle','o','color',[1 0 1]);
plotprof(gca,t,x,anqi,'linestyle','-','color',[0 1 1]);
plotprof(gca,t,x,ranqi,'linestyle','o','color',[0 1 1]);
title('flux de chaleur ionique')
xlabel('sqrt(Phi/pi/B0)');
ylabel('(W*m^-^2*s^-^1)');
%legend('Neo','Neo recalcule','Total','Total recalcule','Anormal','Anormal recalcule');

