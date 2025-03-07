x = param.gene.x;
t = data.gene.temps;
vprd1     = rpdederive(x,data.equi.vpr,2,2,2,1);
ped1     = rpdederive(x,data.prof.pe,2,2,2,1);
dpedt     = zdxdt(data.prof.pe,t);
dvprdt     = zdxdt(data.equi.vpr,t);
sel       =data.source.totale.el;
inte = data.equi.vpr .* data.source.totale.el - (3/2).* data.equi.vpr .* data.prof.dpedt -  ...
       (5/2)  .* dvprdt .* data.prof.pe + (ones(size(t))*x) ./ (data.equi.rhomax*ones(size(x))) .* (data.equi.drhomaxdt*ones(size(x))) .* ...
       ((3/2) .* data.equi.vpr .* ped1 + (5/2) .* vprd1 .* data.prof.pe);



qe = (data.equi.rhomax*ones(size(x))) ./ data.equi.vpr ./ data.equi.grho2 .* cumtrapz(x,inte,2);
qe = zcentre(qe);    

v0   = vprd1(:,3)*ones(size(x));
lqe  =  data.source.totale.el - (3/2).* data.prof.dpedt;
inte2  = cumtrapz(x,lqe,2);
xx     =  ones(size(t)) * x;
qec  = (data.equi.rhomax*ones(size(x))) ./ data.equi.vpr ./ data.equi.grho2 .* v0 .* ( xx .* inte2  - cumtrapz(x,inte2,2));





%intep     = zcentre(rpdederive(x,zcentre(qe./((data.equi.rhomax*ones(size(x))) ./ data.equi.vpr ./ data.equi.grho2)),2,2,2,1));
%xx        = linspace(0,1,1001);
%intex     = interp1(x',inte',xx','spline')'; 
%qex = cumtrapz(xx,intex,2);
gte = rpdederive(x,data.prof.te,2,2,2,1) ./ (data.equi.rhomax*ones(size(x))) ;
kieff = - qe ./ gte ./ param.phys.e; 
chieeff = kieff ./ data.prof.ne;
% calcul du flux avec le coef de transport
qeee = -(data.coef.ee+data.neo.coef.ee) .* gte .* param.phys.e - (data.coef.ve + data.neo.coef.ve) .* data.prof.pe; ;
h1 =findobj(0,'type','figure','tag','flux1');
if isempty(h1)
    h1 = figure('tag','flux1');
else
    figure(h1);
    clf
end

zplotprof(gca,t,x,qe);
zplotprof(gca,t,x,qec,'color','m');
zplotprof(gca,t,x,qeee,'color','r');



h2 = findobj(0,'type','figure','tag','flux2');
if isempty(h2)
    h2 = figure('tag','flux2');
else
    figure(h2);
    clf
end
zplotprof(gca,t,x,data.source.totale.el .* data.equi.vpr,'color','b');
zplotprof(gca,t,x,-(3/2) .* data.prof.dpedt .* data.equi.vpr,'color','r');
zplotprof(gca,t,x,-(5/2) .* data.prof.pe .* data.equi.dvprdt,'color','c');
dia = (ones(size(t))*x) ./ (data.equi.rhomax*ones(size(x))) .* (data.equi.drhomaxdt*ones(size(x)));
zplotprof(gca,t,x,-(3/2) .* data.prof.ped1 .* data.equi.vpr.*dia,'color','m');
zplotprof(gca,t,x,-(5/2) .* data.prof.pe .* vprd1 .* dia,'color','k');




gqe = qeee ./ (data.equi.rhomax*ones(size(x))) .* data.equi.vpr .* data.equi.grho2;
gqe = rpdederive(x,gqe,2,2,2,1);
h3 = findobj(0,'type','figure','tag','flux3');
if isempty(h3)
    h3 = figure('tag','flux3');
else
    figure(h3);
    clf
end
zplotprof(gca,t,x,gqe);
zplotprof(gca,t,x,inte,'color','r');
