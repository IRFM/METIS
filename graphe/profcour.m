% script de courant
coup = cgcgettrait(param.from.shot.num,'tcoupol');
tt = input('Le temps ? ');
ind =iround(data.gene.temps,tt);
indc = iround(coup.times,tt);
x    = param.gene.x;
rho  = coup.rhofit;
johm = data.prof.jmoy - data.source.hyb.j - data.neo.jboot;
figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname','times','defaultlinelinewidth',3);
plot(rho,coup.jmoy(indc,:),'ro',x,data.equi.jmoy(ind,:)/1e6,'k',x,data.source.hyb.j(ind,:)/1e6,'b');
hold on
plot(x,data.neo.jboot(ind,:)/1e6,'color',[0.7,0.5,0.7]);
plot(x,johm(ind,:)/1e6,'k:');
legend('J exp.','J equi.','J LH','J boot','J \Omega');
xlabel('sqrt(\Phi / \pi / B_0) normalise')
ylabel('MA*m^-^2')
st =sprintf('choc %s #%d, t = %g, nbrho = %d',param.from.machine, ...   
             fix(param.from.shot.num),tt,param.gene.nbrho);
title(st);
