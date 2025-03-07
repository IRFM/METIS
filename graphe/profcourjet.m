% script de courant
load /usr/drfc/cgc/matlab5/tcron/JET/51782/efit51782.mat
tt = input('Le temps ? ');
ind =iround(data.gene.temps,tt);
indc = iround(tefit,tt);

x    = param.gene.x;
johm = data.prof.jmoy - data.source.idn.j - data.neo.jboot;
figure('color',[1 1 1],'defaultaxesfontsize',18,'defaultaxesfontname','times','defaultlinelinewidth',3);
plot(xefit,-Jx(indc,:)/1e6,'ro',x,data.equi.jmoy(ind,:)/1e6,'k',x,data.source.idn.j(ind,:)/1e6,'b');
hold on
plot(x,data.neo.jboot(ind,:)/1e6,'color',[0.7,0.5,0.7]);
plot(x,johm(ind,:)/1e6,'k:');
plot([0 1],[0 0],'color',[0 0.7 0]);
legend('J efit','J equi.','J IdN','J boot','J \Omega');
xlabel('sqrt(\Phi / \pi / B_0) normalise')
ylabel('MA*m^-^2')
st =sprintf('choc %s #%d, t = %g s, nbrho = %d',param.from.machine, ...   
             fix(param.from.shot.num),tt,param.gene.nbrho);
title(st);
