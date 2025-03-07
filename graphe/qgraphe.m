% script de plot du profil d q
t=data.gene.temps;
x=param.gene.x;
figure;
plotprof(gca,t,x,data.prof.q,'color',[1 0 0]);
plotprof(gca,t,x,data.equi.q,'color',[1 0 1]);
hold on
plot([0 1],[1 1],'g');
plot([0 1],[1.5 1.5],'c');
plot([0 1],[2 2],'b');
hold off
xlabel('sqrt(phi/pi/b0) (m)')
ylabel('q')
st =sprintf('choc %s #%d, t = %g:%g:%g, nbrho = %d',param.from.machine, ...   
             fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
title(st);
