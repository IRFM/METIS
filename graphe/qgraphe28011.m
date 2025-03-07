% script de plot du profil d q
t=data.gene.temps;
x=param.gene.x;
ind =iround(t,[4 7 14 14.5 15]);
figure('color',[1 1 1]);
cc=get(gca,'colororder');
plot(x,data.prof.q(ind,:));
legend('4','7','14','14.5','15');
hold on
set(gca,'colororder',cc);
plot(x,data.equi.q(ind,:),':');
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
figure('color',[1 1 1]);
cc=get(gca,'colororder');
plot(x,data.prof.jmoy(ind,:));
legend('4','7','14','14.5','15');
hold on
set(gca,'colororder',cc);
plot(x,data.equi.jmoy(ind,:),':');
hold off
xlabel('sqrt(phi/pi/b0) (m)')
ylabel('Jmoy (A*m^-^2)')
st =sprintf('choc %s #%d, t = %g:%g:%g, nbrho = %d',param.from.machine, ...   
             fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
title(st);
