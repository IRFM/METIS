% script de plot de la polarimetrie sur ts

t =data.gene.temps;
nl = post.polar.nl;
nlmes = post.polar.nlmes;
af = post.polar.af.*180./pi;
afmes = post.polar.afmes.*180./pi;

% correction de l'ecart d'alignement vertical
ind = find(isfinite( prod(af' .* afmes')'));
daf = mean( af(ind,:) -afmes(ind,:));
ind = find(isfinite( prod(nl' .* nlmes')'));
dnl = std(nl(ind,:)-nlmes(ind,:))./mean(nlmes(ind,:));

figure('color',[1 1 1])
cc=get(gca,'colororder');
plot(t,af-ones(size(t))*daf);
hold on
set(gca,'colororder',cc);
plot(t,afmes,'o');
st =sprintf('choc %s #%d, t = %g:%g:%g, nbrho = %d',param.from.machine, ...   
             fix(param.from.shot.num),param.gene.tdeb,mean(diff(data.gene.temps)),param.gene.tfin, ...
             param.gene.nbrho);
title(st);
grid
xlabel('temps (s)');
ylabel('angle Faraday (degres)');
legend(num2str(daf(1)),num2str(daf(2)),num2str(daf(3)),num2str(daf(4)),num2str(daf(5)));
gtext(mat2str(dnl));
