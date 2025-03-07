% script pour le calcul des erreurs pour l'article su la MHD avec DDS monstre
zuiload('/usr/drfc/zineb/data/ts28207_maget.mat.gz')
ind = find(data.gene.temps >= 18 & data.gene.temps <= 19);
epfci = mean(post.ripple.frlost(ind))
dte = std(data.prof.te(ind,:),0,1);
ete = dte ./ mean(data.prof.te(ind,:));
gte = abs(pdederive(param.gene.x,data.prof.te,0,2,2,1));
dgte = sqrt(std(gte(ind,:),0,1) .^ 2 + abs(pdederive(param.gene.x,dte,0,2,2,1)) .^2);
egte = dgte ./ mean(gte(ind,:));
indx = find(param.gene.x >= 0.4 & param.gene.x <= 0.6);
errdelta = epfci + mean(egte(:,indx)) + 5/2 * mean(ete(:,indx))
errte = mean(egte(:,indx)) + 5/2 * mean(ete(:,indx))
mean(egte(:,indx))
mean(ete(:,indx))

% verification des sources
rhomax = mean(data.equi.rhomax(ind),1);
vpr    = mean(data.equi.vpr(ind,:),1);
sefci  = mean(data.source.fci.el(ind,:),1);
sqefci = rhomax .* cumtrapz(param.gene.x,sefci .*vpr);
sei   = mean(data.source.qei(ind,:),1);
sqei  = rhomax .* cumtrapz(param.gene.x,sei .*vpr);
x = param.gene.x;

figure
plot(x,1-sqefci./sqefci(end),'b',x,abs(sqei)./sqefci,'r');
axis([0.4 0.6 0 0.2])
mean(1-sqefci(indx)./sqefci(end))
mean(abs(sqei(indx))./sqefci(indx))
