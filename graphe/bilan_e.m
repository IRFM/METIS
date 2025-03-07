% script de plot du bilan de particules
t=data.gene.temps;
gene=data.gene;
sortant = data.prof.flux.sortant;
bord=data.bord;
ind= param.gene.kmin:param.gene.k;

figure('color',[1 1 1],'defaultaxesfontsize',12);
subplot(2,2,1 )
plot(t(ind),gene.netot(ind),'o')
ylabel('electrons');
set(gca,'ylim',[0 max(gene.netot(ind)).*1.2]);
title('evolution du contenu de la decharge')
subplot(2,2,2)
plot(t(ind),gene.bilan_e(ind)./gene.netot(ind));
ylabel('s^-1');
title('bilan instantane')
subplot(2,2,3)
plot(t(ind),gene.zeffm(ind),t(ind),gene.zeff0d(ind),'o');
title('<zeff>')
legend('zeffm','zeff0dt');
xlabel('temps (s)')
ylabel('Zeff');
subplot(2,2,4)
semilogy(t(ind),gene.nions(ind,:),t(ind),gene.evolution(ind,:),'o');
title('composition ');
xlabel('temps(s)');
ylabel('atomes')

