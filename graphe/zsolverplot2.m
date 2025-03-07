    
% plot de control
h=findobj(0,'type','figure','tag','fp');
if isempty(h)
	h=figure('color',[ 1 1 1],'defaultaxesfontsize',12,'tag','fp');
	hc=uicontrol(h,'style','radio','tag','stop');
else
	figure(h);
	hc = findobj(h,'type','uicontrol');
end

subplot(2,2,1)
ii =pe;
if premier_passage == 1
   hold off
   plot(gene.x,fp1(:,ii)-F(:,ii),':r',gene.x,fp2(:,ii)-F(:,ii),':g',gene.x,fp3(:,ii)-F(:,ii),':b', ...
       gene.x,fpini(:,ii)-F(:,ii),':m')
   title('r -> adia, g -> lin, b -> expl., m -> ini, c -> fin')
   xlabel(['dt = ',num2str(dt)]);  
   ind = 1:iround(gene.x,1);
   ymax = max([max(fp1(ind,ii)-F(ind,ii)),max(fp2(ind,ii)-F(ind,ii)), ...
           max(fpini(ind,ii)-F(ind,ii)),max(fp(ind,ii)-F(ind,ii))]);
   ymin = min([min(fp1(ind,ii)-F(ind,ii)),min(fp2(ind,ii)-F(ind,ii)), ...
           min(fpini(ind,ii)-F(ind,ii)),min(fp(ind,ii)-F(ind,ii))]);
   set(gca,'ylim',[ymin,ymax])
   hold on
end
hold on
plot(gene.x,fp(:,ii)-F(:,ii));
           
subplot(2,2,2)
if premier_passage == 1
   hold off
%   plot(gene.x,datak.prof.gte,'r',gene.x,datakp1.prof.gte,'c')
   plot(gene.x,datak.source.qei,'r',gene.x,datakp1.source.qei,'c')
   title('r -> ini, c -> 1er k+1, b -> autre')
   ylabel('Qei')
   hold on
end
hold on
plot(gene.x,datakp1.source.qei,'b')
%plot(gene.x,datakp1.prof.gte)

subplot(2,2,3)
if premier_passage == 1
   hold off
%   plot(gene.x,datak.prof.lte,'or',gene.x,datakp1.prof.lte,'+c')
   plot(gene.x,datak.source.totale.el,'r',gene.x,datakp1.source.totale.el,'c')
   title('r -> ini, c -> 1er k+1, b -> autre')
   ylabel('Pel ')
   hold on
end
hold on
%plot(gene.x,datakp1.prof.lte)
plot(gene.x,datakp1.source.totale.el,'b')

subplot(2,2,4)
semilogy(eee,'or');
title(['dt = ', num2str(gene.dt)]);
xlabel(nb)
while(get(hc,'value') == 1)
   pause(1)
end
