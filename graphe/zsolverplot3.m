    
% plot de control
h=findobj(0,'type','figure','tag','fp');
if isempty(h)
	h=figure('color',[ 1 1 1],'defaultaxesfontsize',12,'tag','fp');
	hc=uicontrol(h,'style','radio','tag','stop');
else
	figure(h);
	hc = findobj(h,'type','uicontrol');
end
ii = pe;
jj = pion;

subplot(3,2,1)
if premier_passage == 1
   hold off
   plot(gene.x,(fpini(:,ii)-F(:,ii)) ./ std(F(:,ii)),'r');
   xlabel(['dt = ',num2str(dt)]);  
   hold on
   axis([0 1 -5e-3 5e-3]);
end
plot(gene.x,(fp(:,ii)-fpini(:,ii)) ./ std(F(:,ii)));
subplot(3,2,4)
if premier_passage == 1
   hold off
   plot(gene.x,(fpini(:,jj)-F(:,jj)) ./ std(F(:,jj)),'r');
   xlabel(['dt = ',num2str(dt)]);  
   hold on
   axis([0 1 -5e-3 5e-3]);
end
plot(gene.x,(fp(:,jj)-fpini(:,jj)) ./ std(F(:,jj)));
          
subplot(3,2,3)
if premier_passage == 1
   hold off
   plot(gene.x,datak.prof.pe,'r',gene.x,datak.prof.pion,'c')
   hold on
end
plot(gene.x,datakp1.prof.pe,'m',gene.x,datakp1.prof.pion,'b')


subplot(3,2,2)
if premier_passage == 1
   hold off
   plot(gene.x,datak.coef.ee,'r',gene.x,datak.neo.coef.ee,'b',gene.x,datak.neo.coef.ee+datak.coef.ee,'g')
   title('r -> An(k) , b -> Neo(k), g -> sum')
   ylabel('Ke ')
end
hold on
plot(gene.x,datakp1.coef.ee,'m',gene.x,datakp1.neo.coef.ee,'c',gene.x,datakp1.neo.coef.ee+datakp1.coef.ee,'k')

subplot(3,2,5)
if premier_passage == 1
   hold off
   plot(gene.x,datak.coef.ii,'r',gene.x,datak.neo.coef.ii,'b',gene.x,datak.neo.coef.ii+datak.coef.ii,'g')
   title('r -> An(k) , b -> Neo(k), g -> sum')
   ylabel('Ki ')
end
hold on
plot(gene.x,datakp1.coef.ii,'m',gene.x,datakp1.neo.coef.ii,'c',gene.x,datakp1.neo.coef.ii+datakp1.coef.ii,'k')

subplot(3,2,6)
semilogy(eee,'or');
title(['dt = ', num2str(gene.dt)]);
xlabel(nb)
while(get(hc,'value') == 1)
   pause(1)
end
