    
% plot de control
h=findobj(0,'type','figure','tag','fp');
if isempty(h)
	h=figure('color',[ 1 1 1],'defaultaxesfontsize',12,'tag','fp');
	hc=uicontrol(h,'style','pop','tag','stop','string','go|halt|keyboard');
else
	figure(h);
	hc = findobj(h,'type','uicontrol');
end
%ii = ne;
ii = pe;

subplot(3,2,1)
if premier_passage == 1
   hold off
   plot(gene.x,(fpini(:,ii)-F(:,ii)) ./ std(F(:,ii)),'r');
   xlabel(['dt = ',num2str(dt)]); 
	title('Delta psi'); 
   hold on
   axis([0 1 -5e-3 5e-3]);
end
plot(gene.x,(fp(:,ii)-fpini(:,ii)) ./ std(F(:,ii)));

subplot(3,2,4)
if premier_passage == 1
   hold off
   plot(gene.x,datak.source.totale.el,'r')
	ylabel('source J et jboot')
   hold on
end
plot(gene.x,datakp1.source.totale.el,'m')
%netot =  datak.equi.rhomax .* trapz(gene.x,datakp1.prof.ne .* datak.equi.vpr,2);
%title(['netot = ',num2str(netot)])			           



subplot(3,2,3)
if premier_passage == 1
   hold off
   plot(gene.x,datak.prof.pe,'r',gene.x,datak.prof.pion,'b')
	ylabel('Pa')
   %hold on
end
plot(gene.x,datakp1.prof.pe,'m',gene.x,datakp1.prof.pion,'k')
	ylabel('Pa')


subplot(3,2,2)
if premier_passage == 1
   hold off
   plot(gene.x,datak.coef.ee./datak.prof.ne,'r')
   ylabel('coeff ee')
end
hold on
   plot(gene.x,datakp1.coef.ee./datakp1.prof.ne,'r')

subplot(3,2,5)
if premier_passage == 1
   hold off
   plot(gene.x,datak.prof.jeff,'r',gene.x,datak.prof.epar./datak.coef.eta+datak.source.totale.j,'b')
  ylabel('Jeff et E/eta+Jni')
end
%hold on
plot(gene.x,datakp1.prof.jeff,'r',gene.x,datakp1.prof.epar./datakp1.coef.eta+datakp1.source.totale.j,'b')
ylabel('Jeff et E/eta+Jni')

subplot(3,2,6)
semilogy(eee,'or');
title(['e = ', num2str(erreur)]);
xlabel(nb)

while(get(hc,'value') > 1)
   pause(1)
	if get(hc,'value') == 3
	   fprintf('Keyboard');
	   keyboard
		set(hc,'value',2);
  end
end
