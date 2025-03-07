    
% plot de control
    h=findobj(0,'type','figure','tag','jmoy');
if isempty(h)
	h=figure('color',[ 1 1 1],'defaultaxesfontsize',12,'tag','jmoy');
	hc=uicontrol(h,'style','radio','tag','stop');
else
	figure(h);
end
hc = findobj(h,'tag','stop');
% if nb >30
% 	set(hc,'value',1);
% end

subplot(2,1,1)
%plot(gene.x,jmoy-datak.prof.jmoy,'r',gene.x,equi.jmoy-datak.prof.jmoy,'m', ...
%      gene.x,datakp1.prof.jmoy-datak.prof.jmoy,'c');
plot(gene.x,jmoy,'r',gene.x,equi.jmoy,'m', ...
      gene.x,datakp1.prof.jmoy,'c');
xlabel('sqrt(Phi/pi/B0) normalise'); 
ylabel('Jmoy (A*m^-^2');
title(['new -> r, equi -> m, last -> c']);

subplot(2,1,2)
%plot(gene.x,equi.psi-datak.prof.psi,'m',gene.x,datakp1.prof.psi-datak.prof.psi,'c')
plot(gene.x,equi.psi,'m',gene.x,datakp1.prof.psi,'c')
xlabel('sqrt(Phi/pi/B0)');
ylabel('Psi');
title(['equi -> m, diff -> c']);

while(get(hc,'value') == 1)
   pause(1)
end
