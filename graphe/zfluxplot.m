    
% plot de control
h=findobj(0,'type','figure','tag','flux');
if isempty(h)
	h=figure('color',[ 1 1 1],'defaultaxesfontsize',12,'tag','flux');
	hc=uicontrol(h,'style','radio','tag','stop');
else
	figure(h);
end
hc = findobj(h,'tag','stop');
% if nb >30
% 	set(hc,'value',1);
% end

% calcul des coefficients a partir des flux
keeff      =  - datak.prof.flux.qe ./ datak.prof.gte ./ phys.e .* datak.equi.grho ; 
kieff      =  - datak.prof.flux.qi  ./ datak.prof.gti ./ phys.e .* datak.equi.grho ; 
chieeff    =  keeff ./ datak.prof.ne;
chiieff    =  kieff ./ datak.prof.ni;
deeff      =  - datak.prof.flux.ge ./  datak.prof.gne ./ datak.equi.grho; 
chiean     =  datak.coef.ee ./ datak.prof.ne;
chiian     =  datak.coef.ii ./ datak.prof.ni;
dean       =  datak.coef.nn;
chieneo     =  datak.neo.coef.ee ./ datak.prof.ne;
chiineo     =  datak.neo.coef.ii ./ datak.prof.ni;
deneo       =  datak.neo.coef.nn;


subplot(2,2,1)
hold off
plot(x,chieeff,'linestyle','-','color',[1 0 0]);
hold on
plot(x,chiean,'linestyle','-','color',[0 1 1]);
plot(x,chieneo,'linestyle','-','color',[0 0 0]);
title('Chie')
xlabel('sqrt(Phi/pi/B0)');
ylabel('Chie (m^2*s^-^1');
title(['eff -> r, an -> c, neo -> k']);

subplot(2,2,2)
hold off
plot(x,chiieff,'linestyle','-','color',[1 0 0]);
hold on
plot(x,chiian,'linestyle','-','color',[1 0 1]);
plot(x,chiineo,'linestyle','-','color',[0 0 0]);
title('Chii')
xlabel('sqrt(Phi/pi/B0)');
ylabel('Chii (m^2*s^-^1');
title(['eff -> r, an -> c, neo -> k']);

subplot(2,2,3)
hold off
plot(x,deeff,'linestyle','-','color',[1 0 0]);
hold on
plot(x,dean,'linestyle','-','color',[1 0 1]);
plot(x,deneo,'linestyle','-','color',[0 0 0]);
title('De')
xlabel('sqrt(Phi/pi/B0)');
ylabel('De (m^2*s^-^1');
title(['eff -> r, an -> c, neo -> k']);

subplot(2,2,4)
hold off
plot(x,abs(datak.coef.ve),'linestyle','-','color',[1 0 0]);
hold on
plot(x,abs(datak.coef.vi),'linestyle','-','color',[1 0 1]);
plot(x,abs(datak.coef.vn),'linestyle','-','color',[0 1 1]);
plot(x,abs(datak.neo.coef.ve),'linestyle',':','color',[1 0 0]);
plot(x,abs(datak.neo.coef.vi),'linestyle',':','color',[1 0 1]);
plot(x,abs(datak.neo.coef.vn),'linestyle',':','color',[0 1 1]);
title('V convection')
xlabel('sqrt(Phi/pi/B0)');
ylabel('V (m*s^-^1');

while(get(hc,'value') == 1)
   pause(1)
end
