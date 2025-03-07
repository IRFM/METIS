% script de plot des source de puissance
x=param.gene.x;
t=data.gene.temps;

figure('color',[1 1 1 ],'defaultaxesyscale','log','defaultaxesfontsize',12);

subplot(3,3,1)
plot(t,sum(abs(data.cons.fci'))','r',t,data.gene.paddfci,'oc');
title('FCI')
xlabel('sqrt(Phi/pi/B0) normalise');
ylabel('Padd (W)');
legend('consigne','calcul');

subplot(3,3,2)
plot(t,sum(abs(data.cons.hyb'))','r',t,data.gene.paddhyb,'oc');
title('Hybride')
xlabel('sqrt(Phi/pi/B0) normalise');
ylabel('Padd (W)');
legend('consigne','calcul');

subplot(3,3,3)
plot(t,sum(abs(data.cons.fce'))','r',t,data.gene.paddfce,'oc');
title('FCE')
xlabel('sqrt(Phi/pi/B0) normalise');
ylabel('Padd (W)');
legend('consigne','calcul');

subplot(3,3,4)
plot(t,sum(abs(data.cons.idn'))','r',t,data.gene.paddidn,'oc');
title('IDN')
xlabel('sqrt(Phi/pi/B0) normalise');
ylabel('Padd (W)');
legend('consigne','calcul');

subplot(3,3,5)
r = data.gene.prad + data.gene.pbrem + data.gene.pcyclo;
n = abs(-data.gene.qei + data.gene.qneo);
semilogy(t,data.gene.paddohm,'r',t,data.gene.paddfus,'c',t,r,'m',t,n,'k',t,data.gene.pn0,'g');
title('autres sources')
xlabel('sqrt(Phi/pi/B0) normalise');
ylabel('Padd (W)');
legend('ohm','fus','rayon','neo','n0');
set(gca,'ylim',[0 inf]);

subplot(3,3,6)
r = data.gene.padde + data.gene.paddion;
n = data.gene.paddtot - data.gene.prad - data.gene.pbrem  ...
    - data.gene.pcyclo + data.gene.pn0;
plot(t,r,'r',t,n,'oc');
title('balance totale')
xlabel('sqrt(Phi/pi/B0) normalise');
ylabel('Padd (W)');
legend('totale','reconstituee');

subplot(3,3,7)
r = data.gene.padde - data.gene.prad - data.gene.pbrem  ...
    - data.gene.pcyclo - data.gene.qei + data.gene.qneo;
plot(t,data.gene.pel,'r',t,r,'oc');
title('balance electronique')
xlabel('sqrt(Phi/pi/B0) normalise');
ylabel('Padd (W)');
legend('totale','reconstituee');

subplot(3,3,8)
r = data.gene.paddion + data.gene.qei - data.gene.qneo;
plot(t,data.gene.pion,'r',t,r,'oc');
title('balance electronique')
xlabel('sqrt(Phi/pi/B0) normalise');
ylabel('Padd (W)');
legend('totale','reconstituee');

subplot(3,3,9)
el = data.prof.flux.qe(:,end) .* data.equi.vpr(:,end);
ion = data.prof.flux.qi(:,end) .* data.equi.vpr(:,end);
elp = data.gene.padde - data.gene.prad - data.gene.pbrem  ...
      - data.gene.pcyclo - data.gene.qei + data.gene.qneo;
ionp = data.gene.paddion + data.gene.qei - data.gene.qneo;
plot(t,el,'r',t,elp,'or',t,data.gene.pel,'r+',t,ion,'c',t,ionp,'oc',t,data.gene.pion,'c+', ...
     t,el+ion,'g',t,data.gene.pel+data.gene.pion,'og');
title('balance flux convecte')
xlabel('sqrt(Phi/pi/B0) normalise');
ylabel('Padd (W)');
legend('flux el','sum el','source el','flux ion','sum el','source ion','flux tot','source tot');

