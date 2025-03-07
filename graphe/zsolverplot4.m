    
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
ii = psi;

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
   plot(gene.x,datak.source.totale.j,'r',gene.x,datak.neo.jboot,'b')
	ylabel('source J et jboot')
   hold on
end
plot(gene.x,datakp1.source.totale.j,'m',gene.x,datakp1.neo.jboot,'k')
%netot =  datak.equi.rhomax .* trapz(gene.x,datakp1.prof.ne .* datak.equi.vpr,2);
%title(['netot = ',num2str(netot)])			           

terme1k_ = datak.equi.grho2r2./phys.mu0.*datak.coef.eta./datak.equi.rhomax.^2./ ...
          datak.equi.r2i.*datak.prof.psid2;
terme1kp1_ = datakp1.equi.grho2r2./phys.mu0.*datakp1.coef.eta./datakp1.equi.rhomax.^2./ ...
          datakp1.equi.r2i.*datakp1.prof.psid2;
lnck_ = log(datak.equi.vpr.*datak.equi.grho2r2./datak.equi.F);
dlnck_ = gradient(lnck_,gene.x);
lnckp1_ = log(datakp1.equi.vpr.*datakp1.equi.grho2r2./datakp1.equi.F);
dlnckp1_ = gradient(lnckp1_,gene.x);

terme2k_ = (datak.equi.grho2r2./phys.mu0.*datak.coef.eta./datak.equi.rhomax.^2./ ...
         datak.equi.r2i.*dlnck_ + gene.x./datak.equi.rhomax.* ...
			datak.equi.drhomaxdt).*datak.prof.psid1;
terme2kp1_ = (datakp1.equi.grho2r2./phys.mu0.*datakp1.coef.eta./datakp1.equi.rhomax.^2./ ...
         datakp1.equi.r2i.*dlnckp1_ + gene.x./datakp1.equi.rhomax.* ...
			datak.equi.drhomaxdt).*datak.prof.psid1;
terme3k_ = datak.geo.b0.*datak.coef.eta./datak.equi.F./datak.equi.r2i.* ...
           datak.source.totale.j;
terme3kp1_ = datakp1.geo.b0.*datakp1.coef.eta./datakp1.equi.F./datakp1.equi.r2i.* ...
           datakp1.source.totale.j;
dpsidtk_   =  terme1k_ +terme2k_ +terme3k_ ;
dpsidtkp1_   =  terme1kp1_ +terme2kp1_ +terme3kp1_ ;


subplot(3,2,3)
if premier_passage == 1
   hold off
   plot(gene.x,datak.prof.dpsidt,'r',gene.x,dpsidtk_,'b')
	ylabel('dpsidt')
   %hold on
end
plot(gene.x,datakp1.prof.dpsidt,'m',gene.x,dpsidtkp1_,'k')
	ylabel('dpsidt')


subplot(3,2,2)
if premier_passage == 1
   hold off
%   plot(gene.x,datak.coef.nn,'r',gene.x,datak.neo.coef.nn,'b',gene.x,datak.neo.coef.nn+datak.coef.nn,'g')
   plot(gene.x,datak.coef.eta,'r')
   ylabel('eta')
end
hold on
plot(gene.x,datakp1.coef.eta,'r')
 

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
