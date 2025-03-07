% script du  plot des puissance 0d
h = findobj(0,'type','figure','tag','z0plotp');
if isempty(h)
       h=figure('tag','z0plotp');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

zs   = post.zerod;
cons = post.z0dinput.cons;
exp0d  = post.z0dinput.exp0d;
t    = zs.temps;
disrup = double(zs.disrup);
disrup(disrup == 0) = NaN;
disrup(isfinite(disrup)) = 0;
k    = 9;
if post.z0dinput.mode_exp == 0

   subplot(k,1,1);	
   plot(t,zs.ploss/1e6,'r',t,zs.pth./1e6,'m',data.gene.temps,data.gene.ploss ./ 1e6,'b');
   ylabel('Ploss (MW)')
    legend('P_{loss}','P_{th}','P_{loss,cronos}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
  title(sprintf('Zerod : %s@%d/power ', ...
             post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,zs.pfus./1e6,'r',t,zs.pfus_th./1e6,'m',t,zs.pfus_nbi./1e6,'g',t,(zs.pfus + zs.pfus_loss) ./1e6,'c',data.gene.temps,data.gene.paddfus./1e6,'b');
   ylabel('Palpha (MW)');
   legend('P_{alpha}','P_{alpha,th}','P_{DD & DT,nbi}','P_{alpha + alpha losses}','P_{alpha,cronos}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,3);	
   plot(t,zs.prad./1e6,'g',t,zs.pradsol./1e6,'m',t,zs.pioniz./1e6,'c',t,(zs.prad+zs.pradsol + zs.pioniz)./1e6,'r', ...
          data.gene.temps,data.gene.prad./1e6,'b',t,disrup,'om');
   ylabel('[MW]');
   legend('P_{line}','P_{SOL}','P_{ioniz}','P_{rad+sol+ioniz}','P_{rad,cronos}','disruption','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('Prad (g), Prad SOL (m)+ disrup (om) [MW]');
   subplot(k,1,4);	
   plot(t,zs.pbrem./1e6,'r',data.gene.temps,data.gene.pbrem./1e6,'b');
   ylabel('Pbrem (MW)');
    legend('P_{brem}','P_{brem,cronos}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
  subplot(k,1,5);	
   plot(t,zs.pcyclo./1e6,'r',data.gene.temps,data.gene.pcyclo./1e6,'b');
   ylabel('Pcyclo (MW)');
   legend('P_{cyclo}','P_{cyclo,cronos}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,6);	
   plot(t,zs.picrh./1e6,'r',t,zs.picrh_th./1e6,'m',data.gene.temps,data.gene.paddfci./1e6,'b');
   ylabel('PICRH (MW)');
    legend('P_{ICRH}','P_{ICRH,th}','P_{ICRH,cronos}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
  subplot(k,1,7);	
   plot(t,real(zs.pnbi)./1e6 + imag(zs.pnbi)./1e6,'r',t,real(zs.pnbi_th)./1e6 + imag(zs.pnbi_th)./1e6,'m', ...
        data.gene.temps,data.gene.paddidn./1e6,'b');
   ylabel('PNBI (MW)');
    legend('P_{NBI}','P_{NBI,th}','P_{NBI,cronos}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
  subplot(k,1,8);	
   plot(t,zs.pecrh./1e6,'r',data.gene.temps,data.gene.paddfce./1e6,'b');
   ylabel('PECRH (MW)');
   legend('P_{ECRH}','P_{ECRH,cronos}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,9);	
   plot(t,zs.plh./1e6,'r',t,zs.plh_th./1e6,'r',data.gene.temps,data.gene.paddhyb./1e6,'b');
   ylabel('PLH (MW)');
   if post.z0dinput.option.lhmode == 5
	legend('P_{ECRH 2}','P_{ECRH 2,th}','P_{LH,cronos}','Location','north','orientation','horizontal');
   else
	legend('P_{LH}','P_{LH,th}','P_{LH,cronos}','Location','north','orientation','horizontal');
   end
   legend(gca,'boxoff')
   xlabel('time (s)');

else
   subplot(k,1,1);	
   plot(t,zs.ploss/1e6,'r',t,zs.pth./1e6,'m',t,exp0d.ploss ./ 1e6,'b');
   legend('P_{loss}','P_{th}','P_{loss,exp}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('Ploss (MW)')
   title(sprintf('Zerod : %s@%d/power', ...
             post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,zs.pfus./1e6,'r',t,zs.pfus_th./1e6,'m',t,zs.pfus_nbi./1e6,'g',t,(zs.pfus + zs.pfus_loss)./1e6,'c',t,exp0d.pfus./1e6,'b');
   ylabel('Palpha (MW)');
   legend('P_{alpha}','P_{alpha,th}','P_{DD & DT,nbi}','P_{alpha + alpha losses}','P_{alpha,exp}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,3);	
   plot(t,zs.prad./1e6,'g',t,zs.pradsol./1e6,'m',t,zs.pioniz./1e6,'c',t,(zs.prad+zs.pradsol + zs.pioniz)./1e6,'r',t,exp0d.prad./1e6,'b',t,disrup,'om');
   ylabel('Prad (MW)');
   legend('P_{line}','P_{SOL}','P_{ioniz}','P_{rad+sol+ioniz}','P_{rad,exp}','disruption','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,4);	
   plot(t,zs.pbrem./1e6,'r',t,exp0d.pbrem./1e6,'b');
   legend('P_{brem}','P_{brem,exp}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('Pbrem (MW)');
   subplot(k,1,5);	
   plot(t,zs.pcyclo./1e6,'r',t,exp0d.pcyclo./1e6,'b');
   ylabel('Pcyclo (MW)');
   legend('P_{cyclo}','P_{cyclo,exp}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,6);	
   plot(t,zs.picrh./1e6,'r',t,zs.picrh_th./1e6,'m',t,exp0d.picrh./1e6,'b');
   ylabel('PICRH (MW)');
   legend('P_{ICRH}','P_{ICRH,th}','P_{ICRH,exp}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,7);	
   plot(t,real(zs.pnbi)./1e6 + imag(zs.pnbi)./1e6,'r',t,real(zs.pnbi_th)./1e6 + imag(zs.pnbi_th)./1e6,'m', ...
        t,exp0d.pnbi./1e6,'b');
   ylabel('PNBI (MW)');
   legend('P_{NBI}','P_{NBI,th}','P_{NBI,exp}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,8);	
   plot(t,zs.pecrh./1e6,'r',t,exp0d.pecrh./1e6,'b');
   if post.z0dinput.option.lhmode == 5
      ylabel('PECRH 1 (MW)');
      legend('P_{ECRH 1}','P_{ECRH 1,exp}','Location','north','orientation','horizontal')
   else
      ylabel('PECRH (MW)');
      legend('P_{ECRH}','P_{ECRH,exp}','Location','north','orientation','horizontal')
   end
   legend(gca,'boxoff')
   subplot(k,1,9);	
   if post.z0dinput.option.lhmode == 5
      plot(t,zs.plh./1e6,'r',t,exp0d.plh./1e6,'b');
      ylabel('PECRH 2 (MW)');
      legend('P_{ECRH}','P_{ECRH,exp}','Location','north','orientation','horizontal')
   else
      plot(t,zs.plh./1e6,'r',t,zs.plh_th./1e6,'m',t,exp0d.plh./1e6,'b');
      ylabel('PLH (MW)');
      legend('P_{LH}','P_{LH,th}','P_{LH,exp}','Location','north','orientation','horizontal')
   end
   legend(gca,'boxoff')
   xlabel('time (s)');

end

hh = findobj(h,'type','axes');
for l=1:k
   lim = get(hh(l),'ylim');
   set(hh(l),'ylim',[0,max(max(lim),0.1)]);
end
joint_axes(h,k);
