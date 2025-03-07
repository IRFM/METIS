% script du  plot de LH 0d
zs   = post.zerod;
cons = post.z0dinput.cons;
exp0d  = post.z0dinput.exp0d;
if isfield(exp0d,'edgeflux')
  exp0d.flux = exp0d.edgeflux;
end
geo = post.z0dinput.geo;
t    = zs.temps;
% precalcul
vloop = zs.vloop;
vloop(find(~isfinite(vloop))) = sqrt(-1);
flux = -cumtrapz(t,vloop,1) ./2 ./ pi;
flux(find(imag(vloop))) = NaN;
flux = real(flux);

h = findobj(0,'type','figure','tag','z0plotlh');
if isempty(h)
       h=figure('tag','z0plotlh');
else
       figure(h);
end
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

k    = 9;
if post.z0dinput.mode_exp == 0
   fluxexp = data.prof.psi(:,end);
   flux0 = mean(fluxexp(isfinite(fluxexp))) - mean(flux(isfinite(flux)));

   subplot(k,1,1);	
   if post.z0dinput.option.lhmode == 5
	plot(t,0 .* zs.plh_th./1e6,'r',data.gene.temps,data.gene.paddhyb./1e6,':b');
   else
	plot(t,zs.plh_th./1e6,'r',data.gene.temps,data.gene.paddhyb./1e6,':b');
   end
   ylabel('Plh (MW)')
   title(sprintf('Zerod : %s@%d/LH ( r -> 0D, m -> eta0 , c -> eta1, b: -> cronos)', ...
             post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,flux,'r',data.gene.temps,data.prof.psi(:,end) - flux0 ,':b');
   ylabel('Flux (Wb/rad)');
   subplot(k,1,3);   
   if post.z0dinput.option.lhmode == 5
      plot(t,0.* zs.ilh ./ zs.ip,'r',data.gene.temps,data.gene.ihyb ./ data.gene.ip,':b');
   else
      plot(t,zs.ilh ./ zs.ip,'r',data.gene.temps,data.gene.ihyb ./ data.gene.ip,':b');
   end
   ylabel('Ilh/Ip ');
   subplot(k,1,4);	
   plot(t,zs.nem./1e19,'r',data.gene.temps,data.gene.nemoy./1e19,':b');
   ylabel('<Ne> (1e19 m^-^3)');
   subplot(k,1,5);	
   plot(t,geo.R,'r',data.gene.temps,data.geo.r0,':b');
   ylabel('R0 (m)');
   subplot(k,1,6);	
   plot(t,zs.tem./1e3,'r',data.gene.temps,data.gene.temoy./1e3,':b');
   ylabel('<Te> (keV)');
   subplot(k,1,7);	
   plot(t,zs.zeff./zs.zmszl,'r',data.gene.temps,data.gene.zeffm,':b');
   ylabel('<Zeff>');
   subplot(k,1,8);	
   plot(t,zs.betaptot + zs.li./2,'r',data.gene.temps,data.gene.betap + data.gene.li./2,':b');
   ylabel('Beta_p + li/2');
   subplot(k,1,9);	
   plot(t,zs.etalh0./1e19,'m',t,zs.etalh1./1e19,'c',t,zs.etalh./1e19,'r');
   ylabel('eta_L_H (1e19 A W^ -1 m^-^2)');
   xlabel('time (s)');
else
   if isfield(exp0d,'flux')
      flux0 = mean(flux(isfinite(flux))) - mean(exp0d.flux(isfinite(exp0d.flux)));
   else
      flux0 = 0;
      exp0d.flux = NaN .* t;
   end

   subplot(k,1,1);	
   if post.z0dinput.option.lhmode == 5
      plot(t,0 .* zs.plh_th./1e6,'r',t,cons.plh./1e6,':b');
   else
      plot(t,zs.plh_th./1e6,'r',t,cons.plh./1e6,':b');
   end
   ylabel('Plh (MW)')
   title(sprintf('Zerod : %s@%d/LH ( r -> 0D, m -> eta0 , c -> eta1, b: -> experiment)', ...
             post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,flux,'r',t,exp0d.flux + flux0,':b');
   ylabel('Flux (Wb/rad)');
   subplot(k,1,3);	
   if post.z0dinput.option.lhmode == 5
      plot(t,0 .* zs.ilh ./ zs.ip,'r',t,exp0d.ilh ./ exp0d.ip,':b');
   else
      plot(t,zs.ilh ./ zs.ip,'r',t,exp0d.ilh ./ exp0d.ip,':b');
   end
   ylabel('Ilh/Ip ');
   subplot(k,1,4);	
   plot(t,zs.nem./1e19,'r',t,exp0d.nem./1e19,':b');
   ylabel('<Ne> (1e19 m^-^3)');
   subplot(k,1,5);	
   plot(t,geo.R,'r');
   ylabel('R0 (m)');
   subplot(k,1,6);	
   plot(t,zs.tem./1e3,'r',t,exp0d.tem./1e3,':b');
   ylabel('<te> (keV)');
   subplot(k,1,7);	
   plot(t,zs.zeff./zs.zmszl,'r',t,exp0d.zeff,':b');
   ylabel('<Zeff>');
   subplot(k,1,8);	
   plot(t,zs.betaptot + zs.li./2,'r',t,exp0d.betap + exp0d.li./2,':b');
   ylabel('Beta_p + li/2');
   subplot(k,1,9);	
   plot(t,zs.etalh0./1e19,'m',t,zs.etalh1./1e19,'c',t,zs.etalh./1e19,'r');
   ylabel('eta_L_H (1e19 A W^ -1 m^-^2)');
   xlabel('time (s)');
   set(gca,'ylim',[0,3]);

end
joint_axes(h,k);
