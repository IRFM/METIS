% script du  plot des temperature 0d
zs   = post.zerod;
cons = post.z0dinput.cons;
exp0d  = post.z0dinput.exp0d;
t    = zs.temps;
% precalcul
if post.z0dinput.mode_exp == 0
	vpr  = data.equi.vpr;
	nvpr = trapz(param.gene.x,vpr,2);
	%   n1m = trapz(param.gene.x,vpr .* squeeze(sum(data.impur.impur(:,:,ind1),3)),2)./ nvpr;
end 

piqn = post.profil0d.nep(:,1) ./ (trapz(post.profil0d.xli,post.profil0d.vpr .* post.profil0d.nep,2) ./ trapz(post.profil0d.xli,post.profil0d.vpr,2));

h = findobj(0,'type','figure','tag','z0plott');
if isempty(h)
       h=figure('tag','z0plott');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

k    = 8;
if post.z0dinput.mode_exp == 0
  
   tite  = trapz(param.gene.x,data.equi.vpr .* data.prof.ti,2) ./trapz(param.gene.x,data.equi.vpr .* data.prof.te,2);

   subplot(k,1,1);	
   plot(t,zs.tem/1e3,'r',t,zs.teped/1e3,'m',data.gene.temps,data.gene.temoy ./ 1e3,':b');
   ylabel('<Te> (r) & Te_{ped} (m) (keV)')
   title(sprintf('Zerod : %s@%d/temperature (r -> metis , b: -> cronos, g -> goal/reference)', ...
             post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,zs.te0/1e3,'r',t,zs.tem .* (1 + zs.ate)/1e3,'m:',data.gene.temps,data.prof.te(:,1) ./ 1e3,':b');
   ylabel('Te0 (keV)')
   subplot(k,1,3);	
   plot(t,zs.tebord,'r',data.gene.temps,data.prof.te(:,end) ,':b');
   ylabel('Tea (eV)')
   subplot(k,1,4);	
   plot(t,zs.pel./1e6,'r',t,(zs.pel-zs.pei)./1e6,':r',data.gene.temps,data.gene.pel./1e6,':b');
   ylabel('P_{el} (w/wo Pei, MW)');
   subplot(k,1,5);	
   plot(t,zs.pion./1e6,'r',t,(zs.pion+zs.pei)./1e6,':r',data.gene.temps,data.gene.pion./1e6,':b');
   ylabel('P_{ion} (w/wo Pei, MW)');
   subplot(k,1,6);	
   %plot(t,zs.tite,'r',t,zs.tite .* (1+zs.ane) ./ (1+zs.ate),'r:',data.gene.temps,data.prof.ti(:,1) ./ data.prof.te(:,1),':b');
   plot(t,zs.tite,'r',data.gene.temps,tite,':b');
   ylabel('Ti/Te');
   subplot(k,1,7);	
   plot(t,zs.ane + 1,'g',data.gene.temps,data.gene.piqne,':b',post.profil0d.temps,piqn,'r');
   ylabel('Ne0 /<Ne>');
   subplot(k,1,8);	
   plot(t,zs.ate + 1,'r',data.gene.temps,data.gene.piqte,':b');
   ylabel('Te0 /<Te>');
   xlabel('time (s)');
  
else
   subplot(k,1,1);	
   plot(t,zs.tem/1e3,'r',t,zs.teped/1e3,'m',t,exp0d.tem ./ 1e3,':b');
   ylabel('<Te> (r) & Te_{ped} (m) (keV)')
   title(sprintf('Zerod : %s@%d/temperature (r -> metis, b: -> experiment, g -> goal/reference)', ...
             post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,zs.te0/1e3,'r',t,zs.tem .* (1 + zs.ate)/1e3,'m:',t,exp0d.te0 ./ 1e3,':b');
   ylabel('Te0 (keV)')
   subplot(k,1,3);	
   plot(t,zs.tebord,'r',t,exp0d.tebord ,':b');
   ylabel('Tea (eV)')
   subplot(k,1,4);	
   plot(t,zs.pel./1e6,'r',t,(zs.pel-zs.pei)./1e6,':r',t,exp0d.pel./1e6,':b');
   ylabel('P_{el} (w/wo Pei, MW)');
   subplot(k,1,5);	
   plot(t,zs.pion./1e6,'r',t,(zs.pion+zs.pei)./1e6,':r',t,exp0d.pion./1e6,':b');
   ylabel('P_{ion} (w/wo Pei, MW)');
   subplot(k,1,6);	
   plot(t,zs.tite,'r',t,exp0d.tite,':b');
   %plot(t,zs.tite,'r',t,zs.tite .* (1+zs.ane) ./ (1+zs.ate),'m',t,exp0d.tite,':b');
   ylabel('<Ti/Te>');
   subplot(k,1,7);	
   plot(t,zs.ane + 1,'g',t,exp0d.ane + 1,':b',post.profil0d.temps,piqn,'r');
   ylabel('Ne0 /<Ne>');
   subplot(k,1,8);	
   plot(t,zs.ate + 1,'r',t,exp0d.ate + 1,':b');
   ylabel('Te0 /<Te>');
   xlabel('time (s)');
  

end
hh = findobj(h,'type','axes');
for l=1:k
   lim = get(hh(l),'ylim');
   set(hh(l),'ylim',[0,max(max(lim),0.1)]);
end
joint_axes(h,k);
