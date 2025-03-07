% script du  plot des confinement 0d
zs   = post.zerod;
cons = post.z0dinput.cons;
exp0d  = post.z0dinput.exp0d;
t    = zs.temps;
disrup = double(zs.disrup);
disrup(disrup == 0) = NaN;
disrup(isfinite(disrup)) = 0;
% rapport pour le facteur d'amplification (Ealpha +En)/Ealpha
rfan = (3.56e6 + 14.03e6) ./ 3.56e6 ;
padd = cons.picrh + cons.pecrh + real(cons.pnbi) + imag(cons.pnbi) + zs.pohm + zs.plh;


% Lin-Liu betaN limit
%ref : Y.R. Lin-Liu & R.D. Stambaugh NF 44 (2004) 548-554
A  = post.z0dinput.geo.R ./ post.z0dinput.geo.a;
betan_linliu = 10.0 .* (-0.7748  + 1.2869.* post.z0dinput.geo.K -0.2921 .*  post.z0dinput.geo.K .^ 2 + 0.0197 .* post.z0dinput.geo.K .^ 3) .* ...
               coth((1.8524 + 0.2319 .* post.z0dinput.geo.K) ./ A .^ 0.6163) ./ A .^ 0.5523;


% precalcul
if post.z0dinput.mode_exp == 0
	vpr  = data.equi.vpr;
	nvpr = trapz(param.gene.x,vpr,2);
	%esup_fus = data.equi.rhomax .* trapz(param.gene.x,data.source.fus.psupra .* vpr,2);
end
tauec = zs.tauee;
tauic = zs.tauii;

h = findobj(0,'type','figure','tag','z0plotc');
if isempty(h)
       h=figure('tag','z0plotc');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

k    = 9;
if post.z0dinput.mode_exp == 0

   subplot(k,1,1);	
   plot(t,zs.modeh,'or');
   ylabel('mode L/H')
   title(sprintf('Zerod : %s@%d/confinment ', ...
          post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);
   plot(t,zs.plhthr/1e6,'r',t,zs.pin/1e6,':c',t,zs.ploss/1e6,':m',t,zs.plossl2h ./ 1e6,'g', ...
        t,post.z0dinput.option.l2hmul + zs.plossl2h ./ 1e6,'k');
   ylabel('Pl2h  (MW)')
   legend('P_{Pl2h }','P_{IN}','P_{loss}','P_{scaling}', 'P_{scaling  + offset }','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,3);
   semilogy(t,zs.taue,'r',t,zs.tauthl,'m',t,zs.tauh,'c',data.gene.temps,data.gene.tauth,'b');
   ylabel('taue (s)')
   legend('tau_{E}','tau_{E,Lmode}','tau_{E,Hmode}','tau_{E,cronos}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   z0loglin(gca);
   subplot(k,1,4);
   semilogy(t,zs.tauhe,'r',t,zs.tauhe_l,'m',t,zs.tauhe_h,'c');
   z0loglin(gca);
   ylabel('tau_H_e (s)')
   legend('tau_{He}','tau_{He,Lmode}','tau_{He,Hmode}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,5);
   plot(t,rfan .* zs.pfus ./ padd,'r', ...
	data.gene.temps,rfan .* data.gene.paddfus ./ (data.gene.paddtot - data.gene.paddfus),'b');
   vmax = max(rfan .* zs.pfus ./ padd);
   if isfinite(vmax) & (vmax >0)
      set(gca,'ylim',[0,vmax]);
   end
   ylabel('Q factor ')
   legend('P_{fusion}','P_{fusion,cronos}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,6);	
   semilogy(t,tauec,'m',t,tauic,'c');
   ylabel('tau scale (c) (s)')
   z0loglin(gca);
   legend('tau_{electron}','tau_{ion}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,7);
   plot(t,zs.wrad./ 2 ./ pi ./ 1e3,'r');
   ylabel('<\Omega _{\phi}> (2\pi 10^3 rad/s)')
   legend('Metis','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   %z0loglin(gca);
   subplot(k,1,8)
   plot(t,zs.hitb,'b',t,zs.hmhd,'-.r',t,disrup,'om');
   set(gca,'ylim',[0,2]);
   ylabel('H_H factor');
   legend('H_{H,itb}','H_{H,mhd}','disruption indicator','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,9)
   plot(t,zs.betaptot + zs.li /2,'b',t,zs.betan .* 100,'r',t,4 .* zs.li,'m',t,6 .* zs.li,'g',t,betan_linliu,'c');
   ylabel('Beta');
   legend('Beta_p + l_i/2','Beta_{N}','4 * l_i','6 * l_i','Lin-Liu limit(scaling)','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   xlabel('time (s)');

else
   subplot(k,1,1);	
   plot(t,zs.modeh,'or',t,exp0d.modeh,':b');
   legend('H mode flag','H_{\alpha , normalized}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('mode L/H')
   title(sprintf('Zerod : %s@%d/confinment ', ...
          post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,zs.plhthr/1e6,'r',t,zs.pin/1e6,':c',t,zs.ploss/1e6,':m',t,zs.plossl2h ./ 1e6,'g', ...
        t,post.z0dinput.option.l2hmul + zs.plossl2h ./ 1e6,'k');
   ylabel('Pl2h (MW)')
   legend('P_{Pl2h }','P_{IN}','P_{loss}','P_{scaling}', 'P_{scaling  + offset }','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,3);	
   semilogy(t,zs.taue,'r',t,zs.tauthl,'m',t,zs.tauh,'c',t,exp0d.taue,'b');
   ylabel('taue (s)')
   legend('tau_{E}','tau_{E,Lmode}','tau_{E,Hmode}','tau_{E,exp}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   z0loglin(gca);
   subplot(k,1,4);
   semilogy(t,zs.tauhe,'r',t,zs.tauhe_l,'m',t,zs.tauhe_h,'c');
   ylabel('tau_H_e (s)')
   legend('tau_{He}','tau_{He,Lmode}','tau_{He,Hmode}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   z0loglin(gca);
   subplot(k,1,5);
   plot(t,rfan .* zs.pfus ./ padd,'r',t,rfan .* exp0d.pfus ./ (exp0d.pin - exp0d.pfus),'b');
   set(gca,'ylim',[0,max(1,max(rfan .* zs.pfus ./ padd))]);
   ylabel('Q factor ')
   legend('P_{fusion}','P_{fusion,exp}','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,6);	
   semilogy(t,tauec,'m',t,tauic,'c');
   ylabel('tau  scale (s)')
   z0loglin(gca);
   legend('tau_{electron}','tau_{ion}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,7);
   switch post.z0dinput.machine
   case 'JET'
   	% c'est la valeur centrale pour l'impurete principale
   	plot(post.profil0d.temps,post.profil0d.vtor(:,1)./ 2 ./pi ./1e3./ post.profil0d.Raxe(:,1),'r', ...
	     t,exp0d.wrad./ 2 ./ pi ./ 1e3,'b',t,zs.wrad./ 2 ./pi ./1e3,'g');	
  	legend('Metis centre','Experimental centre','Metis volume averaged','Location','north','orientation','horizontal')
   	legend(gca,'boxoff')
   otherwise
   	plot(t,zs.wrad./ 2 ./pi ./1e3,'r',t,exp0d.wrad./ 2 ./ pi ./ 1e3,'b:');
   	legend('Metis','Experimental','Location','north','orientation','horizontal')
   	legend(gca,'boxoff')
   end
   ylabel('<\Omega _{\phi}> (2\pi 10^3 rad/s)')
   subplot(k,1,8)
   plot(t,zs.hitb,'b',t,zs.hmhd,'-.r',t,disrup,'om');
   set(gca,'ylim',[0,2]);
   ylabel('H_H factor');
   legend('H_{H,itb}','H_{H,mhd}','disruption indicator','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   subplot(k,1,9)
   plot(t,zs.betaptot + zs.li /2,'b',t,zs.betan .* 100,'r',t,4 .* zs.li,'m',t,6 .* zs.li,'g',t,betan_linliu,'c');
   ylabel('Beta');
   legend('Beta_p + l_i/2','Beta_{N}','4 * l_i','6 * l_i','Lin-Liu limit(scaling)','Location','north','orientation','horizontal')
   legend(gca,'boxoff')
   xlabel('time (s)');

end
joint_axes(h,k);
