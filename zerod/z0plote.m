% script du  plot des energy 0d
zs   = post.zerod;
cons = post.z0dinput.cons;
exp0d  = post.z0dinput.exp0d;
t    = zs.temps;
% precalcul
if post.z0dinput.mode_exp == 0
	vpr  = data.equi.vpr;
	nvpr = trapz(param.gene.x,vpr,2);
	esup_fus = (3/2) .* data.equi.rhomax .* trapz(param.gene.x,data.source.fus.psupra .* vpr,2) +  ...
	           	data.equi.rhomax .* trapz(param.gene.x,data.source.fus.paniso .* vpr,2);
	esup_nbi = (3/2) .* data.equi.rhomax .* trapz(param.gene.x,data.source.idn.psupra .* vpr,2) + ...
			data.equi.rhomax .* trapz(param.gene.x,data.source.idn.paniso .* vpr,2);
	esup_icrh = (3/2) .* data.equi.rhomax .*trapz(param.gene.x,data.source.fci.psupra .* vpr,2)+  ...
	           	data.equi.rhomax .* trapz(param.gene.x,data.source.fci.paniso .* vpr,2);
end


% calcul de pped
wped   = 3./ 2 .* zs.pped .* trapz(linspace(0,1,21),cat(2,ones(1,20),0),2) .* zs.vp;
wpedmax = 3./ 2 .* zs.ppedmax .* trapz(linspace(0,1,21),cat(2,ones(1,20),0),2) .* zs.vp;

h = findobj(0,'type','figure','tag','z0plote');
if isempty(h)
       h=figure('tag','z0plote');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

k    = 9;
if post.z0dinput.mode_exp == 0
   subplot(k,1,1);	
   plot(t,zs.wdia/1e6,'r',t,zs.w./1e6,'m',data.gene.temps,data.gene.wdia ./ 1e6,'b');
   legend('W_{dia}','W_{tot}','W_{dia, CRONOS}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MJ')
   title(sprintf('Zerod : %s@%d/energy', ...
          post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,zs.wth/1e6,'r',data.gene.temps,data.gene.wth ./ 1e6,'b',t,wped./1e6,'m',t,wpedmax/1e6,'c');
   legend('W_{th}','W_{th, CRONOS}','W_{ped}','W_{ped, max}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MJ')
   subplot(k,1,3);	
   plot(t,zs.wbp/1e6,'r',data.gene.temps,data.gene.wbp ./ 1e6,'b');
   legend('W_{Bpol}','W_{Bpol, CRONOS}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MJ')
   subplot(k,1,4);	
   plot(t,zs.dwdt/1e6,'r',data.gene.temps,data.gene.dwdiadt ./ 1e6,':b');
   legend('dWdt_{tot}','dWdt_{tot, CRONOS}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MW')
   subplot(k,1,5);	
   plot(t,zs.dwthdt/1e6,'r',data.gene.temps,data.gene.dwthdt ./ 1e6,'b');
   legend('dWdt_{th}','dWdt_{th, CRONOS}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MW')
   subplot(k,1,6);	
   plot(t,zs.dwbpdt/1e6,'r',data.gene.temps,data.gene.dwbpdt ./ 1e6,'b');
   legend('dWdt_{Bpol}','dWdt_{Bpol, CRONOS}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MW')
   subplot(k,1,7);	
   plot(t,zs.esup_fus/1e6,'r',data.gene.temps,esup_fus ./ 1e6,'b');
   legend('Esupra_{ \alpha}','Esupra_{ \alpha, CRONOS}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MJ')
   subplot(k,1,8);	
   plot(t,zs.esup_icrh/1e6,'r',data.gene.temps,esup_icrh ./ 1e6,'b');
   legend('Esupra_{ICRH}','Esupra_{ICRH, CRONOS}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MJ')
   subplot(k,1,9);	
   plot(t,real(zs.esup_nbi)/1e6 + imag(zs.esup_nbi)/1e6,'r',data.gene.temps,esup_nbi ./ 1e6,'b');
   legend('Esupra_{NBI}','Esupra_{NBI, CRONOS}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MJ')
   xlabel('time (s)');
else
   subplot(k,1,1);
   plot(t,zs.wdia/1e6,'r',t,zs.w/1e6,'b-.',t,exp0d.wdia ./ 1e6,'m',t,exp0d.w ./ 1e6,'c-.');
   legend('W_{dia}','W_{tot}','W_{dia, exp}','W_{MHD, exp}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MJ')
   title(sprintf('Zerod : %s@%d/energy', ...
          post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,zs.wth/1e6,'r',t,exp0d.wth ./ 1e6,'b',t,wped./1e6,'m',t,wpedmax/1e6,'c');
   legend('W_{th}','W_{th, exp}','W_{ped}','W_{ped, max}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MJ')
   subplot(k,1,3);	
   plot(t,zs.wbp/1e6,'r',t,exp0d.wbp ./ 1e6,'b');
   legend('W_{Bpol}','W_{Bpol, exp}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MJ')
   subplot(k,1,4);	
   plot(t,zs.dwdt/1e6,'r',t,exp0d.dwdt ./ 1e6,'b');
   legend('dWdt_{tot}','dWdt_{tot, exp}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MW')
   subplot(k,1,5);	
   plot(t,zs.dwthdt/1e6,'r',t,exp0d.dwthdt ./ 1e6,'b');
   legend('dWdt_{th}','dWdt_{th, exp}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MW')
   subplot(k,1,6);	
   plot(t,zs.dwbpdt/1e6,'r',t,exp0d.dwbpdt ./ 1e6,'b');
   legend('dWdt_{Bpol}','dWdt_{Bpol, exp}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MW')
   subplot(k,1,7);	
   plot(t,zs.esup_fus/1e6,'r',t,exp0d.esup_fus ./ 1e6,'b');
   legend('Esupra_{ \alpha}','Esupra_{ \alpha, exp}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MJ')
   subplot(k,1,8);	
   plot(t,zs.esup_icrh/1e6,'r',t,exp0d.esup_icrh ./ 1e6,'b');
   legend('Esupra_{ICRH}','Esupra_{ICRH, exp}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MJ')
   subplot(k,1,9);	
   plot(t,real(zs.esup_nbi)/1e6 + imag(zs.esup_nbi)/1e6,'r',t,exp0d.esup_nbi ./ 1e6,'b');
   legend('Esupra_{NBI}','Esupra_{NBI, exp}','Location','NorthWest','orientation','horizontal')
   legend(gca,'boxoff')
   ylabel('MJ')
   xlabel('time (s)');

end

joint_axes(h,k);
