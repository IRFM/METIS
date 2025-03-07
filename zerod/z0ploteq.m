% script du  plot de l'equilibre 0d
zs   = post.zerod;
cons = post.z0dinput.cons;
geo  = post.z0dinput.geo;
exp0d  = post.z0dinput.exp0d;
t    = zs.temps;
% precalcul
if post.z0dinput.mode_exp == 0
	spr  = data.equi.spr;
	vpr  = data.equi.vpr;
        jmoy = data.equi.jmoy;
	pmax  =data.equi.rhomax .* trapz(param.gene.x,vpr .* data.coef.eta .*jmoy .^ 2,2);
        imax  = data.equi.rhomax .* trapz(param.gene.x,spr .* jmoy,2);
	res  =  pmax./ imax .^ 2;
	tauip_cronos  =  2 .* data.gene.wbp ./ res ./ data.gene.ip .^ 2;
end
tauip = zs.tauip;
vloopmax  = ermax(zs.tem .* (1+zs.ate),zs.nem .* (1+ zs.ane),zs.zeff,1) .* geo.R .* 2 .* pi;

% maximum beta pour la stabilite verticale (Bz) : regle heuristique
beta_max = 0.5 .* geo.R ./ geo.a;
beta_GS = geo.R ./ geo.a + 1 - zs.li ./ 2;

% self interne a partir du flux
mu0  = 4 .* pi .* 1e-7;
tps  = post.profil0d.temps;
lips = (post.profil0d.psi(:,1) - post.profil0d.psi(:,end)) ./ interp1(t,zs.ip,tps,'pchip',Inf) ./ ...
       mu0 .* post.profil0d.Raxe(:,1);

h = findobj(0,'type','figure','tag','z0ploteq');
if isempty(h)
       h=figure('tag','z0ploteq');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

k    = 9;
if post.z0dinput.mode_exp == 0
   subplot(k,1,1);	
   plot(t,post.z0dinput.geo.b0,'r',data.gene.temps,data.geo.b0 );
   ylabel('B0 (T)')
   title(sprintf('Zerod : %s@%d/equilibrium ( r -> 0D, b: -> cronos)', ...
          post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,zs.ip./1e6,'r',data.gene.temps,data.gene.ip./1e6,':b');
   ylabel('Ip (MA)');
   subplot(k,1,3);	
   plot(t,zs.betap,'r',t,zs.betaptot,':r',data.gene.temps,data.gene.betap,':b',t,beta_max,'g-.',t,beta_GS,'g:');
   ylabel('Betap (thermal -r and total r:)');
   subplot(k,1,4);	
   plot(t,zs.li,'r',data.gene.temps,data.gene.li,':b');
   ylabel(' li ');
   subplot(k,1,5);	
   plot(t,zs.qa,'r',t,zs.qeff,'m',t,zs.q95,'-.r',data.gene.temps,data.prof.q(:,end),':b', ...
	t,ones(size(t)),':g',t,3 .* ones(size(t)),'-g',t,2 .* ones(size(t)),'-.g');
   set(gca,'ylim',[0,min(20,max(zs.qa))])
   ylabel('q_a, q_e_f_f & q_9_5');
   subplot(k,1,6);	
   plot(t,zs.q0,'r',t,zs.qmin,'r-.',data.gene.temps,data.prof.q(:,1),':b',data.gene.temps,min(data.prof.q,[],2),'-.b', ...
	t,ones(size(t)),'-g',t,2 .* ones(size(t)),':g',t,3/2 .* ones(size(t)),'-.g');
   set(gca,'ylim',[0,min(5,max(zs.q0))])
   ylabel('q_0 & q_{min}');
   subplot(k,1,7);	
   asser = zs.asser;
   asser(asser==0) = NaN;
   plot(t,zs.vloop,'r',data.gene.temps,data.gene.vres,':b',t,0.*t,'g',t,vloopmax,'g-.',t,-vloopmax,'g-.',t,asser .* zs.vloop,'om');
   ylabel('Vloop (V)');
   set(gca,'ylim',sort([max(-0.1,min(zs.vloop-eps)),min(3,max(zs.vloop+eps))]))
   subplot(k,1,8);	
   semilogy(t,zs.RR.*1e6,'r',data.gene.temps,res.*1e6,':b');
   ylabel('Res (1e-6 ohm)');
   z0loglin(gca);
   subplot(k,1,9);	
   plot(t,tauip,'r',data.gene.temps,tauip_cronos,':b');
   ylabel('tau_i_p (s) ');
   xlabel('time (s)');
else
   subplot(k,1,1);	
   plot(t,post.z0dinput.geo.b0,'r');
   ylabel('B0 (T)')
   title(sprintf('Zerod : %s@%d/equilibrium ( r -> 0D, b: -> experiment)', ...
          post.z0dinput.machine,post.z0dinput.shot));
   subplot(k,1,2);	
   plot(t,zs.ip./1e6,'r',t,exp0d.ip./1e6,':b');
   ylabel('Ip (MA)');
   subplot(k,1,3);	
   plot(t,zs.betap,'r',t,zs.betaptot,'r:',t,exp0d.betap,':b',t,beta_max,'g-.',t,beta_GS,'g:');
   ylabel('Betap (thermal -r and total :r)');
   subplot(k,1,4);	
   plot(t,zs.li,'r',t,exp0d.li,':b');
   ylabel(' li ');
   subplot(k,1,5);	
   plot(t,zs.qa,'r',t,zs.qeff,'m',t,zs.q95,'-.r',t,exp0d.qa,':b', ...
	t,ones(size(t)),'g',t,3 .* ones(size(t)),'-.g',t,2 .* ones(size(t)),':g');
   set(gca,'ylim',[0,min(20,max(zs.qa))])
   ylabel('q_a , q_e_f_f & q_9_5');
   subplot(k,1,6);	
   plot(t,zs.q0,'r',t,zs.qmin,'r-.',t,exp0d.q0,':b', ...
	t,ones(size(t)),'g',t,2 .* ones(size(t)),'-.g',t,3/2 .* ones(size(t)),'-.g');
   set(gca,'ylim',[0,min(5,max(zs.q0))])
   ylabel('q_0 & q_{min}');
   subplot(k,1,7);	
   asser = zs.asser;
   asser(asser==0) = NaN;
   plot(t,zs.vloop,'r',t,exp0d.vloop,':b',t,0.*t,'g',t,vloopmax,'g-.',t,-vloopmax,'g-.',t,asser .*  zs.vloop,'om');
   ylabel('Vloop (V)');
   set(gca,'ylim',sort([max(-0.1,min(zs.vloop-eps)),min(3,max(zs.vloop+eps))]))
   subplot(k,1,8);	
   semilogy(t,zs.RR.*1e6,'r');
   z0loglin(gca);
   ylabel('Res (1e-6 ohm)');
   subplot(k,1,9);	
   plot(t,tauip,'r');
   ylabel('tau_i_p (s)');
   xlabel('time (s)');

end
joint_axes(h,k);
