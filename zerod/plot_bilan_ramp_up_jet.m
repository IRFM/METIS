% plot bilan ramp-up JET
% time intervalle
tmin = min(post.zerod.temps);
switch post.z0dinput.shot
case 73221
  tmax = 44.2;
case 73224
  tmax = 44.5;
case 78834
  tmax = 46.6;
case {78842,738842}
  post.z0dinput.shot = 78842;
  tmax = 46.6;
case 79649
  tmax = 44.5;
case 89723
  tmax = 53.4;
otherwise
  tmax = max(post.zerod.temps);
end
indt = find((post.zerod.temps >= tmin) & (post.zerod.temps <= tmax));

h = findobj(0,'type','figure','tag','bilan_rampup_jet');
if isempty(h)
       h=figure('tag','bilan_rampup_jet');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

subplot(3,2,1)
plot(post.zerod.temps,post.zerod.ip ./ 1e6, ...
     post.z0dinput.cons.temps,post.z0dinput.geo.b0, ...
     post.z0dinput.cons.temps,post.z0dinput.cons.nbar ./ 1e19, ...
     post.zerod.temps,post.zerod.zeff,post.zerod.temps,post.zerod.pin ./ 1e7);
legend({'I_{p} MA','B_{t,0} T','n_{bar} 10^{19} m^{-3}','Z_{eff}','P_{in} (10 MW)'},'location','best')
set(gca,'ylim',[0,Inf]);    
set(gca,'xlim',[tmin,tmax]);
title(sprintf('shot JET #%d',post.z0dinput.shot));

subplot(3,2,2)
plot(post.zerod.temps,post.zerod.vloop,'b',post.zerod.temps,post.z0dinput.exp0d.vloop,'r');
legend({'V_{loop} METIS','V_{loop} JET'},'location','best')
ylabel('V');
set(gca,'ylim',[floor(min(post.zerod.vloop(indt))),ceil(max(post.zerod.vloop(indt)))]);    
set(gca,'xlim',[tmin,tmax]);

subplot(3,2,3)
plot(post.zerod.temps,post.zerod.li,'r',post.zerod.temps,post.z0dinput.exp0d.li,'r-.', ...
     post.zerod.temps,post.zerod.betap,'b',post.zerod.temps,post.z0dinput.exp0d.betap,'b-.', ...
     post.zerod.temps,post.zerod.betap + post.zerod.li /2,'k',post.zerod.temps,post.z0dinput.exp0d.betap + post.z0dinput.exp0d.li / 2,'k-.');
legend({'l_{i} METIS','l_{i} EFIT','\beta_{p} METIS','\beta_{p} EFIT','\beta_{p}+l_{i}/2 METIS','\beta_{p}+l_{i}/2 EFIT'},'location','best') 
set(gca,'ylim',[0,ceil(max(post.zerod.li(indt)))]);    
  set(gca,'xlim',[tmin,tmax]);
  
subplot(3,2,4)
plot(post.zerod.temps,medfilt1(post.zerod.te0,3) ./ 1e3,'r',post.zerod.temps,post.z0dinput.exp0d.te0 ./ 1e3,'r-.', ...
     post.zerod.temps,post.zerod.tem ./ 1e3,'b',post.zerod.temps,post.z0dinput.exp0d.tem ./ 1e3,'b-.');
ylabel('keV');
xlabel('time (s)');
legend({'T_{e,0} METIS','T_{e,0} JET','<T_e> METIS','<T_e> JET'},'location','best');
set(gca,'xlim',[tmin,tmax]);
set(gca,'ylim',[0,ceil(max(medfilt1(post.zerod.te0(indt),3) ./ 1e3))]);    
   
subplot(3,2,5)
plot(post.zerod.temps,post.zerod.ne0 ./ 1e19,'r',post.zerod.temps,post.z0dinput.exp0d.ne0 ./ 1e19,'r-.', ...
     post.zerod.temps,post.zerod.nem ./ 1e19,'b',post.zerod.temps,post.z0dinput.exp0d.nem ./ 1e19,'b-.');
ylabel('10^{19} m^{-3}');
xlabel('time (s)');
legend({'n_{e,0} METIS','n_{e,0} JET','<n_e> METIS','<n_e> JET'},'location','best');
set(gca,'ylim',[0,ceil(max(post.zerod.ne0(indt) ./ 1e19))]);    
set(gca,'xlim',[tmin,tmax]);
   
subplot(3,2,6)
plot(post.zerod.temps,post.zerod.prad/1e6,'b',post.zerod.temps,post.z0dinput.exp0d.prad/1e6,'b-.',...
     post.zerod.temps,post.zerod.w/1e6,'r',post.zerod.temps,post.z0dinput.exp0d.w/1e6,'r-.');
xlabel('time (s)');
legend({'P_{rad} METIS (MW)','P_{rad} JET (MW)','W_{MHD} METIS (MJ)','W_{MHD} JET (MJ)'},'location','best');
set(gca,'ylim',[0,ceil(max(post.zerod.prad(indt)./1e6))]);    
set(gca,'xlim',[tmin,tmax]);

edition2
