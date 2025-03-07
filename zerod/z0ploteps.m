zs   = post.zerod;
cons = post.z0dinput.cons;
exp0d  = post.z0dinput.exp0d;
t    = zs.temps;
disrup = double(zs.disrup);
disrup(disrup == 0) = NaN;
disrup(disrup == 1) = 0;
% rapport pour le facteur d'amplification (Ealpha +En)/Ealpha
rfan = (3.56e6 + 14.03e6) ./ 3.56e6 ;
padd = cons.picrh + cons.pecrh + cons.pnbi + zs.pohm + zs.plh;


h = findobj(0,'type','figure','tag','z0ploteps');
if isempty(h)
       h=figure('tag','z0ploteps');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

% decodage de simout
temps_sim          = simout.signals.values(:,9);
vt                 = ones(size(temps_sim));
iboot_sim          = simout.signals.values(:,1);
frboot_sim         = simout.signals.values(:,2);
picrh_fus_sim      = simout.signals.values(:,3);
pecrh_fus_sim      = simout.signals.values(:,4);
Q_sim              = simout.signals.values(:,5);
plh_control_sim    = simout.signals.values(:,6);
picrh_control_sim  = simout.signals.values(:,7);
pecrh_control_sim  = simout.signals.values(:,8);

subplot(5,1,1)
plot(post.zerod.temps(3:end),post.zerod.ip(3:end)./1e6,post.zerod.temps(3:end),post.zerod.iboot(3:end)./1e6, ...
     post.zerod.temps(3:end),(post.zerod.ilh(3:end) - iboot_sim + post.zerod.iboot(3:end)) ./1e6, temps_sim,iboot_sim./1e6);
legend('I_P','I_{boot,real}','I_{LH,control}','I_{boot,simulation}');
set(gca,'ylim',[0 ,1]);
set(gca,'xlim',[0 ,60]);
ylabel('MA');

subplot(5,1,2)
plot(temps_sim,Q_sim,'b',temps_sim,vt * 15,'r', ...
     post.zerod.temps(3:end),post.zerod.ndd(3:end)./1e12,'g');
legend('Q_{simulation}','Q_{target}','R_{DD,real} (10^{12} s^{-1})');
set(gca,'ylim',[0 ,20]);
set(gca,'xlim',[0 ,60]);

subplot(5,1,3)
plot(temps_sim,picrh_control_sim./1e6,post.zerod.temps,post.zerod.picrh_th./1e6,temps_sim,vt*12,'r');
legend('P_{ICRH,control}','P_{ICRH,total}','P_{ICRH,limit}')
ylabel('MW');
set(gca,'xlim',[0 ,60]);
set(gca,'Ylim',[0 ,15]);

subplot(5,1,4)
plot(temps_sim,pecrh_control_sim./1e6,post.zerod.temps,post.zerod.pecrh./1e6,temps_sim,vt*6,'r');
legend('P_{ECRH,offset}','P_{ECRH,total}','P_{ECRH,limit}')
ylabel('MW');
set(gca,'xlim',[0 ,60]);
set(gca,'Ylim',[0 ,8]);

subplot(5,1,5)
plh_control_sim = min(plh_control_sim,interp1(post.zerod.temps,post.zerod.plh_th,temps_sim,'pchip','extrap'));
plot(temps_sim,plh_control_sim./1e6,post.zerod.temps,post.zerod.plh_th./1e6,temps_sim,vt*8,'r');
legend('P_{LH,control}','P_{LH,total}','P_{LH,limit}')
ylabel('MW');
xlabel('time (s)');
set(gca,'xlim',[0 ,60]);
set(gca,'Ylim',[0 ,10]);

edition2

pfusion_sim = picrh_fus_sim + pecrh_fus_sim;
li_sim = post.zerod.li(3:end);
nbar_sim = post.zerod.nbar(3:end);
pecrh_sim = post.z0dinput.cons.pecrh(3:end);
picrh_sim = post.z0dinput.cons.picrh(3:end);
plh_sim = post.z0dinput.cons.plh(3:end);


save data4plot temps_sim iboot_sim frboot_sim Q_sim pfusion_sim li_sim nbar_sim  pecrh_sim picrh_sim plh_sim -V6