% script pour le plot des neutrons
scenario ='';
while isempty(scenario) 
	scenario = upper(input('Scenario ? H1, H2, HYB1 (HYB11) , HYB2 (HYB12), SS1, SS2 or SS3 ?','s'));
end
load(sprintf('METIS_FT3_%s',scenario));
zs = post.zerod;
profli = post.profil0d;
z0dinput = post.z0dinput;

% 1 - plot du scenario
h = findobj(0,'type','figure','tag','z0dsc');
if isempty(h)
       h=figure('tag','z0dsc');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

%zs   = post.zerod;
cons = post.z0dinput.cons;
exp0d= post.z0dinput.exp0d;
t    = zs.temps;
subplot(4,1,3);	
plot(t,zs.ilh./1e6,t,zs.ieccd./1e6,t,zs.iboot./1e6,t,zs.iohm./1e6,t,zs.ip./1e6);
legend('LH','ECCD','Boot','Ohm','Ip');
ylabel('MA')
subplot(4,1,1);	
plot(t,zs.picrh./1e6,t,zs.plh./1e6,t,zs.pecrh./1e6,t,zs.pohm./1e6,t,zs.prad./1e6 + zs.pbrem./1e6+ zs.pcyclo./1e6);
legend('ICRH','LH','ECRH','Ohm','Rad');
ylabel('MW');
title(sprintf('Metis : FT3 (%s)',scenario));
subplot(4,1,4)
plot(t,zs.vloop,'b',t,zs.li,'r',t,zs.betap,'g',t,zs.betaptot,'g:');
legend('Vloop','li','betap (th)','betap (tot)');
set(gca,'ylim',[-0.1,3]);
ylabel('V, su ,su')
xlabel('time (s)');
subplot(4,1,2);	
plot(t,zs.ne0./1e19,t,zs.nem./1e19,t,zs.nDm./1e19);
legend('ne_0','<n_e>','<n_D>');
ylabel('10^1^9 m ^-^3');
edition2

% 2 D plot 
z0dmovie(post,profli.temps(end-2));
title(sprintf('Metis : FT3 (%s)',scenario));


% quelques profils (ne,nD Te,Ti)
profli = z0neutron_profil(zs,post.profil0d);
xli    = profli.xli;
h = findobj(0,'type','figure','tag','profil0d');
if isempty(h)
       h=figure('tag','profil0d');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])
% visulaisation 3D des profils 0D

k = 3;
subplot(k,k-1,1)
zplotprof(gca,t,xli,profli.tep./1e3,'color','r');
zplotprof(gca,t,xli,profli.tip./1e3,'color','b');
title('T_e (r) & T_i (b) [keV]')
subplot(k,k-1,2)
zplotprof(gca,t,xli,profli.nep./1e19,'color','r');
zplotprof(gca,t,xli,profli.nD./1e19,'color','b');
title('n_e (r), n_D (b) [10^1^9 m^-^3]')
subplot(k,k-1,3)
zplotprof(gca,t,xli,profli.qjli,'color','b');
s  = pdederive(xli,profli.qjli,0,2,2,1) .* (ones(size(t)) * xli) ./ profli.qjli;
zplotprof(gca,t,xli,s,'color','r');
hold on
hhp = plot([0,1],[1,1],'g',[0,1],[1.5,1.5],'g-.',[0,1],[2,2],'g:');
set(hhp,'linewidth',0.5);
title('q (b) & s (r)')
subplot(k,k-1,4)
zplotprof(gca,t,xli,profli.jboot./1e6,'color','b');
zplotprof(gca,t,xli,profli.jlh./1e6,'color','r');
zplotprof(gca,t,xli,profli.jeccd./1e6,'color','m');
title('J_N_I (b:boot, r:LHCD, m: ECCD ) [MA m^-^2]')
subplot(k,k-1,5)
zplotprof(gca,t,xli,profli.sn0_th./1e16,'color','r');
title('neutron dd ')
xlabel('x (su)');
subplot(k,k-1,6)
zplotprof(gca,t,xli,profli.picrh./1e6,'color','r');
zplotprof(gca,t,xli,profli.plh./1e6,'color','m');
zplotprof(gca,t,xli,profli.pecrh./1e6,'color','c');
zplotprof(gca,t,xli,profli.pohm./1e6,'color','b');
zplotprof(gca,t,xli,profli.prad./1e6 + profli.pbrem./1e6 + profli.pcyclo./1e6,'color','k');
title('P_I_C_R_H (r),P_L_H (m), P_E_C_R_H (c), Pohm (b), Prad (k) [MW m^-^3]')
xlabel('x (su)');
edition2

% donnee de profils 
h = findobj(0,'type','figure','tag','profilsimple');
if isempty(h)
       h=figure('tag','profilsimple');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

plot(xli,profli.tep(end-2,:)./1e3,'r',xli,profli.tip(end-2,:)./1e3,'m');
hold on
plot(xli,profli.nep(end-2,:)./1e19,'b',xli,profli.nD(end-2,:)./1e19,'c');
xlabel('x (su)');
title(sprintf('Metis : FT3 (%s)',scenario));
legend('T_e (keV)','T_i (keV)','n_e (10^{19} m^{-3})','n_D (10^{19} m^{-3})');
edition2



% les sources de neutrons
h = findobj(0,'type','figure','tag','neutron');
if isempty(h)
       h=figure('tag','neutron');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

semilogy(zs.temps,zs.ndd,'r')
      
xlabel('time (s)');
ylabel('neutrons/s');
title(sprintf('FT3 (%s)',scenario));
set(gca,'ylim',[1e12,1e18])     
edition2

fprintf('DD total  = %g\n',zs.ndd(end-2));


vs = trapz(zs.temps,zs.vloop)
V = zs.vloop(end-2)
neutron = trapz(zs.temps,zs.ndd)
fn = zs.ndd(end-2)

