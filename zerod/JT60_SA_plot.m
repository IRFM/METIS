%  shape_liste={'Iter like CDR','Demo like CDR','Iter like reduce field','Demo like reduce field'};
%  scenario_liste ={'H mode @ fgr = 0.6','H mode @ fgr = 0.85', ...
%               'Hybrid @ fgr = 0.85','Hybrid @ fgr = 0.39', ...
%               'Steady State @ fgr = 0.86','Steady State @ fgr = 0.58'};
%  gaz_liste ={'H','D','DT','He4'};
%  option_phy ={'Johner','Natural'};

%  post.z0dinput.config_jt60.shape = shape_liste{shape};
%  post.z0dinput.config_jt60.scenario = scenario_liste{scenar};
%  post.z0dinput.config_jt60.gaz      =  gaz_liste{gaz};
%  post.z0dinput.config_jt60.op_phy = option_phy{op_phy};
zs = post.zerod;
profli = post.profil0d;
z0dinput = post.z0dinput;

h = findobj(0,'type','figure','tag','neutron');
if isempty(h)
       h=figure('tag','neutron');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

 
plot(zs.temps,zs.ndd_total,'k', ... 
     zs.temps,zs.ndd_th,'g', ...
      zs.temps,zs.ndd_nbi_thp,'b',zs.temps,zs.ndd_nbi_thn,'c', ...
      zs.temps,zs.ndd_nbi_nbi_p,'r',zs.temps,zs.ndd_nbi_nbi_n,'m',zs.temps,zs.ndd_nbi_nbi_pn,'k:')
      
xlabel('time (s)');
ylabel('neutron/s');
title(sprintf('JT60 SA (%s ->  %s  in %s)', ...
post.z0dinput.config_jt60.shape,post.z0dinput.config_jt60.scenario,post.z0dinput.config_jt60.gaz));
legend('Total','Thermal','Posifif : beam-plasma','Negative: beam-plasma','Posifif : beam-beam','Negative: beam-beam','Negative-Positive: beam-beam')
edition2

% affichage type Johner
ind = find(post.zerod.temps > 80,1);
ipn  = post.zerod.ip./max(post.zerod.ip);
[n,v] = hist(ipn,0:0.05:1);
indflat = find(n == max(n),1);
indflat = find(abs(ipn - v(indflat)) <= 0.01);
fprintf('\n\n-------------------------------------------------------\n');
%fprintf('Plasma shape \t, %s \n',post.z0dinput.config_jt60.shape);
%fprintf('Scenario \t, %s \n',post.z0dinput.config_jt60.scenario);
%fprintf('Gas \t, %s \n',post.z0dinput.config_jt60.gaz);
%fprintf('Model \t, %s \n',post.z0dinput.config_jt60.op_phy);
fprintf('Bt (T) \t, %g\n',post.z0dinput.geo.b0(ind));
fprintf('Ip (MA) \t, %g\n',post.zerod.ip(ind) ./ 1e6);
fprintf('k95 \t, %g \n',post.profil0d.kx(ind,end-1));
fprintf('d95 \t, %g \n',post.profil0d.dx(ind,end-1));
fprintf('R (m) \t, %g\n',post.z0dinput.geo.R(ind));
fprintf('A   \t, %g\n',post.z0dinput.geo.R(ind) ./ post.z0dinput.geo.a(ind));
fprintf('a   \t, %g\n',post.z0dinput.geo.a(ind));
fprintf('V (m^3)   \t, %g\n',post.zerod.vp(ind));
%fprintf('alpha_T / beta_T    \t, %g / %g \n',NaN,NaN);
fprintf('q95    \t, %g\n',post.zerod.q95(ind));
fprintf('Flux (Wb) \t, %g \n',- 2 .* pi .* diff(post.profil0d.psi([1,end],end)));
fprintf('Flux flat top (Wb) \t, %g \n',- 2 .* pi .* diff(post.profil0d.psi(indflat([1,end]),end)));
fprintf('Vloop (V)    \t, %g\n',post.zerod.vloop(ind));
fprintf('li (3)    \t, %g\n',post.zerod.li(ind));
fprintf('ne_bar / ne_gr    \t, %g\n',post.zerod.nbar(ind)./post.zerod.negr(ind));
fprintf('ne_bar (1e19 m^-3)    \t, %g\n',post.zerod.nbar(ind)./1e19);
fprintf('Padd (MW)    \t, %g\n',post.zerod.pin(ind)./1e6);
if all(post.z0dinput.cons.picrh < 1e6)
	fprintf('Padd CD (MW)    \t, %g\n',post.zerod.pecrh(ind)./1e6);
else
	fprintf('Padd CD (MW)    \t, %g\n',post.zerod.pnbi_th(ind)./1e6 + post.zerod.pecrh(ind)./1e6);
end
fprintf('rho_ped    \t, %g\n',0.95);
fprintf('Te_ped  (keV)  \t, %g\n',post.profil0d.tep(ind,end - 1) ./ 1e3);
fprintf('<Te>  (keV)  \t, %g\n',post.zerod.tem(ind) ./ 1e3);
fprintf('Te_0  (keV)  \t, %g\n',post.profil0d.tep(ind,1) ./ 1e3);
fprintf('betan  (%%)   \t, %g\n',post.zerod.betan(ind).* 100);
fprintf('betan_th   (%%)   \t, %g\n',post.zerod.betan(ind) ./ post.zerod.w(ind) .* post.zerod.wth(ind) .* 100);
fprintf('P_OH (MW)    \t, %g\n',post.zerod.pohm(ind)./1e6);
fprintf('P_brem (MW)    \t, %g\n',post.zerod.pbrem(ind)./1e6);
fprintf('P_line (MW)    \t, %g\n',post.zerod.prad(ind)./1e6);
%fprintf('P_cyclo (MW)    \t, %g\n',post.zerod.pcyclo(ind)./1e6);
fprintf('tau_E (s)    \t, %g\n',post.zerod.taue(ind));
%fprintf('Ti0 * ni0 * tau_E (1e20 s keV / m^3)    \t, %g\n',post.zerod.taue(ind) .* post.profil0d.tip(ind,1) ./1e3 .*  ...
%          post.profil0d.nip(ind,1) ./ 1e20);

fprintf('f_BS (%%)   \t, %g  - %g \n',post.zerod.iboot(ind) ./ post.zerod.ip(ind) .* 100,post.zerod.iboot(ind) ./ post.zerod.ipar(ind) .* 100);
fprintf('f_cd  (%%)  \t, %g - %g \n',post.zerod.icd(ind) ./ post.zerod.ip(ind) .* 100,post.zerod.icd(ind) ./ post.zerod.ipar(ind) .* 100);
fprintf('f_ini  (%%)  \t, %g - %g \n',post.zerod.ini(ind) ./ post.zerod.ip(ind) .* 100,post.zerod.ini(ind) ./ post.zerod.ipar(ind) .* 100);
fprintf('Sn total (1e16 neutron/s)  \t, %g\n',post.zerod.ndd_total(ind) ./ 1e16);      
%  fprintf('Sn thermal (1e16 neutron/s)  \t, %g\n',post.zerod.ndd_th(ind) ./ 1e16);      
%  fprintf('Sn beam-plasma p (1e16 neutron/s)  \t, %g\n',post.zerod.ndd_nbi_thp(ind) ./ 1e16);      
%  fprintf('Sn beam-plasma n(1e16 neutron/s)  \t, %g\n',post.zerod.ndd_nbi_thn(ind) ./ 1e16);      
%  fprintf('Sn beam-beam p (1e16 neutron/s)  \t, %g\n',post.zerod.ndd_nbi_nbi_p(ind) ./ 1e16);      
%  fprintf('Sn beam-beam n (1e16 neutron/s)  \t, %g\n',post.zerod.ndd_nbi_nbi_n(ind) ./ 1e16);      
%  fprintf('Sn beam-beam p/n (1e16 neutron/s)  \t, %g\n',post.zerod.ndd_nbi_nbi_pn(ind) ./ 1e16);      

fprintf('grad_P_ped  (1e5 Pa/m)  \t, %g\n',post.zerod.pped(ind) ./ ...
         (post.zerod.rm(ind) .* 0.05) ./ 1e5 .* post.profil0d.grho(ind,end));



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
plot(t,zs.inbicd./1e6,t,zs.ieccd./1e6,t,zs.iboot./1e6,t,zs.iohm./1e6,t,zs.ip./1e6);
legend('N-NBI','ECCD','Boot','Ohm','Ip');
ylabel('MA')
subplot(4,1,1);	
if all(post.z0dinput.cons.picrh < 1e6)
	plot(t,zs.pnbi./1e6,t,zs.pecrh./1e6,t,zs.pohm./1e6,t,zs.prad./1e6 + zs.pbrem./1e6+ zs.pcyclo./1e6);
	legend('P-NBI','ECRH','Ohm','Rad');
else
	plot(t,zs.picrh./1e6,t,zs.pnbi./1e6,t,zs.pecrh./1e6,t,zs.pohm./1e6,t,zs.prad./1e6 + zs.pbrem./1e6+ zs.pcyclo./1e6);
	legend('P-NBI','N-NBI','ECRH','Ohm','Rad');
end
ylabel('MW');
title(sprintf('JT60 SA (%s ->  %s  in %s)', ...
post.z0dinput.config_jt60.shape,post.z0dinput.config_jt60.scenario,post.z0dinput.config_jt60.gaz));
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
z0dmovie(post,profli.temps(ind));
title(sprintf('JT60 SA (%s ->  %s  in %s)', ...
post.z0dinput.config_jt60.shape,post.z0dinput.config_jt60.scenario,post.z0dinput.config_jt60.gaz));


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

plot(xli,profli.tep(ind,:)./1e3,'r',xli,profli.tip(ind,:)./1e3,'m');
hold on
plot(xli,profli.nep(ind,:)./1e19,'b',xli,profli.nD(ind,:)./1e19,'c');
xlabel('x (su)');
title(sprintf('JT60 SA (%s ->  %s  in %s)', ...
post.z0dinput.config_jt60.shape,post.z0dinput.config_jt60.scenario,post.z0dinput.config_jt60.gaz));
legend('T_e (keV)','T_i (keV)','n_e (10^{19} m^{-3})','n_D (10^{19} m^{-3})');
edition2



% les sources de neutrons
h = findobj(0,'type','figure','tag','neutron2');
if isempty(h)
       h=figure('tag','neutron2');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1,'color',[1 1 1])

semilogy(zs.temps,zs.ndd,'r')
      
xlabel('time (s)');
ylabel('neutrons/s');
title(sprintf('JT60 SA (%s ->  %s  in %s)', ...
post.z0dinput.config_jt60.shape,post.z0dinput.config_jt60.scenario,post.z0dinput.config_jt60.gaz));
set(gca,'ylim',[1e12,1e18])     
edition2

fprintf('DD total  = %g\n',zs.ndd(ind));


vs = trapz(zs.temps,zs.vloop)
V = zs.vloop(end-2)
neutron = trapz(zs.temps,zs.ndd)
fn = zs.ndd(end-2)





