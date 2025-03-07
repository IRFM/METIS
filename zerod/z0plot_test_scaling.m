% script du  plot de test de metis pour le respect des lois de confinement
% constante physique (phys)
phys.c           =   2.99792458e8;             % speed of light in vacuum (m/s)  (definition)
phys.h           =   6.62606876e-34;           % Planck constant (J*s) (+/- 0.0000052e-34)
phys.e           =   1.602176462e-19;          % electron charge (C)   (+/- 0.000000063e-19)
phys.mu0         =   4*pi*1e-7;                % permeablity of vacuum (H/m) (definition)
phys.epsi0       =   1./phys.c.^2./phys.mu0;   % permitivity of vacuum (F/m)  (definition)
phys.g           =   6.673e-11;                % gravitation constant (N*m^2/kg^2) (+/- 0.010e-11)
phys.k           =   1.3806503e-23;            % Boltzmann constant (J/K)  (+/- 0.0000024e-23)
phys.alpha       =   7.297352533e-3 ;          % fine structure constant (+/- 0.000000027e-3 )
phys.me          =   9.10938188e-31;           % electron mass (at rest) (kg) (+/- 0.00000079e-31)
phys.mp          =   1.6726485e-27;            % proton mass (at rest) (kg)
phys.ua          =   1.66053873e-27;           % Atomic mass unit (kg) (+/- 1.00000013e-27)
phys.avo         =   6.02214199e23;            % Avogadro number (mol^-1) (+/- 0.00000047e23)
phys.sigma       =   5.670400e-8;              % Stephan constant ( W*m^-2*K^-4) (+/- 0.000040e-8)



% variables pour simplifier l'ecriture
zs   = post.zerod;
cons = post.z0dinput.cons;
geo  = post.z0dinput.geo;
exp0d  = post.z0dinput.exp0d;
t    = zs.temps;
disrup = double(zs.disrup);
disrup(disrup == 0) = NaN;
disrup(isfinite(disrup)) = 0;
% rapport pour le facteur d'amplification (Ealpha +En)/Ealpha
rfan = (3.56e6 + 14.03e6) ./ 3.56e6 ;
padd = (cons.picrh + cons.pecrh + real(cons.pnbi) + imag(cons.pnbi) + zs.pohm + cons.plh) / 1e6;
pin  = padd + zs.pfus ./ 1e6;
ploss = max(1e-6,pin - (zs.pbrem + zs.pcyclo + 0.33 .* zs.prad + zs.pioniz) ./ 1e6);
ip = zs.ip ./ 1e6;
Bt = geo.b0;
ne = real(cons.nbar)./ 1e19;
R  =  geo.R;
a  =  geo.a;
ep  = a ./ R;
K   = geo.K; 
Vp = zs.vp;
Ka  = Vp ./ (2*pi^2.*R.*a.^2);
meff = zs.meff;
% precalcul

h = findobj(0,'type','figure','tag','z0plottsc');
if isempty(h)
       h=figure('tag','z0plottsc');
else
       figure(h);
end   
clf
set(h,'defaultaxesfontsize',12,'defaultaxesfontweight','bold','defaultaxesfontname','times', ...
	'defaultlinelinewidth',1.5,'color',[1 1 1])

k    = 4;
% respect des lois interne (sans modification)
subplot(k,1,1);
plot(t,zs.taue,'r',t,zs.modeh .* zs.tauh + (~zs.modeh .* zs.tauthl),'--b')	
legend('tau_E','tau_{internal scaling }','Location','north','orientation','horizontal')
%legend(gca,'boxoff')
ylabel('tau_E(s)')
title(sprintf('Zerod : %s@%d/test METIS scaling ', ...
          post.z0dinput.machine,post.z0dinput.shot));

% test de la coherence entre la pression, le contenu en energie et le temps de confinement
pression = phys.e .* (post.profil0d.nep .* post.profil0d.tep + post.profil0d.nip .* post.profil0d.tip);
wtest = 3/2 * trapz(post.profil0d.xli,post.profil0d.vpr .* pression,2);
subplot(k,1,2);	
plot(t,zs.wth./1e6,'r',post.profil0d.temps,wtest./1e6,'--b',t,zs.taue .* zs.ploss./1e6,'-.c',t,zs.taue .* zs.pth./1e6,':m');
ylabel('W_{th} (MJ)')
legend('METIS','inte(Pression)','tau_E * P_{loss}','tau_E * P_{th}','Location','north','orientation','horizontal')
%legend(gca,'boxoff')

% tests des données des lois d'echelle
subplot(k,1,3);	
% loi standard ITER
% ITERH-96P(th)        
tauthl  = 23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* pin .^ -0.73 .* ...
      R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2; % s       

% ITERH-98P(y,2)        
tauh   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* ploss .^ -0.69 .* ...
	R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s    

plot(t,zs.taue,'r',t,zs.modeh .* tauh + (~zs.modeh .* tauthl),'--b')	
legend('tau_E','tau_{scaling ITERH-96P(th) or ITERH-98P(y,2)}','Location','north','orientation','horizontal')
%legend(gca,'boxoff')
ylabel('tau_E(s)')
set(gca,'ylim',[0,max(zs.taue)]);

subplot(k,1,4);	
%  % loi standard ITER
%  % ITERH-96P(th)        
%  tauthl  = 23e-3  .* ip .^ 0.96 .* Bt .^ 0.03 .* ne .^ 0.4 .* (zs.pin/1e6) .^ -0.73 .* ...
%        R .^ 1.83 .* K .^ 0.64 .* ep .^ -0.06 .* meff .^ 0.2; % s       
%  
%  % ITERH-98P(y,2)        
%  tauh   = 56.2e-3  .* ip .^ 0.93 .* Bt .^ 0.15 .* ne .^ 0.41 .* (zs.ploss/1e6) .^ -0.69 .* ...
%  	R .^ 1.97 .* Ka .^ 0.78 .* ep .^ 0.58 .* meff .^ 0.19;    % s    

indic = zs.modeh;
indic(indic < 1) = NaN;
indic(indic == 1) = 0;
plot(t,zs.taue ./  tauh ,'r',t,zs.taue ./ tauthl,'--b',t,zs.wth ./  (tauh .* ploss) / 1e6,'-.m',t,zs.wth ./ (tauthl .* pin) / 1e6,':c',t,indic,'^g')
hold on
plot(t([1,end]),[1,1],'color','k','linewidth',0.5,'linestyle',':')	
legend('H_H (ITERH-98P(y,2), \tau_E/\tau_{sc})','H_L (ITERH-96P(th), \tau_E/\tau_{sc})','H_H (ITERH-98P(y,2),W_{th}/W_{sc})','H_L (ITERH-96P(th),W_{th}/W_{sc})','H mode flag','Location','north','orientation','horizontal')
%legend(gca,'boxoff')
ylabel('H factor')
xlabel('time (s)');
set(gca,'ylim',[0,4]);

joint_axes(h,k);
